import os

import pandas as pd
import requests
from Bio import Entrez
from io import BytesIO
from tqdm import tqdm

from .paper_retriever import PaperRetriever
from .image_extractor import ImageExtractor
from .reference_retriever import ReferenceRetriever

tqdm.pandas()

class PubMedSearcher:
    """
    A class to search PubMed for articles, retrieve metadata, download full texts, and process references.

    Attributes:
    - search_string (str): Query string for searching PubMed.
    - df (DataFrame): Stores retrieved article metadata, references, and processing statuses.
    - email (str): Email address required for querying PubMed via Entrez.

    Core Methods:
    - search(count, min_date, max_date, order_by, only_open_access, only_case_reports): 
      Retrieves PubMed articles based on the search query with filtering options.
    - fetch_abstracts(): 
      Fetches missing abstracts from PubMed for articles in the DataFrame.
    - download_articles(allow_scihub, download_directory, max_articles): 
      Downloads full-text PDFs, prioritizing open access sources.
    - extract_images(image_directory): 
      Extracts images from downloaded PDFs and stores them.
    - fetch_references(): 
      Retrieves references for each article using multiple sources (PubMed, PMC, CrossRef, Europe PMC).
    - fetch_cited_by(): 
      Retrieves articles that cite each article using Europe PMC.
    - download_xml_fulltext(download_directory): 
      Downloads XML full-text versions when available (PubMed Open Access, Europe PMC).
    - save(csv_path): 
      Saves the DataFrame to a CSV file.

    Internal Helper Methods:
    - _validate_dataframe(df): Ensures DataFrame contains necessary columns.
    - _parse_records_to_df(records_xml_bytes): Converts retrieved PubMed records to a DataFrame.
    - _get_xml_for_pmcid(pmcid): Fetches full-text XML for an article using a PMCID.
    
    This class integrates PubMed's Entrez API, Europe PMC, and CrossRef to facilitate systematic literature retrieval, 
    citation analysis, and document processing.
    """

    def __init__(self, search_string=None, df=None, email=""):
        """
        Initializes the PubMedSearcher object with a search string and an optional DataFrame.

        Parameters:
        search_string (str): The search string to use when querying PubMed.
        df (DataFrame): An optional DataFrame that may already contain articles, previous search results, etc.
        email (str): The email address to use when querying PubMed.
        """
        self.search_string = search_string
        self.df = df if df is not None else pd.DataFrame()
        if df is not None:
            self._validate_dataframe(df)
        elif search_string is not None:
            self.df = pd.DataFrame()
        if email:
            self.email = email
        else:
            print("Please provide an email address to use for querying PubMed.")
            raise ValueError("Email address is required for PubMed queries.")

    def search(self, count=10, min_date=None, max_date=None, order_by='chronological', only_open_access=False, only_case_reports=False):
        """
        Searches PubMed for articles based on the search string and retrieves the specified number of articles.

        Parameters:
        count (int): The number of articles to retrieve.
        min_date (int, optional): The minimum publication year to consider.
        max_date (int, optional): The maximum publication year to consider.
        order_by (str, optional): The order in which to retrieve articles. Can be 'chronological' or 'relevance'. Defaults to 'chronological'.
        only_open_access (bool, optional): Whether to retrieve only open access articles. Defaults to False.
        only_case_reports (bool, optional): Whether to retrieve only case reports. Defaults to False.
        """
        if not self.search_string:
            raise ValueError("Search string is not provided")
        
        additional_filters = []
        if only_open_access:
            # Combine open access filters with an OR condition
            additional_filters.append("(open access[filter] OR free full text[sb])")
        
        if only_case_reports:
            # Add filter for case reports
            additional_filters.append("case reports[pt]")
        
        # Join any additional filters to the search string with AND
        if additional_filters:
            self.search_string += " AND " + " AND ".join(additional_filters)

        Entrez.email = self.email
        search_params = {
            'db': "pubmed",
            'term': self.search_string,
            'retmax': count,
            'sort': 'relevance' if order_by == 'relevance' else 'pub date',
        }

        if min_date is not None:
            search_params['mindate'] = str(min_date)
        if max_date is not None:
            search_params['maxdate'] = str(max_date)

        search_handle = Entrez.esearch(**search_params)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        id_list = search_results['IdList']
        fetch_handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        records_xml_bytes = fetch_handle.read()
        fetch_handle.close()

        records_df = self._parse_records_to_df(records_xml_bytes)
        self.df = pd.concat([self.df, records_df], ignore_index=True)

    def download_articles(self, allow_scihub=False, download_directory="pdf_downloads", max_articles=None):
        """
        Downloads full-text PDFs for articles in the DataFrame, prioritizing open-access sources.

        Parameters:
        - allow_scihub (bool, optional): Whether to attempt Sci-Hub as a fallback source if open access is unavailable. Defaults to True.
        - download_directory (str, optional): Directory where downloaded PDFs will be stored. Defaults to "pdf_downloads".
        - max_articles (int, optional): Maximum number of articles to download in a single execution. If None, all available articles will be processed.

        Process:
        1. Checks if the DataFrame is populated and contains necessary columns ('pmid' and 'doi').
        2. Initializes tracking columns ('download_complete' and 'pdf_filepath') if they do not exist.
        3. Iterates through articles, skipping those already marked as 'complete' or 'unavailable'.
        4. Attempts to retrieve the full-text PDF using the PaperRetriever class.
        5. Updates the DataFrame:
           - If the PDF is successfully downloaded, stores the file path and marks the download as 'complete'.
           - If unavailable, marks it as 'unavailable'.
        6. Saves progress after each successful download.

        Returns:
        - self: Returns the updated PubMedSearcher instance with the download status updated.

        Notes:
        - Sci-Hub should only be enabled if legally permissible in your jurisdiction.
        - Articles with missing or invalid DOIs are automatically skipped.
        - Uses tqdm for progress tracking.
        """
        if self.df.empty:
            print("DataFrame is empty.")
            return
        if 'pmid' not in self.df.columns or 'doi' not in self.df.columns:
            print(f"DataFrame is missing required columns for article download (pmid and doi)")
            return
        if 'download_complete' not in self.df.columns:
            self.df['download_complete'] = 'not_started'
        if 'pdf_filepath' not in self.df.columns:
            self.df['pdf_filepath'] = None
        articles_processed = 0
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0], desc="Downloading articles"):
            if row.get('download_complete') == 'complete' or row.get('download_complete') == 'unavailable':
                continue
            if max_articles and articles_processed >= max_articles:
                break
            pmid = row.get('pmid')
            doi = row.get('doi')
            if not doi or len(str(doi)) < 5:
                self.df.at[index, 'download_complete'] = 'unavailable'
                continue
            pdf_filepath = PaperRetriever(pmid=pmid,
                                        doi=doi,
                                        email=self.email,
                                        allow_scihub=allow_scihub,
                                        download_directory=download_directory
                                    ).download().filepath
            if pdf_filepath in [None, '', 'unavailable']:
                self.df.at[index, 'download_complete'] = 'unavailable'
                self.df.at[index, 'pdf_filepath'] = None
            else:
                self.df.at[index, 'download_complete'] = 'complete'
                self.df.at[index, 'pdf_filepath'] = pdf_filepath
            articles_processed += 1
            self.save()
        return self
    
    def extract_images(self):
        """
        Extracts images from the PDFs of the articles in the DataFrame using the ImageExtractor class.
        Only applies to successfully downloaded PDFs.
        
        Parameters:
        - image_directory (str): The base directory where extracted images will be stored.
        """
        if self.df.empty:
            print("DataFrame is empty. No articles to extract images from.")
            return
        
        required_columns = ['download_complete', 'pdf_filepath']
        for col in required_columns:
            if col not in self.df.columns:
                print(f"Error: DataFrame is missing required column '{col}'.")
                return
        
        # Initialize columns for image paths and extraction status if they don't exist
        if 'image_paths' not in self.df.columns:
            self.df['image_paths'] = [[] for _ in range(len(self.df))]
        if 'image_extraction_complete' not in self.df.columns:
            self.df['image_extraction_complete'] = 'not_started'
        
        # Iterate over the DataFrame with a progress bar
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0], desc="Extracting Images"):
            # Skip if download was not complete
            if row.get('download_complete') != 'complete':
                continue
            
            # Skip if image extraction is already complete
            if row.get('image_extraction_complete') == 'complete':
                continue
            
            pdf_filepath = row.get('pdf_filepath')
            if not pdf_filepath or not os.path.isfile(pdf_filepath):
                self.df.at[index, 'image_extraction_complete'] = 'pdf_not_found'
                continue
            
            try:             
                # Initialize the ImageExtractor with the PDF file path
                extractor = ImageExtractor(pdf_file_path=pdf_filepath)
                
                # Extract images
                extractor.extract_images()
                
                # Retrieve the list of extracted image paths
                extracted_images = extractor.img_paths
                
                # Update the DataFrame with the extracted image paths
                self.df.at[index, 'image_paths'] = extracted_images
                self.df.at[index, 'image_extraction_complete'] = 'complete' if extracted_images else 'no_images_found'
            
            except Exception as e:
                print(f"Error extracting images for PMID {row.get('pmid', 'Unknown')}: {e}")
                self.df.at[index, 'image_extraction_complete'] = 'failed'
            
            # Optionally, save progress after each extraction
            self.save()
        
        print("Image extraction process completed.")
        return self

    def fetch_references(self):
        """
        Fetches references for each article in the DataFrame using the ReferenceRetriever class.

        Process:
        1. Checks if the DataFrame is populated and contains the necessary identifiers ('doi' or 'pmid').
        2. Iterates through articles, skipping those that already have references.
        3. Uses ReferenceRetriever to fetch references based on DOI or PMID.
        4. Updates the DataFrame with retrieved references.
        5. Saves progress after processing each article.

        Returns:
        - self: The updated PubMedSearcher instance with references added.

        Notes:
        - Prioritizes DOI-based retrieval when available.
        - Uses multiple sources (PubMed, PMC, Europe PMC, CrossRef).
        - References are stored as lists of dictionaries in the 'references' column.
        """
        if self.df.empty:
            print("DataFrame is empty. No articles to fetch references for.")
            return

        if 'references' not in self.df.columns:
            self.df['references'] = None

        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0], desc="Fetching References"):
            if pd.notna(row['references']):  # Skip if references already exist
                continue
            
            # Initialize ReferenceRetriever with available identifiers
            retriever = ReferenceRetriever(email=self.email, doi=row.get('doi'), pmid=row.get('pmid'), standardize=True)
            references = retriever.fetch_references()
            
            # Store references or mark as "Not found" if empty
            self.df.at[index, 'references'] = references if references else "Not found"

            # Save progress after each update
            self.save()

        return self

    def fetch_cited_by(self):
        """
        Fetches articles that cite each article in the DataFrame using ReferenceRetriever.

        Process:
        1. Ensures the DataFrame contains necessary identifiers ('pmid' or 'doi').
        2. Iterates through each article, skipping those that already have citing articles.
        3. Uses ReferenceRetriever to fetch citing articles from Europe PMC and PubMed.
        4. Updates the DataFrame with citing article information.
        5. Saves progress after processing each article.

        Returns:
        - self: The updated PubMedSearcher instance with citing articles added.

        Notes:
        - Citing articles are retrieved from Europe PMC (default) and PubMed if available.
        - Works best for articles with an existing record in Europe PMC.
        """
        if self.df.empty:
            print("DataFrame is empty. No articles to fetch cited-by data for.")
            return

        if 'cited_by' not in self.df.columns:
            self.df['cited_by'] = None

        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0], desc="Fetching Cited By"):
            if pd.notna(row['cited_by']):  # Skip if already populated
                continue
            
            # Initialize ReferenceRetriever with available identifiers
            retriever = ReferenceRetriever(email=self.email, doi=row.get('doi'), pmid=row.get('pmid'))
            cited_by = retriever.fetch_cited_by()  # Now uses both Europe PMC & PubMed
            
            # Store citing articles or mark as "Not found" if empty
            self.df.at[index, 'cited_by'] = cited_by if cited_by else "Not found"

            # Save progress after each update
            self.save()

        return self

    def fetch_abstracts(self):
        """
        Fetches abstracts for each article in the DataFrame using the Entrez API.
        Unnecessary if you used the 'search' method to retrieve articles, as abstracts are already included.
        """
        if not hasattr(self, 'df') or self.df.empty:
            print("DataFrame does not exist or is empty.")
            return
        
        if 'abstract' not in self.df.columns:
            self.df['abstract'] = None
        
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0], desc="Fetching Abstracts"):
            if pd.notna(row['abstract']):
                continue
            pmid = row.get('pmid')
            if pd.notna(pmid):
                abstract = self.get_abstract(pmid)
                self.df.at[index, 'abstract'] = abstract
            
            self.save()
    
    def get_abstract(self, pmid):
        """Fetches the abstract for an article identified by its PMID using the Entrez API."""
        Entrez.email = self.email
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        article_details = Entrez.read(handle)
        handle.close()
        abstract = article_details['PubmedArticle'][0]['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', '')
        return " ".join(abstract)

    def download_xml_fulltext(self, download_directory="downloads"):
        """
        Downloads the XML full text for each article in the DataFrame to the specified directory.

        Parameters:
        download_directory (str): The directory where the XML files should be saved.

        XML full text download is not very common (it has to be in the pubmed OA subset, either USA or European)
        """
        if not hasattr(self, 'df') or self.df.empty:
            print("DataFrame does not exist or is empty.")
            return
        
        if 'xml_download_complete' not in self.df.columns:
            self.df['xml_download_complete'] = 'Not started'
        if 'xml_filepath' not in self.df.columns:
            self.df['xml_filepath'] = None

        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0], desc="Downloading XML full texts"):
            if row.get('xml_download_complete') == 'Complete':
                continue
            
            custom_download_dir = self._determine_download_directory(row, download_directory, index)
            os.makedirs(custom_download_dir, exist_ok=True)

            last_name = row.get('first_author', '').split(',')[0].strip() if 'first_author' in row else None
            year = str(row.get('publication_year')) if 'publication_year' in row else None
            filename_suffix = f"{last_name}_{year}.xml" if last_name and year else None

            # Attempt to download XML full text
            if row.get('is_oa', False):
                file_path = None
                if pd.notna(row.get('europe_pmc_url')):
                    file_path = self.download_article_xml_europe(row.get('pmid'), custom_download_dir, filename_suffix)
                elif pd.notna(row.get('pmcid')):
                    file_path = self.download_article_xml_pubmed_oa_subset(row.get('pmcid'), custom_download_dir, filename_suffix)
                
                if file_path:
                    self.df.at[index, 'xml_download_complete'] = 'Complete'
                    self.df.at[index, 'xml_filepath'] = file_path
                else:
                    self.df.at[index, 'xml_download_complete'] = "Unavailable"
                    self.df.at[index, 'xml_filepath'] = None
            else:
                self.df.at[index, 'xml_download_complete'] = "Not OA or no XML available"
                self.df.at[index, 'xml_filepath'] = None

    def save(self, csv_path="master_list.csv"):
        """
        Saves the DataFrame to a CSV file.
        """
        self.df.to_csv(csv_path, index=False)

    def save_abstracts_as_csv(self, filename="abstracts.csv"):
        """Saves a DataFrame containing only the 'pmid' and 'abstract' columns to a CSV file."""
        abstracts_df = self.df[['pmid','abstract']].copy()
        abstracts_df.to_csv(filename, index=False)
    
    def download_article_xml_europe(self, pmid, download_directory="downloads", filename_suffix=None):
        """
        Downloads the XML data for an article identified by its PMID from Europe PMC and saves it to a specified directory.
        """
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/MED/{pmid}/fullTextXML"
        print(url)

        try:
            response = requests.get(url)
            if response.status_code == 200:
                xml_content = response.text
                os.makedirs(download_directory, exist_ok=True)
                # Use filename_suffix if provided, else default to PMID
                filename = f"{filename_suffix}.xml" if filename_suffix else f"{pmid}.xml"
                file_path = os.path.join(download_directory, filename)
                
                with open(file_path, 'w', encoding='utf-8') as file:
                    file.write(xml_content)
                
                print(f"Article XML downloaded successfully to {file_path}.")
                return file_path
            else:
                print(f"Failed to fetch article XML. Status code: {response.status_code}")
                return None
        except Exception as e:
            print(f"Exception occurred while fetching article XML: {e}")
            return None

    def download_article_xml_pubmed_oa_subset(self, pmcid, download_directory="downloads", filename_suffix=None):
        """
        Downloads the XML data for an article identified by its PMCID from PubMed OA subset and saves it to a specified directory.
        """
        xml_content = self._get_xml_for_pmcid(pmcid)
        if xml_content:
            os.makedirs(download_directory, exist_ok=True)
            # Use filename_suffix if provided, else default to PMCID
            filename = f"{filename_suffix}.xml" if filename_suffix else f"{pmcid}.xml"
            file_path = os.path.join(download_directory, filename)
            with open(file_path, 'wb') as file:
                file.write(xml_content)
            print(f"Article XML downloaded successfully to {file_path}.")
            return file_path
        else:
            print(f"Failed to download article XML for PMCID {pmcid}.")
            return None

    def _validate_dataframe(self, df):
        """Validates an input DataFrame to ensure it contains the required columns.
        Required columns: 'title', 'doi' # We may need to adjust this in the future
        """
        required_columns = ['title', 'doi']  # Adjusted required columns
        df.rename(columns={col: col.lower() for col in df.columns}, inplace=True)
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"DataFrame is missing required columns: {', '.join(missing_columns)}")

    def _parse_records_to_df(self, records_xml_bytes):
        """ Parses the API result Entrez into a DataFrame."""
        records_io = BytesIO(records_xml_bytes)
        records = Entrez.read(records_io)
        records_data = []

        for record in records['PubmedArticle']:
            article_data = {}
            medline = record['MedlineCitation']
            article = medline['Article']
            article_data['title'] = article.get('ArticleTitle', '')
            authors_list = article.get('AuthorList', [])
            authors = [f"{a.get('LastName', '')}, {a.get('ForeName', '')}" for a in authors_list]
            article_data['authors'] = "; ".join(authors)
            article_data['first_author'] = authors[0].split(',')[0] if authors else ''
            article_data['abstract'] = " ".join(article.get('Abstract', {}).get('AbstractText', []))
            publication_date = article.get('ArticleDate', [])
            article_data['publication_date'] = publication_date[0] if publication_date else {}
            article_data['publication_year'] = publication_date[0]['Year'] if publication_date else None
            article_data['journal_info'] = article.get('Journal', {}).get('Title', '')
            article_id_dict = {article_id.attributes['IdType']: str(article_id) for article_id in record.get('PubmedData', {}).get('ArticleIdList', [])}
            doi = article_id_dict.get('doi', "")
            pmcid = article_id_dict.get('pmc', "")
            pmid = medline.get('PMID', "")
            article_data['doi'] = doi
            article_data['pmcid'] = pmcid
            article_data['pmid'] = pmid

            keywords = medline.get('KeywordList', [])
            article_data['keywords'] = "; ".join([kwd for sublist in keywords for kwd in sublist]) if keywords else ""
            
            publication_types = article.get('PublicationTypeList', [])
            article_data['article_type'] = "; ".join([ptype for ptype in publication_types])
            medline_journal_info = medline.get('MedlineJournalInfo', {})
            article_data['country'] = medline_journal_info.get('Country', "")
            article_data['language'] = "; ".join(article.get('Language', []))

            records_data.append(article_data)

        df =  pd.DataFrame(records_data)
        cols = df.columns.tolist()
        cols.insert(0, cols.pop(cols.index('pmid')))
        df = df[cols]
        return df

    def _get_xml_for_pmcid(self, pmcid):
        """
        Fetches XML content for a given PMCID from PubMed Central.
        """
        pmcid = pmcid.replace("PMC", "")
        base_url = "https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi"
        params = {
            "verb": "GetRecord",
            "identifier": "oai:pubmedcentral.nih.gov:" + pmcid,
            "metadataPrefix": "pmc"
        }
        response = requests.get(base_url, params=params)

        if response.status_code == 200:
            if 'is not supported by the item or by the repository.' in response.content.decode():
                print("Error: This PMC is not open access through PubMed Central, or the ID is invalid.")
                return None
            return response.content
        else:
            print("Error:", response.status_code)
            return None
    