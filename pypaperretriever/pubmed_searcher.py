"""High level utilities for interacting with PubMed."""

from __future__ import annotations

import os
from io import BytesIO

import pandas as pd
import requests
from Bio import Entrez
from tqdm import tqdm
from typing import Self

from .image_extractor import ImageExtractor
from .paper_retriever import PaperRetriever
from .reference_retriever import ReferenceRetriever

tqdm.pandas()


class PubMedSearcher:
    """Search PubMed and manage retrieved articles.

    Args:
        search_string (str | None): Query used for PubMed search.
        df (pandas.DataFrame | None): Existing table of articles.
        email (str): Email address required by Entrez.

    Attributes:
        df (pandas.DataFrame): Table of article metadata and processing flags.
        search_string (str | None): Stored search query.
        email (str): Email address used for API calls.
    """

    def __init__(self, search_string=None, df=None, email=""):
        """Initialize the searcher.

        Args:
            search_string (str | None): Query to submit to PubMed.
            df (pandas.DataFrame | None): Existing table of articles.
            email (str): Email address required by Entrez.
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

    def search(
        self,
        count: int = 10,
        min_date: int | None = None,
        max_date: int | None = None,
        order_by: str = "chronological",
        only_open_access: bool = False,
        only_case_reports: bool = False,
    ) -> Self:
        """Search PubMed for articles.

        Args:
            count (int): Number of articles to retrieve.
            min_date (int | None): Minimum publication year.
            max_date (int | None): Maximum publication year.
            order_by (str): ``"chronological"`` or ``"relevance"``.
            only_open_access (bool): If ``True``, restrict to open-access articles.
            only_case_reports (bool): If ``True``, restrict to case reports.

        Returns:
            Self: This instance.

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
        return self

    def download_articles(
        self,
        allow_scihub: bool = False,
        download_directory: str = "pdf_downloads",
        max_articles: int | None = None,
    ) -> Self:
        """Download full-text PDFs for articles in ``df``.

        Args:
            allow_scihub (bool): Use Sci-Hub as a fallback source.
            download_directory (str): Directory to store downloaded PDFs.
            max_articles (int | None): Maximum number of articles to process.

        Returns:
            Self: The updated instance.

        """
        if self.df.empty:
            print("DataFrame is empty.")
            return self
        if 'pmid' not in self.df.columns or 'doi' not in self.df.columns:
            print("DataFrame is missing required columns for article download (pmid and doi)")
            return self
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
    
    def extract_images(self) -> Self:
        """Extract images from downloaded PDFs using :class:`ImageExtractor`.

        Only rows marked as successfully downloaded are processed.

        Returns:
            Self: The updated instance.

        """
        if self.df.empty:
            print("DataFrame is empty. No articles to extract images from.")
            return self

        required_columns = ['download_complete', 'pdf_filepath']
        for col in required_columns:
            if col not in self.df.columns:
                print(f"Error: DataFrame is missing required column '{col}'.")
                return self
        
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

    def fetch_references(self) -> Self:
        """Fetch references for each article in ``df``.

        Returns:
            Self: The updated instance.

        """
        if self.df.empty:
            print("DataFrame is empty. No articles to fetch references for.")
            return self

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

    def fetch_cited_by(self) -> Self:
        """Fetch citing articles for each entry in ``df``.

        Returns:
            Self: The updated instance.

        """
        if self.df.empty:
            print("DataFrame is empty. No articles to fetch cited-by data for.")
            return self

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
        """Retrieve abstracts for articles missing them in ``df``."""
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
    
    def get_abstract(self, pmid: str) -> str:
        """Fetch the abstract for a PMID.

        Args:
            pmid (str): Identifier of the article.

        Returns:
            str: Abstract text.
        """
        Entrez.email = self.email
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        article_details = Entrez.read(handle)
        handle.close()
        abstract = article_details['PubmedArticle'][0]['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', '')
        return " ".join(abstract)

    def download_xml_fulltext(self, download_directory: str = "downloads") -> Self:
        """Download XML full text for open-access articles.

        Args:
            download_directory (str): Destination directory for XML files.

        Returns:
            Self: The updated instance.

        """
        if not hasattr(self, 'df') or self.df.empty:
            print("DataFrame does not exist or is empty.")
            return self
        
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

        return self

    def save(self, csv_path: str = "master_list.csv") -> Self:
        """Persist the internal DataFrame to CSV.

        Args:
            csv_path (str): Output path for the CSV file.

        Returns:
            Self: This instance.

        """
        self.df.to_csv(csv_path, index=False)
        return self

    def save_abstracts_as_csv(self, filename: str = "abstracts.csv") -> Self:
        """Save only PMIDs and abstracts to a CSV file.

        Args:
            filename (str): Output filename.

        Returns:
            Self: This instance.

        """
        abstracts_df = self.df[['pmid', 'abstract']].copy()
        abstracts_df.to_csv(filename, index=False)
        return self

    def _determine_download_directory(
        self, row: pd.Series, base_directory: str, index: int
    ) -> str:
        """Build a per-article download directory.

        Args:
            row (pandas.Series): Row representing the article.
            base_directory (str): Root directory for downloads.
            index (int): Row index used as a fallback.

        Returns:
            str: Path to the directory where files should be saved.

        """
        identifier = row.get('pmid') or row.get('doi') or f"article_{index}"
        return os.path.join(base_directory, str(identifier))
    
    def download_article_xml_europe(
        self,
        pmid: str,
        download_directory: str = "downloads",
        filename_suffix: str | None = None,
    ) -> str | None:
        """Download XML full text from Europe PMC.

        Args:
            pmid (str): PubMed identifier.
            download_directory (str): Output directory.
            filename_suffix (str | None): Optional filename override.

        Returns:
            str | None: Path to the saved file or ``None`` on failure.
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

    def download_article_xml_pubmed_oa_subset(
        self,
        pmcid: str,
        download_directory: str = "downloads",
        filename_suffix: str | None = None,
    ) -> str | None:
        """Download XML from the PubMed Open Access subset.

        Args:
            pmcid (str): PubMed Central identifier.
            download_directory (str): Output directory.
            filename_suffix (str | None): Optional filename override.

        Returns:
            str | None: Path to the saved file or ``None`` on failure.
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

    def _validate_dataframe(self, df: pd.DataFrame) -> None:
        """Ensure an input DataFrame contains required columns."""
        required_columns = ['title', 'doi']  # Adjusted required columns
        df.rename(columns={col: col.lower() for col in df.columns}, inplace=True)
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"DataFrame is missing required columns: {', '.join(missing_columns)}")

    def _parse_records_to_df(self, records_xml_bytes: bytes) -> pd.DataFrame:
        """Convert Entrez XML bytes into a DataFrame."""
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

    def _get_xml_for_pmcid(self, pmcid: str) -> bytes | None:
        """Fetch XML content for a PMCID from PubMed Central.

        Args:
            pmcid (str): PubMed Central identifier.

        Returns:
            bytes | None: XML content if available.
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
