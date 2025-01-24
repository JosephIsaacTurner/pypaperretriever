import os
import re

import numpy as np
import pandas as pd
import requests
from Bio import Entrez
from io import BytesIO
from itertools import chain
from lxml import etree
from tqdm import tqdm

import fitz  # PyMuPDF
# from calvin_utils.gpt_sys_review.image_utils import ImageExtractor
from pypaperretriever import PaperRetriever

tqdm.pandas()

class PubMedSearcher:
    """
    A class to search PubMed for articles based on a search string and retrieve article information.

    Attributes:
    - search_string (str): The search string to use when querying PubMed.
    - df (DataFrame): A DataFrame containing the retrieved articles.
    - email (str): The email address to use when querying PubMed.

    User-Facing Methods:
    - search: Searches PubMed for articles based on the search string and retrieves the specified number of articles.
    - download_articles: Downloads articles from the DataFrame to the specified directory (open access is prioritized, but may use PyPaperBot as a fallback).
    - fetch_references: Fetches references for each article in the DataFrame using multiple methods.
    - standardize_references: Standardizes the references column in the DataFrame to only contain the following keys: ['doi', 'pmid', 'pmcid', 'title', 'authors']
    - fetch_cited_by: Fetches list of articles that cite each article in the DataFrame using Europe PMC (only works for articles with a record in Europe PMC)
    - download_xml_fulltext: Downloads the XML full text for each article in the DataFrame to the specified directory (rarely available, but can be useful).
    - save: Saves the df to a CSV file.
    """

    def __init__(self, search_string=None, df=None, email="lesion_bank@gmail.com"):
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
        self.email = email

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

    def download_articles(self, allow_scihub=True, download_directory="pdf_downloads", max_articles=None):
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
            pdf_filepath = PaperRetriever(
                                        pmid=pmid,
                                        doi=doi,
                                        email=self.email,
                                        allow_scihub=allow_scihub,
                                        download_directory=download_directory
                                    ).find_and_download().filepath
            if pdf_filepath in [None, '', 'unavailable']:
                self.df.at[index, 'download_complete'] = 'unavailable'
                self.df.at[index, 'pdf_filepath'] = None
            else:
                self.df.at[index, 'download_complete'] = 'complete'
                self.df.at[index, 'pdf_filepath'] = pdf_filepath
            articles_processed += 1
            self.save()
        return self

    def fetch_references(self):
        """
        Fetches references for each article in the DataFrame using the find_references method.

        The find_references method will attempt to fetch references in the following order:
        1. PubMed
        2. PMC
        3. Europe PMC
        4. CrossRef
        If no references are found using these methods, it will return "Not found".
        """
        if not hasattr(self, 'df') or self.df.empty:
            print("DataFrame does not exist or is empty.")
            return
        
        if 'references' not in self.df.columns:
            self.df['references'] = None
        
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0], desc="Fetching References"):
            if pd.isna(row['references']):
                references = self._find_references_for_row(row)
                if not references:  # If references list is empty or None
                    references = "Not found"
                self.df.at[index, 'references'] = references

    def standardize_references(self):
        """
        Standardizes the references column in the DataFrame to only contain the following keys:
        ['doi', 'pmid', 'pmcid', 'title', 'authors']
        Populates a new column 'references_standardized' with the standardized references (list of dicts)
        """
        def standardize_references_for_row(references):
            standard_keys = ['doi', 'pmid', 'pmcid', 'title', 'authors']
            return [
                {key: ref.get(key, None) for key in standard_keys}
                for ref in references
                if isinstance(ref, dict)
            ]

        if 'references' not in self.df.columns:
            print('Error: No references column found in DataFrame.')
            return

        if 'references_standardized' not in self.df.columns:
            self.df['references_standardized'] = None

        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0], desc="Standardizing references"):
            ref_standardized = row['references_standardized']
            if pd.notna(ref_standardized) and (isinstance(ref_standardized, (list, np.ndarray)) and len(ref_standardized) > 0):
                continue
            
            references = row['references']
            if isinstance(references, list) and not references:
                continue
            if isinstance(references, np.ndarray) and references.size == 0:
                continue
            if not isinstance(references, (list, np.ndarray)) and pd.isna(references):
                continue

            self.df.at[index, 'references_standardized'] = standardize_references_for_row(references)

    def fetch_cited_by(self):
        """
        Fetches list of articles that cite each article in the DataFrame using Europe PMC.
        Currently only works for articles with a record in Europe PMC.
        """
        if not hasattr(self, 'df') or self.df.empty:
            print("DataFrame does not exist or is empty.")
            return
        
        if 'cited_by' not in self.df.columns:
            self.df['cited_by'] = None
        
        for index, row in tqdm(self.df.iterrows(), total=self.df.shape[0], desc="Fetching Cited By"):
            if pd.notna(row['cited_by']):
                continue
            if pd.notna(row.get('is_oa')) and row.get('is_oa') and pd.notna(row.get('europe_pmc_url')):
                cited_by = self.get_citing_articles_europe(row.get('pmid'))
                self.df.at[index, 'cited_by'] = cited_by

    def fetch_abstracts(self, save_progress=True):
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
            
            if save_progress:
                self.save()
    
    def get_abstract(self, pmid):
        """Fetches the abstract for an article identified by its PMID using the Entrez API."""
        Entrez.email = self.email
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        article_details = Entrez.read(handle)
        handle.close()
        abstract = article_details['PubmedArticle'][0]['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', '')
        return " ".join(abstract)

    # def extract_images(self):
    #     """
    #     Extracts images from PDFs for each article in the DataFrame and saves the image paths in a new column 'img_paths'.
    #     Able to handle both native and image-based PDFs.
    #     """
    #     if 'pdf_filepath' not in self.df.columns:
    #         print("PDF file paths column not found in DataFrame.")
    #         return
    #     self.df["img_paths"] = [[] for _ in range(len(self.df))]
    #     for idx, row in tqdm(self.df.iterrows(), total=len(self.df), desc="Extracting images"):
    #         # Check if img_paths already exists, or is None, and skip if so
    #         if pd.notna(row["img_paths"]) and len(row["img_paths"]) > 0:
    #             continue
    #         pdf_path = row["pdf_filepath"]
    #         if pd.isna(pdf_path) or pdf_path in [0, "0", None, "", "<NA>"]:
    #             continue
    #         img_extractor = ImageExtractor(pdf_path)
    #         img_extractor.extract_images()
    #         self.df.at[idx, "img_paths"] = img_extractor.img_paths if len(img_extractor.img_paths) > 0 else ["Unavailable"]
    #     return self

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

    def get_references_europe(self, pmid):
        """
        Fetches references for an article identified by its PMID from Europe PMC API (hits the MEDLINE database).
        """
    
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/MED/{pmid}/references?page=1&pageSize=1000&format=json"

        try:
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                data =  data.get('referenceList', [])
                if data:
                    parsed_data = data.get('reference', [])
                    if "id" in parsed_data:
                        parsed_data["pmid"] = parsed_data.pop("id")
                    if "authorString" in data:
                        parsed_data["authors"] = parsed_data.pop("authorString")
                return parsed_data.get('reference', []) if data else []
            else:
                # print(f"Failed to fetch references. Status code: {response.status_code}")
                return None
        except Exception as e:
            # print(f"Exception occurred while fetching references: {e}")
            return None
        
    def get_references_entrez_pmc(self, pmcid):
        """Finds references for a given PMCID using Entrez API.
        Seems to return identical results to the get_references_pubmed_oa_subset method.
        """
        Entrez.email = self.email
        handle = Entrez.efetch(db="pmc", id=pmcid, retmode="xml")
        xml_data = handle.read()
        handle.close()
        references = self._parse_pubmed_references(xml_data)
        return references
        
    def get_references_pubmed_oa_subset(self, pmcid):
        xml_content = self._get_xml_for_pmcid(pmcid)
        if xml_content:
            references = self._parse_pubmed_references(xml_content)
            return references
        else:
            return None
        
    def get_references_entrez_pubmed(self, pmid):
        """
        Returns a list of PMCIDs for the references of a given pmid. 
        Doesn't seem to work for all PMIDs, so use with caution.
        """
        Entrez.email = self.email
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        article_details = Entrez.read(handle)
        handle.close()

        if article_details['PubmedArticle'][0]['PubmedData']['ReferenceList'] and len(article_details['PubmedArticle'][0]['PubmedData']['ReferenceList']) > 0:
            references = []
            try:
                # Attempt to navigate to the ReferenceList
                result_list = article_details['PubmedArticle'][0]['PubmedData']['ReferenceList'][0]['Reference']
                authors_pattern = r"^(.*?)\s+et al\."
                doi_pattern = r"doi\s*:\s*([^\s.]+)\.?"
                doi_pattern2 = r"doi\.org/([^\s,;]+)"
                doi_pattern3 = r"doi\.wiley\.org/([^\s,;]+)"

                for ref in result_list:
                    article_id_list = ref.get('ArticleIdList', [])
                    citation = ref.get('Citation', '')
                    ref_dict = {'citation': citation}

                    if article_id_list:
                        for element in article_id_list:
                            value = str(element)
                            id_type = element.attributes['IdType']
                            ref_dict[id_type] = value

                    if 'doi' not in ref_dict:
                        match = re.search(doi_pattern, citation, re.IGNORECASE)
                        if match:
                            ref_dict['doi'] = match.group(1)
                        elif 'doi' not in ref_dict:
                            match2 = re.search(doi_pattern2, citation, re.IGNORECASE)
                            if match2:
                                ref_dict['doi'] = match2.group(1)
                        else:
                            match3 = re.search(doi_pattern3, citation, re.IGNORECASE)
                            if match3:
                                ref_dict['doi'] = match3.group(1)

                    authors_match = re.search(authors_pattern, citation, re.IGNORECASE)
                    if authors_match:
                        ref_dict['authors'] = authors_match.group(1)

                    if 'pubmed' in ref_dict:
                        ref_dict['pmid'] = ref_dict.pop('pubmed')
                    if 'pmc' in ref_dict:
                        ref_dict['pmcid'] = ref_dict.pop('pmc')
                    
                    references.append(ref_dict)
                return references

            except (KeyError, IndexError, TypeError) as e:
                print(f"Error navigating article details: {e}")
                return None
        return None
        
    def get_references_crossref(self, doi):
        """
        Fetches references for a given DOI using the CrossRef REST API and formats them into a pandas DataFrame.
        
        Parameters:
            doi (str): The DOI of the article for which to fetch references.
            
        Returns:
            DataFrame: A pandas DataFrame containing the references, with each column representing a common key.
        """
        base_url = "https://api.crossref.org/works/"
        full_url = f"{base_url}{doi}"
        
        try:
            response = requests.get(full_url)
            if response.status_code == 200:
                data = response.json()
                references = data['message'].get('reference', [])
                if references:
                    for reference in references:
                        if 'DOI' in reference:
                            reference['doi'] = reference.pop('DOI')
                        if 'author' in reference:
                            reference['authors'] = reference.pop('author')
                        if 'article-title' in reference:
                            reference['title'] = reference.pop('article-title')
                    return references
                else:
                    return None
            else:
                print(f"Failed to fetch references, HTTP status code: {response.status_code}")
                return None
        except requests.exceptions.RequestException as e:
            print(f"Request failed: {str(e)}")
            return None

    def get_citing_articles_europe(self, pmid):
        """
        Fetches references for an article identified by its PMID from Europe PMC.
        Tries two different search methods and returns the results from the first successful one.
        """

        def try_restful_search(pmid):
            url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/MED/{pmid}/citations?format=json"
            try:
                response = requests.get(url)
                if response.status_code == 200:
                    data = response.json()
                    if data.get("citationList", {}).get("citation"):
                        parsed_data = data['citationList']
                        for citation in parsed_data["citation"]:
                            if "id" in citation:
                                citation["pmid"] = citation.pop("id")
                            if "authorString" in citation:
                                citation["authors"] = citation.pop("authorString")
                        return parsed_data["citation"]
                else:
                    print(f"Failed to fetch references. Status code: {response.status_code}")
            except Exception as e:
                print(f"Exception occurred while fetching references: {e}")
            return None

        def try_query_search(pmid):
            url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=cites:{pmid}_MED&format=json"
            try:
                response = requests.get(url)
                if response.status_code == 200:
                    data = response.json()
                    if data.get("resultList", {}).get("result"):
                        return data
                else:
                    print(f"Failed to fetch references. Status code: {response.status_code}")
            except Exception as e:
                print(f"Exception occurred while fetching references: {e}")
            return None

        # Try the RESTful search first
        restful_result = try_restful_search(pmid)
        if restful_result:
            return restful_result
        
        # If the RESTful search fails to provide results, try the query search
        query_result = try_query_search(pmid)
        if query_result:
            return query_result
        
        # If both searches fail, return a message indicating no results were found
        return "No citing articles found."
    
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
    
    def _find_references_for_row(self, row):
        pmcid = row.get('pmcid', None)  # Correct usage for a pandas Series
        pmid = row.get('pmid', None)
        doi = row.get('doi', None)

        def try_pubmed(pmid):
            if pmid:
                references = self.get_references_entrez_pubmed(pmid)
                if references is not None and len(references) > 0:
                    return references

        def try_pmc(pmcid):
            if pmcid:
                references = self.get_references_entrez_pmc(pmcid)
                if references is not None and len(references) > 0:
                    return references
                else:
                    references = self.get_references_pubmed_oa_subset(pmcid)
                    if references is not None and len(references) > 0:
                        return references

        def try_europe(pmid):
            if pmid:
                references = self.get_references_europe(pmid)
                if references is not None and len(references) > 0:
                    return references

        def try_crossref(doi):
            if doi:
                references = self.get_references_crossref(doi)
                if references is not None and len(references) > 0:
                    return references

        # Try all methods until one works
        references = try_pubmed(pmid)
        if references:
            return references

        references = try_pmc(pmcid)
        if references:
            return references

        references = try_europe(pmid)
        if references:
            return references

        references = try_crossref(doi)
        if references:
            return references

        # Return None or an empty list if no references found
        return []
        
    def _stringify_children(self, node):
        """
        Filters and removes possible Nones in texts and tails
        ref: http://stackoverflow.com/questions/4624062/get-all-text-inside-a-tag-in-lxml
        """
        parts = (
            [node.text]
            + list(chain(*([c.text, c.tail] for c in node.getchildren())))
            + [node.tail]
        )
        return "".join(filter(None, parts))

    def _parse_article_meta(self, tree):
        """
        Parse PMID, PMC and DOI from given article tree
        """
        article_meta = tree.find(".//article-meta")
        if article_meta is not None:
            pmid_node = article_meta.find('article-id[@pub-id-type="pmid"]')
            pmc_node = article_meta.find('article-id[@pub-id-type="pmc"]')
            pub_id_node = article_meta.find('article-id[@pub-id-type="publisher-id"]')
            doi_node = article_meta.find('article-id[@pub-id-type="doi"]')
        else:
            pmid_node = None
            pmc_node = None
            pub_id_node = None
            doi_node = None

        pmid = pmid_node.text if pmid_node is not None else ""
        pmc = pmc_node.text if pmc_node is not None else ""
        pub_id = pub_id_node.text if pub_id_node is not None else ""
        doi = doi_node.text if doi_node is not None else ""

        dict_article_meta = {"pmid": pmid, "pmc": pmc, "doi": doi, "publisher_id": pub_id}

        return dict_article_meta

    def _remove_namespace(self, tree):
        """
        Strip namespace from parsed XML
        """
        for node in tree.iter():
            try:
                has_namespace = node.tag.startswith("{")
            except AttributeError:
                continue  # node.tag is not a string (node is a comment or similar)
            if has_namespace:
                node.tag = node.tag.split("}", 1)[1]

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
        
    def _parse_pubmed_references(self, xml_content):
        """
        Parse reference articles from XML content to list of dictionaries, 
        independent of namespace.
        
        Parameters
        ----------
        xml_content: bytes
            XML content as bytes.
        
        Returns
        -------
        DataFrame
            A DataFrame containing references made in the given XML.
        """
        
        tree = etree.fromstring(xml_content)
        self._remove_namespace(tree)
        dict_article_meta = self._parse_article_meta(tree)
        pmid = dict_article_meta["pmid"]
        pmc = dict_article_meta["pmc"]
        references = tree.xpath(".//ref-list/ref")
        dict_refs = []

        for reference in references:
            ref_id = reference.attrib.get("id")
            ref_type = reference.xpath(".//citation/@citation-type")
            journal_type = ref_type[0] if ref_type else ""

            # Extract names
            names = reference.xpath(".//person-group[@person-group-type='author']/name/surname/text()") + \
                    reference.xpath(".//person-group[@person-group-type='author']/name/given-names/text()")
            names = [" ".join(names[i:i+2]) for i in range(0, len(names), 2)]
            
            # Extract article title, source, year, DOI, PMID
            article_title = reference.xpath(".//article-title/text()")
            article_title = article_title[0].replace("\n", " ").strip() if article_title else ""
            
            journal = reference.xpath(".//source/text()")
            journal = journal[0] if journal else ""
            
            year = reference.xpath(".//year/text()")
            year = year[0] if year else ""
            
            doi_cited = reference.xpath(".//pub-id[@pub-id-type='doi']/text()")
            doi_cited = doi_cited[0] if doi_cited else ""
            
            pmid_cited = reference.xpath(".//pub-id[@pub-id-type='pmid']/text()")
            pmid_cited = pmid_cited[0] if pmid_cited else ""
            
            dict_ref = {
                "pmid": pmid,
                "pmc": pmc,
                "ref_id": ref_id,
                "pmid_cited": pmid_cited,
                "doi_cited": doi_cited,
                "title": article_title,
                "authors": "; ".join(names),
                "year": year,
                "journal": journal,
                "journal_type": journal_type,
            }
            
            dict_refs.append(dict_ref)
            
        return dict_refs