import requests
from Bio import Entrez
import re

from .utils import doi_to_pmid

class ReferenceRetriever:
    """
    A class to retrieve academic paper references and citations using DOI or PMID identifiers.
    
    Fetches data from multiple sources including PubMed, Europe PMC, and CrossRef.
    """
    
    def __init__(self, email, doi=None, pmid=None):
        """
        Initialize the retriever with either DOI or PMID.

        Args:
            email (str): Email for API access
            doi (str, optional): Digital Object Identifier
            pmid (str, optional): PubMed ID
        """
        self.email = email
        self.doi = doi
        self.pmid = pmid

        print(f"[ReferenceRetriever] Initializing with DOI: {doi} and PMID: {pmid}")

        if self.doi and not self.pmid:
            print(f"[ReferenceRetriever] Converting DOI to PMID for DOI: {self.doi}")
            self.pmid = doi_to_pmid(self.doi, self.email)
            print(f"[ReferenceRetriever] Converted DOI {self.doi} to PMID: {self.pmid}")

    def fetch_references(self):
        """
        Fetch references for the given DOI or PMID using multiple sources.

        Returns:
            list: List of dictionaries containing reference metadata. Empty list if none found.

        Raises:
            ValueError: If neither DOI nor PMID is provided.
        """
        print(f"[ReferenceRetriever] Fetching references for DOI: {self.doi}, PMID: {self.pmid}")
        if not self.doi and not self.pmid:
            raise ValueError("Either DOI or PMID must be provided.")
        
        references = self._find_references()
        if not references:
            print("[ReferenceRetriever] No references found.")
            return []
        print(f"[ReferenceRetriever] Found {len(references)} references.")
        return references

    def fetch_cited_by(self):
        """
        Fetch articles that cite the given DOI or PMID.

        Returns:
            list: List of dictionaries containing citing article metadata. Empty list if none found.

        Raises:
            ValueError: If PMID conversion fails or is not provided.
        """
        print(f"[ReferenceRetriever] Fetching citing articles for DOI: {self.doi}, PMID: {self.pmid}")
        if not self.pmid:
            if self.doi:
                print(f"[ReferenceRetriever] Converting DOI to PMID for DOI: {self.doi}")
                self.pmid = doi_to_pmid(self.doi, self.email)
                print(f"[ReferenceRetriever] Converted DOI {self.doi} to PMID: {self.pmid}")
                if not self.pmid:
                    raise ValueError("Unable to convert DOI to PMID.")
            else:
                raise ValueError("PMID must be provided to fetch cited_by articles.")
        
        cited_by = self._find_cited_by()
        if not cited_by:
            print("[ReferenceRetriever] No citing articles found.")
            return []
        print(f"[ReferenceRetriever] Found {len(cited_by)} citing articles.")
        return cited_by

    def get_paper_metadata(self):
        """
        Fetch metadata for the paper identified by DOI or PMID.

        Returns:
            dict: Paper metadata including DOI, PMID, title, authors, and year. Empty dict if not found.
        """
        print(f"[ReferenceRetriever] Fetching metadata for DOI: {self.doi}, PMID: {self.pmid}")
        if self.pmid:
            print(f"[ReferenceRetriever] Fetching metadata using PMID: {self.pmid}")
            articles = self._fetch_articles_details([self.pmid])
            if articles:
                print(f"[ReferenceRetriever] Metadata found for PMID: {self.pmid}")
                return articles[0]
            else:
                print(f"[ReferenceRetriever] No metadata found for PMID: {self.pmid}")
        elif self.doi:
            print(f"[ReferenceRetriever] Fetching metadata using DOI: {self.doi}")
            pmid = doi_to_pmid(self.doi, self.email)
            if pmid:
                print(f"[ReferenceRetriever] Converted DOI {self.doi} to PMID: {pmid}")
                articles = self._fetch_articles_details([pmid])
                if articles:
                    print(f"[ReferenceRetriever] Metadata found for PMID: {pmid}")
                    return articles[0]
                else:
                    print(f"[ReferenceRetriever] No metadata found for PMID: {pmid}")
            else:
                print(f"[ReferenceRetriever] Unable to convert DOI {self.doi} to PMID.")
        return {}

    def _find_references(self):
        """
        Find references using multiple sources (PubMed, Europe PMC, CrossRef).

        Returns:
            list: Combined list of references from all sources.
        """
        print("[ReferenceRetriever] Finding references using multiple sources.")
        references = []
        if self.pmid:
            print("[ReferenceRetriever] Fetching references from Entrez PubMed.")
            refs = self.get_references_entrez_pubmed(self.pmid)
            if refs:
                print(f"[ReferenceRetriever] Retrieved {len(refs)} references from Entrez PubMed.")
                references.extend(refs)
            
            print("[ReferenceRetriever] Fetching references from Europe PMC.")
            refs = self.get_references_europe(self.pmid)
            if refs:
                print(f"[ReferenceRetriever] Retrieved {len(refs)} references from Europe PMC.")
                references.extend(refs)
        
        if self.doi:
            print("[ReferenceRetriever] Fetching references from CrossRef.")
            refs = self.get_references_crossref(self.doi)
            if refs:
                print(f"[ReferenceRetriever] Retrieved {len(refs)} references from CrossRef.")
                references.extend(refs)
        
        print(f"[ReferenceRetriever] Total references found: {len(references)}")
        return references

    def _find_cited_by(self):
        """
        Find citing articles using multiple sources (PubMed, Europe PMC).

        Returns:
            list: Combined list of citing articles from all sources.
        """
        print("[ReferenceRetriever] Finding citing articles using multiple sources.")
        cited_by = []
        print("[ReferenceRetriever] Fetching citing articles from PubMed.")
        pubmed_cites = self.get_citing_articles_pubmed(self.pmid)
        if pubmed_cites:
            print(f"[ReferenceRetriever] Retrieved {len(pubmed_cites)} citing articles from PubMed.")
            cited_by.extend(pubmed_cites)

        print("[ReferenceRetriever] Fetching citing articles from Europe PMC.")
        europe_cites = self.get_citing_articles_europe(self.pmid)
        if europe_cites:
            print(f"[ReferenceRetriever] Retrieved {len(europe_cites)} citing articles from Europe PMC.")
            cited_by.extend(europe_cites)
        
        print(f"[ReferenceRetriever] Total citing articles found: {len(cited_by)}")
        return cited_by

    def get_references_europe(self, pmid):
        """
        Fetch references from Europe PMC API.

        Args:
            pmid (str): PubMed ID to fetch references for

        Returns:
            list: References found in Europe PMC. Empty list if none found or on error.
        """
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/MED/{pmid}/references?page=1&pageSize=1000&format=json"
        print(f"[ReferenceRetriever] Requesting Europe PMC references from URL: {url}")
        try:
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                references = data.get('referenceList', {}).get('reference', [])
                print(f"[ReferenceRetriever] Europe PMC returned {len(references)} raw references.")
                return self._parse_europe_references(references)
            else:
                print(f"[ReferenceRetriever] Europe PMC request failed with status code: {response.status_code}")
                return []
        except Exception as e:
            print(f"[ReferenceRetriever] Error fetching references from Europe PMC for PMID {pmid}: {e}")
            return []

    def get_references_entrez_pubmed(self, pmid):
        """
        Fetch references from PubMed using Entrez.

        Args:
            pmid (str): PubMed ID to fetch references for

        Returns:
            list: References found in PubMed. Empty list if none found or on error.
        """
        Entrez.email = self.email
        print(f"[ReferenceRetriever] Requesting Entrez PubMed references for PMID: {pmid}")
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            article_details = Entrez.read(handle)
            handle.close()
            
            references = []
            if 'PubmedArticle' in article_details:
                ref_list = article_details['PubmedArticle'][0].get('PubmedData', {}).get('ReferenceList', [])
                if ref_list:
                    references = ref_list[0].get('Reference', [])
                    parsed_refs = self._parse_pubmed_references(references)
                    print(f"[ReferenceRetriever] Entrez PubMed returned {len(parsed_refs)} parsed references.")
                    return parsed_refs
            print("[ReferenceRetriever] Entrez PubMed returned no references.")
            return []
        except Exception as e:
            print(f"[ReferenceRetriever] Error fetching references from PubMed for PMID {pmid}: {e}")
            return []

    def get_references_crossref(self, doi):
        """
        Fetch references from CrossRef API.

        Args:
            doi (str): DOI to fetch references for

        Returns:
            list: References found in CrossRef. Empty list if none found or on error.
        """
        url = f"https://api.crossref.org/works/{doi}"
        print(f"[ReferenceRetriever] Requesting CrossRef references from URL: {url}")
        try:
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                references = data['message'].get('reference', [])
                print(f"[ReferenceRetriever] CrossRef returned {len(references)} raw references.")
                return self._parse_crossref_references(references)
            else:
                print(f"[ReferenceRetriever] CrossRef request failed with status code: {response.status_code}")
                return []
        except Exception as e:
            print(f"[ReferenceRetriever] Error fetching references from CrossRef for DOI {doi}: {e}")
            return []

    def get_citing_articles_europe(self, pmid):
        """
        Fetch citing articles from Europe PMC API.

        Args:
            pmid (str): PubMed ID to fetch citations for

        Returns:
            list: Citing articles found in Europe PMC. Empty list if none found or on error.
        """
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/MED/{pmid}/citations?format=json"
        print(f"[ReferenceRetriever] Requesting Europe PMC citing articles from URL: {url}")
        try:
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                citations = data.get('citationList', {}).get('citation', [])
                print(f"[ReferenceRetriever] Europe PMC returned {len(citations)} raw citing articles.")
                return self._parse_europe_cited_by(citations)
            else:
                print(f"[ReferenceRetriever] Europe PMC citing articles request failed with status code: {response.status_code}")
                return []
        except Exception as e:
            print(f"[ReferenceRetriever] Error fetching citing articles from Europe PMC for PMID {pmid}: {e}")
            return []

    def get_citing_articles_pubmed(self, pmid):
        """
        Fetch citing articles from PubMed using Entrez.

        Args:
            pmid (str): PubMed ID to fetch citations for

        Returns:
            list: Citing articles found in PubMed. Empty list if none found or on error.
        """
        Entrez.email = self.email
        print(f"[ReferenceRetriever] Requesting PubMed citing articles for PMID: {pmid}")
        try:
            handle = Entrez.elink(dbfrom="pubmed", id=pmid, linkname="pubmed_pubmed_citedin")
            record = Entrez.read(handle)
            handle.close()
            
            citing_articles = []
            if 'LinkSetDb' in record[0]:
                for linksetdb in record[0]['LinkSetDb']:
                    if linksetdb.get('LinkName') == 'pubmed_pubmed_citedin':
                        citing_articles = linksetdb.get('Link', [])
                        break
            
            pmids = [link['Id'] for link in citing_articles]
            print(f"[ReferenceRetriever] PubMed returned {len(pmids)} citing PMIDs.")
            return self._fetch_articles_details(pmids)
        except Exception as e:
            print(f"[ReferenceRetriever] Error fetching citing articles from PubMed for PMID {pmid}: {e}")
            return []

    def _fetch_articles_details(self, pmids):
        """
        Fetch detailed metadata for multiple PMIDs from PubMed.

        Args:
            pmids (list): List of PubMed IDs to fetch details for

        Returns:
            list: Article details for each PMID. Empty list if none found or on error.
        """
        Entrez.email = self.email
        pmid_string = ",".join(pmids)
        print(f"[ReferenceRetriever] Fetching article details for PMIDs: {pmid_string}")
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid_string, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            citing_articles = self._parse_pubmed_articles(records)
            print(f"[ReferenceRetriever] Retrieved details for {len(citing_articles)} articles.")
            return citing_articles
        except Exception as e:
            print(f"[ReferenceRetriever] Error fetching details for PMIDs: {e}")
            return []


    def _parse_europe_references(self, references):
        """
        Parses a list of references from Europe PMC and extracts relevant information.
        Args:
            references (list): A list of reference dictionaries obtained from Europe PMC.
        Returns:
            list: A list of dictionaries, where each dictionary contains parsed
                    information about a single reference, including authors, title,
                    journal, year, and DOI.
        """
        print("[ReferenceRetriever] Parsing Europe PMC references.")
        parsed_references = []
        for ref in references:
            parsed_ref = {
                'pmid': ref.get('id', ''),
                'authors': ref.get('authorString', ''),
                'title': ref.get('title', ''),
                'journal': ref.get('journal', ''),
                'year': ref.get('year', ''),
                'doi': ref.get('doi', '')
            }
            parsed_references.append(parsed_ref)
        print(f"[ReferenceRetriever] Parsed {len(parsed_references)} Europe PMC references.")
        return parsed_references

    def _parse_pubmed_references(self, references):
        """
        Parses a list of PubMed references and extracts relevant information.
        This method processes a list of references obtained from Entrez PubMed,
        extracting citation details, DOI, authors, and various identifiers such as
        PubMed ID (PMID) and PubMed Central ID (PMCID).
        Args:
            references (list): A list of reference dictionaries, where each dictionary
                               contains information about a single reference.
        Returns:
            list: A list of dictionaries, where each dictionary contains parsed
                  information about a single reference, including citation, DOI,
                  authors, PMID, and PMCID if available.
        """
        
        print("[ReferenceRetriever] Parsing Entrez PubMed references.")
        parsed_references = []
        for ref in references:
            citation = ref.get('Citation', '')
            authors_pattern = r"^(.*?)\s+et al\."
            doi_pattern = r"doi\s*:\s*([^\s.]+)\.?"
            
            ref_dict = {'citation': citation}
            article_id_list = ref.get('ArticleIdList', [])
            if article_id_list:
                for element in article_id_list:
                    value = str(element)
                    id_type = element.attributes.get('IdType')
                    if id_type:
                        ref_dict[id_type] = value

            match = re.search(doi_pattern, citation, re.IGNORECASE)
            if match:
                ref_dict['doi'] = match.group(1)
                
            authors_match = re.search(authors_pattern, citation, re.IGNORECASE)
            if authors_match:
                ref_dict['authors'] = authors_match.group(1)

            if 'pubmed' in ref_dict:
                ref_dict['pmid'] = ref_dict.pop('pubmed')
            if 'pmc' in ref_dict:
                ref_dict['pmcid'] = ref_dict.pop('pmc')
                
            parsed_references.append(ref_dict)
        print(f"[ReferenceRetriever] Parsed {len(parsed_references)} Entrez PubMed references.")
        return parsed_references

    def _parse_crossref_references(self, references):
        """
        Parses a list of references from CrossRef and extracts relevant information.
        Args:  
            references (list): A list of reference dictionaries obtained from CrossRef.
        Returns:
            list: A list of dictionaries, where each dictionary contains parsed
                    information about a single reference, including authors, title,
                    journal, year, and DOI.
        """
        print("[ReferenceRetriever] Parsing CrossRef references.")
        parsed_references = []
        for ref in references:
            parsed_ref = {
                'doi': ref.get('DOI', ''),
                'authors': self._format_crossref_authors(ref.get('author')),
                'title': ref.get('article-title', ''),
                'journal': ref.get('journal-title', ''),
                'year': ref.get('year', '')
            }
            parsed_references.append(parsed_ref)
        print(f"[ReferenceRetriever] Parsed {len(parsed_references)} CrossRef references.")
        return parsed_references

    def _parse_pubmed_articles(self, records):
        """
        Parses a list of PubMed articles and extracts relevant information.
        Args:
            records (dict): A dictionary containing PubMed articles data.
        Returns:    
            list: A list of dictionaries, where each dictionary contains parsed
                    information about a single article, including authors, title,
                    journal, year, and DOI.
        """
        print("[ReferenceRetriever] Parsing PubMed articles details.")
        parsed_articles = []
        for article in records.get('PubmedArticle', []):
            article_data = {
                'pmid': str(article['MedlineCitation']['PMID']),
                'title': article['MedlineCitation']['Article']['ArticleTitle'],
                'authors': self._get_author_list(article['MedlineCitation']['Article'].get('AuthorList', [])),
                'journal': article['MedlineCitation']['Article']['Journal']['Title'],
                'year': self._get_pub_date_year(article['MedlineCitation']['Article']['Journal']['JournalIssue'].get('PubDate', {})),
                'doi': None
            }
            for id in article['PubmedData']['ArticleIdList']:
                if id.attributes.get('IdType') == 'doi':
                    article_data['doi'] = str(id)
            parsed_articles.append(article_data)
        print(f"[ReferenceRetriever] Parsed details for {len(parsed_articles)} PubMed articles.")
        return parsed_articles

    def _get_author_list(self, authors):
        """
        Generate a comma-separated list of author names.
        Args:
            authors (list): List of dictionaries containing author details.
        Returns:
            str: Comma-separated string of author names in "LastName Initials" format.
        """

        author_list = []
        for author in authors:
            if 'LastName' in author and 'Initials' in author:
                name = f"{author['LastName']} {author['Initials']}"
                author_list.append(name)
        return ', '.join(author_list)

    def _get_pub_date_year(self, pub_date):
        """
        Extract the year from the publication date data.
        """
        if 'Year' in pub_date:
            return pub_date['Year']
        elif 'MedlineDate' in pub_date:
            # Extract year from MedlineDate
            match = re.search(r'\d{4}', pub_date['MedlineDate'])
            if match:
                return match.group(0)
        return ''

    def _parse_europe_cited_by(self, cited_by):
        """
        Parses a list of citing articles from Europe PMC and extracts relevant information.
        """
        print("[ReferenceRetriever] Parsing Europe PMC citing articles.")
        parsed_citations = []
        for citation in cited_by:
            parsed_citation = {
                'pmid': citation.get('id', ''),
                'authors': citation.get('authorString', ''),
                'title': citation.get('title', ''),
                'journal': citation.get('journalAbbreviation', ''),
                'year': citation.get('pubYear', ''),
                'doi': citation.get('doi', '')  # Assuming DOI might be provided
            }
            parsed_citations.append(parsed_citation)
        print(f"[ReferenceRetriever] Parsed {len(parsed_citations)} Europe PMC citing articles.")
        return parsed_citations

    def _format_crossref_authors(self, authors):
        """
        Format authors list from CrossRef API response.
        """
        if not authors:
            return ''
        formatted_authors = []
        for author in authors:
            if 'given' in author and 'family' in author:
                formatted_authors.append(f"{author['family']} {author['given']}")
        print(f"[ReferenceRetriever] Formatted {len(formatted_authors)} CrossRef authors.")
        return ', '.join(formatted_authors)