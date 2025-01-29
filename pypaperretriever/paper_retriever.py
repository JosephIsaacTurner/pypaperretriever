import os
import requests
import time
import random
from bs4 import BeautifulSoup
import json
from glob import glob
import fitz as pymupdf
import argparse
import re
from urllib.parse import urljoin, urlparse

from .utils import entrez_efetch, pmid_to_doi, doi_to_pmid, encode_doi, decode_doi

class PaperRetriever:
    """
    A class to find and download scientific papers.
    Attributes:
        email (str): The email address used for API requests.
        doi (str): The DOI of the paper.
        pmid (str): The PubMed ID of the paper.
        allow_scihub (bool): Whether to allow searching Sci-Hub for the paper.
        download_directory (str): Directory to save downloaded PDFs.
        filename (str): Filename for the downloaded PDF.
        is_oa (bool): Indicates if the paper is open access.
        on_scihub (bool): Indicates if the paper is available on Sci-Hub.
        pdf_urls (list): List of URLs to the PDF versions of the paper.
        is_downloaded (bool): Indicates if the paper has been downloaded.
        filepath (str): Path to the downloaded PDF file.
        override_previous_attempt (bool): Whether to override a previous download attempt.
        user_agents (list): List of user agents for making requests.
    Methods:
        download(): Finds and downloads the paper using the DOI or PMID.
        check_open_access(): Checks if the paper is open access using Unpaywall.
        check_pubmed_central_access(): Checks if the paper is available in PubMed Central.
        check_crossref_access(doi): Checks if the paper is available on Crossref.
        check_scihub_access(): Checks if the paper is available on Sci-Hub.
        _download_pdf(): Downloads the PDF from the first available URL.
        _create_json_sidecar(download_success, pdf_filepath, json_filepath, url=None): Creates a JSON sidecar file with download information.
        _determine_paths(): Determines the file paths for the download.
        _check_if_downloaded(download_directory_or_path, filetype=".pdf"): Checks if a non-corrupted file exists in the specified directory.
        _look_for_previous_download(): Looks for a previous download of the PDF file.
        _get_pdf_element(html_text, mirror): Extracts the PDF link from the HTML text of a Sci-Hub page.
    """

    def __init__(self, email, doi=None, pmid=None, allow_scihub=True, download_directory='PDFs', filename=None, override_previous_attempt=False):
        self.email = email
        if not doi and not pmid:
            raise ValueError("Either a DOI or PMID must be provided")
        if not doi and pmid:
            doi = pmid_to_doi(pmid, email)
        self.doi = encode_doi(doi)
        self.pmid = pmid
        self.allow_scihub = allow_scihub
        self.is_oa = False
        self.on_scihub = False
        self.pdf_urls = []
        self.is_downloaded = False
        self.filepath = None
        self.override_previous_attempt = override_previous_attempt
        self.download_directory = download_directory
        self.filename = filename
        self.user_agents = [
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3",
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/14.0.3 Safari/605.1.15",
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:88.0) Gecko/20100101 Firefox/88.0",
            "Mozilla/5.0 (Linux; Android 6.0; Nexus 5 Build/MRA58N) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/88.0.4324.182 Mobile Safari/537.36",
        ]
 
    def download(self):
        """
        Finds and downloads a paper using the DOI or PMID provided.
        """
        if not self.override_previous_attempt and self._look_for_previous_download():
            return self
        
        self.check_open_access()
        self.check_pubmed_central_access()
        self.check_crossref_access(decode_doi(self.doi))
        if len(self.pdf_urls) > 0:
            print("[PyPaperRetriever] Found Open-Access PDF link(s). Attempting download...")
            if self._download_pdf():
                return self
        if self.allow_scihub:
            self.pdf_urls = []
            self.check_scihub_access()
            if len(self.pdf_urls) > 0:
                print("[PyPaperRetriever] Found PDF on Sci-Hub. Attempting download...")
                self._download_pdf()
            else:
                print(f"[PyPaperRetriever] No PDFs found for {decode_doi(self.doi)}")
        else:
            print(f"[PyPaperRetriever] No Open-Access PDF found for {decode_doi(self.doi)}. Sci-Hub access is disabled.")
        return self
    
    def check_open_access(self):
        """
        Checks if an article is open access using Unpaywall and updates the instance with access information.
        """

        url = f"https://api.unpaywall.org/v2/{decode_doi(self.doi)}?email={self.email}"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            pdf_urls = [None, None, None, None]
            pdf_locations = [loc.get("url_for_pdf") for loc in data.get("oa_locations", []) if loc.get("url_for_pdf")]
            pdf_urls[:len(pdf_locations)] = pdf_locations[:4]
            pdf_urls = [url for url in pdf_urls if url]

            pubmed_europe_info = next((
                (loc.get("url").split("?")[0], loc.get("url").split("pmc")[-1].split("/")[0])
                for loc in data.get("oa_locations", [])
                if "europepmc.org/articles/pmc" in loc.get("url", "")
            ), (None, None))
            pubmed_europe_url, pmcid = pubmed_europe_info # Not used in current implementation

            if len(pdf_urls) > 0:
                self.is_oa = True
                self.pdf_urls += pdf_urls
            return self
        
        else:
            print("error", f"Unpaywall API request failed with status code {response.status_code}")
            return self
        
    def check_pubmed_central_access(self):
        """
        Checks if an article is available in PubMed Central, and finds associated PDF links.
        """
        pmc_id = None
        id = self.pmid if self.pmid else doi_to_pmid(decode_doi(self.doi), self.email)
        records = entrez_efetch(self.email, id)
        try:
            id_list = records['PubmedArticle'][0]['PubmedData']['ArticleIdList']
            for element in id_list:
                if element.attributes.get('IdType') == 'pmc':
                    pmc_id = str(element)

        except Exception as e:
            print(f"Error processing while checking PMC access for id {id}: {e}")

        if pmc_id is not None:

            article_link = f'https://pmc.ncbi.nlm.nih.gov/articles/{pmc_id}/'

            response = requests.get(article_link, headers={"User-Agent": random.choice(self.user_agents)})

            if response.status_code == 200:
                soup = BeautifulSoup(response.content, 'html.parser')
                # Find PDF links
                pdf_links = [a['href'] for a in soup.find_all('a', href=True) if a['href'].endswith('.pdf')]
                pdf_links = [f"{article_link}{pdf_link}" if pdf_link.startswith('/') else pdf_link for pdf_link in pdf_links]
                pdf_links = list(set(pdf_links))
                for link in pdf_links:
                    self.pdf_urls.append(link)
            else:
                print(f"Failed to fetch the PubMed Central link. Status code: {response.status_code}")

    def check_crossref_access(self, doi):
        """
        Checks if an article is available on Crossref, and finds associated PDF links.
        """
        base_url = "https://api.crossref.org/works/"
        full_url = f"{base_url}{doi}"
        urls = []
        pdf_urls = []
        
        try:
            response = requests.get(full_url)
            if response.status_code == 200:
                data = response.json()
                primary_url = data.get('message', {}).get('URL', None)
                if primary_url:
                    urls.append(primary_url)
                doi_link = f"https://doi.org/{doi}"
                urls.append(doi_link)
            
            for url in urls:
                try:
                    response = requests.get(url, headers={
                        "User-Agent":  "Mozilla/5.0 (Windows NT 10.0; Win64; x64) " +
                                    "AppleWebKit/537.36 (KHTML, like Gecko) " +
                                    "Chrome/58.0.3029.110 Safari/537.3"
                    }, timeout=10)  # Added timeout for better error handling
                    
                    if response.status_code == 200:
                        final_url = response.url  # The final resolved URL after redirects
                        soup = BeautifulSoup(response.content, 'html.parser')
                        
                        pdf_links = set()

                        # 1. Extract PDF links from <a> tags
                        for a in soup.find_all('a', href=True):
                            href = a['href']
                            if href.lower().endswith('.pdf'):
                                absolute_url = urljoin(final_url, href)
                                pdf_links.add(absolute_url)
                        
                        # 2. Extract PDF links from JavaScript
                        for script in soup.find_all('script'):
                            if script.string:
                                # Regex to find patterns like window.open('/path/to/file.pdf') or href = "/path/to/file.pdf"
                                matches = re.findall(r'''(?:window\.open|href\s*=\s*)\(['"]([^'"]+\.pdf)['"]\)''', script.string, re.IGNORECASE)
                                for match in matches:
                                    absolute_url = urljoin(final_url, match)
                                    pdf_links.add(absolute_url)
                                
                                # Another regex pattern based on the example provided
                                matches = re.findall(r'''location\s*=\s*['"]([^'"]+\.pdf)['"]''', script.string, re.IGNORECASE)
                                for match in matches:
                                    absolute_url = urljoin(final_url, match)
                                    pdf_links.add(absolute_url)

                        # 3. Optionally, search for direct links in data attributes or other patterns
                        # Example: data-pdf-url="/path/to/file.pdf"
                        data_pdf_urls = re.findall(r'data-pdf-url=["\']([^"\']+\.pdf)["\']', response.text, re.IGNORECASE)
                        for match in data_pdf_urls:
                            absolute_url = urljoin(final_url, match)
                            pdf_links.add(absolute_url)

                        # Remove any invalid URLs (optional)
                        valid_pdf_links = set()
                        for link in pdf_links:
                            parsed = urlparse(link)
                            if parsed.scheme in ['http', 'https']:
                                valid_pdf_links.add(link)

                        if valid_pdf_links:
                            for link in valid_pdf_links:
                                pdf_urls.append(link)

                    else:
                        print(f"Failed to access URL: {url} with status code {response.status_code}")
                
                except requests.exceptions.RequestException as e:
                    print(f"Error accessing URL: {url}")
                    print(e)

            final_pdf_urls = list(set(pdf_urls))
            for link in final_pdf_urls:
                self.pdf_urls.append(link)
        
        except requests.exceptions.RequestException as e:
            print("Something went wrong while trying to access Crossref API")
            print(e)

    def check_scihub_access(self):
        """
        Checks access to Sci-Hub mirrors for a given DOI and retrieves the PDF URL if available.
        This method iterates over a list of Sci-Hub mirrors, constructs URLs using the DOI,
        and sends HTTP GET requests to check for access. It handles rate limiting by adding
        random delays between requests and uses random user agents to avoid detection.
        If a PDF is found, it updates the instance attributes `on_scihub` and `pdf_urls`.
        """
        mirror_list = ["https://sci-hub.st", "https://sci-hub.ru", "https://sci-hub.se"]
        urls = [f"{mirror}/{decode_doi(self.doi)}" for mirror in mirror_list]

        for i, url in enumerate(urls):
            time.sleep(random.randint(1, 3)) # Delay between requests, avoids being blocked
            headers = {
                "User-Agent": random.choice(self.user_agents),
                "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9",
                "Accept-Language": "en-US,en;q=0.9",
                "Referer": "https://www.google.com/",
            }
            try:
                r = requests.get(url, headers=headers)
                if r.status_code == 200:
                    if len(r.text) < 1:
                        print("""We probably got blocked by Sci-Hub for too many requests. 
                            For being a source of free scientific knowledge, they sure are stingy with their bandwidth.
                            Although they don't specify rate limit and have no robots.txt, they still block IPs with too many requests.
                            Try connecting to a different proxy IP with a VPN.""")
                        break
                    result = self._get_pdf_element(r.text, mirror_list[i])
                    if result == "unavailable":
                        continue
                    elif result:
                        self.on_scihub = True
                        self.pdf_urls.append(result)
                        break
            except requests.RequestException as e:
                print(f"Failed to scrape {url} due to {e}")
                print("If this error includes 'Connection reset by peer', your ISP may be blocking Sci-Hub. Try using a VPN, like ProtonVPN.")
        return self

    def _download_pdf(self):
        """
        Downloads a PDF from the first successful URL found in the list of PDF URLs.
        Stores the PDF in a specified directory, verifies the download, and creates a JSON sidecar.
        """
        file_directory, pdf_path, json_path = self._determine_paths()
        os.makedirs(file_directory, exist_ok=True)

        headers = {
            "User-Agent": random.choice(self.user_agents),
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9",
            "Accept-Language": "en-US,en;q=0.9",
            "Referer": "https://www.google.com/",
        }

        if not self.pdf_urls:
            self.filepath = "unavailable"
            self._create_json_sidecar(download_success=False, pdf_filepath=pdf_path, json_filepath=json_path)
            return False

        for pdf_url in self.pdf_urls:
            try:
                response = requests.get(pdf_url, headers=headers, stream=True)
                if response.status_code == 200:
                    with open(pdf_path, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            f.write(chunk)
                    # Check if the file is downloaded and not corrupted
                    if self._check_if_downloaded(file_directory, '.pdf'):
                        self._create_json_sidecar(download_success=True, pdf_filepath=pdf_path, json_filepath=json_path, url=pdf_url)
                        print(f"[PyPaperRetriever] PDF downloaded successfully to {pdf_path} for {decode_doi(self.doi)} from {pdf_url}")
                        return True
            except requests.RequestException as e:
                continue

        # If no URLs resulted in a successful download
        self.filepath = "unavailable"
        self._create_json_sidecar(download_success=False, pdf_filepath=pdf_path, json_filepath=json_path)
        print(f"[PyPaperRetriever] Failed to download PDF for {decode_doi(self.doi)}")
        return False

    def _create_json_sidecar(self, download_success, pdf_filepath, json_filepath, url=None):
        """
        Creates a JSON sidecar file with information about the downloaded paper.
        Args:
            download_success (bool): Indicates whether the PDF was successfully downloaded.
            pdf_filepath (str): The file path to the downloaded PDF.
            json_filepath (str): The file path where the JSON sidecar will be saved.
            url (str, optional): The URL from which the PDF was downloaded. Defaults to None.
        The JSON sidecar file contains the following fields:
            - doi (str): The DOI of the paper.
            - pmid (str): The PubMed ID of the paper.
            - id (str): The ID of the paper (either the DOI or the PMID).
            - source_url (str): The URL from which the PDF was downloaded.
            - all_urls (list): A list of all URLs attempted for downloading the PDF.
            - download_success (bool): Indicates whether the PDF was successfully downloaded.
            - pdf_filepath (str): The file path to the downloaded PDF, or "unavailable" if the download failed.
            - open_access (bool): Indicates whether the paper is open access.
        """

        open_access = self.is_oa
        if (url and ('scihub' in url or 'sci-hub' in url)) :
            open_access = False
        info = {
            'doi': decode_doi(self.doi),
            'encoded_doi': self.doi,
            'pmid': self.pmid,
            'id': self.pmid if self.pmid else self.doi,
            'source_url':url,
            'all_urls': self.pdf_urls,
            'download_success': download_success,
            'pdf_filepath': pdf_filepath if download_success else "unavailable",
            'open_access': open_access
        }
        with open(json_filepath, 'w') as f:
            json.dump(info, f, indent=4)

    def _determine_paths(self):
        """
        Determines the file directory, PDF path, and JSON path for the download.
        """
        if self.filename:
            file_directory = self.download_directory
            pdf_path = os.path.join(self.download_directory, self.filename)
            json_path = os.path.join(self.download_directory, self.filename.replace('.pdf', '.json'))
        else:
            subdir_name = f"pmid-{self.pmid}" if self.pmid else f"doi-{self.doi}"
            file_directory = os.path.join(self.download_directory, subdir_name)
            filename = f"pmid-{self.pmid}.pdf" if self.pmid else f"doi-{self.doi}.pdf"
            pdf_path = os.path.join(file_directory, filename)
            json_path = os.path.join(file_directory, filename.replace('.pdf', '.json'))

        return file_directory, pdf_path, json_path

    def _check_if_downloaded(self, download_directory_or_path, filetype=".pdf"):
        """Check if a non-corrupted file with the given extension exists in the specified directory or path and delete corrupted files."""
        files_with_type = glob(os.path.join(download_directory_or_path, f'*{filetype}'))

        non_corrupted_files = 0
        
        for file_path in files_with_type:
            try:
                with pymupdf.open(file_path) as doc:
                    if len(doc) > 0:  
                        non_corrupted_files += 1
            except Exception as e:
                os.remove(file_path)
        if non_corrupted_files > 0:
            self.is_downloaded = True
            self.filepath = file_path
        return self
    
    def _look_for_previous_download(self):
        """
        Looks for a previous download of the PDF file in the download directory.
        Uses the paths determined by _determine_paths() to ensure consistency with current settings.
        """
        file_directory, pdf_path, json_path = self._determine_paths()

        if os.path.exists(json_path):
            json_data = json.load(open(json_path))
            self.filepath = json_data.get("pdf_filepath", "unavailable")
            self.is_downloaded = json_data.get("download_success", False)
            return True
        return False

    def _get_pdf_element(self, html_text, mirror):
        """
        Extracts the PDF link from the HTML text of a Sci-Hub page.
        """
        soup = BeautifulSoup(html_text, 'lxml')
        pdf_link = ""
        if soup.find('p', string="Unfortunately, Sci-Hub doesn't have the requested document:"):
            return "unavailable"
        embed_tag = soup.find('embed', {'type': 'application/pdf'})
        if embed_tag and embed_tag.has_attr('src'):
            pdf_link_raw = embed_tag['src']
            if pdf_link_raw.startswith('/downloads') or pdf_link_raw.startswith('/tree') or pdf_link_raw.startswith('/uptodate'):
                pdf_link = f"{mirror}{pdf_link_raw}"
            elif pdf_link_raw.startswith('//'):
                pdf_link = 'https:' + pdf_link_raw
            else:
                pdf_link = pdf_link_raw
        return pdf_link

def main():
    """
    Command-line interface for fetching scientific papers.
    """
    parser = argparse.ArgumentParser(description='Download scientific papers automatically.')
    parser.add_argument('--email', required=True, help='Email address for API usage.')
    parser.add_argument('--doi', help='Digital Object Identifier of the paper.')
    parser.add_argument('--pmid', help='PubMed ID of the paper.')
    parser.add_argument('--dwn-dir', default='PDFs', help='Directory to download the PDFs into. Defaults to "PDFs".')
    parser.add_argument('--filename', help='Custom filename for the downloaded PDF.')
    parser.add_argument('--override', action='store_true', help='Override previous download attempts.')
    parser.add_argument('--allow-scihub', action='store_true', default=True, help='Allow downloading from Sci-Hub if available.')

    args = parser.parse_args()

    retriever = PaperRetriever(
        email=args.email,
        doi=args.doi,
        pmid=args.pmid,
        download_directory=args.dwn_dir,
        filename=args.filename,
        override_previous_attempt=args.override,
        allow_scihub=args.allow_scihub
    )

    retriever.download()