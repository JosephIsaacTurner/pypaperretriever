import os
import requests
import time
import random
from Bio import Entrez
from bs4 import BeautifulSoup
import json
from glob import glob
import fitz
import argparse

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
        user_agents (list): List of user agents for making requests.
    Methods:
        check_open_access(): Checks if the paper is open access using Unpaywall.
        check_scihub_access(): Checks if the paper is available on Sci-Hub.
        download_pdf(): Downloads the PDF from the first available URL.
        pmid_to_doi(pmid): Converts a PMID to a DOI using the Entrez API.
        doi_to_pmid(): Converts a DOI to a PMID using the Entrez API.
        find_and_download(override_previous_attempt=False): Finds and downloads the paper.
        download_from_list(url_list): Downloads PDFs from a list of URLs.
        create_json_sidecar(download_success, pdf_filepath, json_filepath, url=None): Creates a JSON sidecar file with download information.
        _determine_paths(): Determines the file paths for the download.
        _check_if_downloaded(download_directory_or_path, filetype=".pdf"): Checks if a non-corrupted file exists in the specified directory.
        _look_for_previous_download(): Looks for a previous download of the PDF file.
        _process_doi(doi): Processes and encodes a DOI.
        _get_pdf_element(html_text, mirror): Extracts the PDF link from the HTML text of a Sci-Hub page.
    """

    def __init__(self, email, doi=None, pmid=None, allow_scihub=True, download_directory='PDFs', filename=None):
        self.email = email
        if not doi and not pmid:
            raise ValueError("Either a DOI or PMID must be provided")
        if not doi and pmid:
            doi = self.pmid_to_doi(pmid)
        self.doi = self._process_doi(doi)
        self.pmid = pmid
        self.allow_scihub = allow_scihub
        self.is_oa = False
        self.on_scihub = False
        self.pdf_urls = []
        self.is_downloaded = False
        self.filepath = None
        self.download_directory = download_directory
        self.filename = filename
        self.user_agents = [
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3",
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/14.0.3 Safari/605.1.15",
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:88.0) Gecko/20100101 Firefox/88.0",
            "Mozilla/5.0 (Linux; Android 6.0; Nexus 5 Build/MRA58N) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/88.0.4324.182 Mobile Safari/537.36",
        ]
    
    def check_open_access(self):
        """
        Checks if an article is open access using Unpaywall and updates the instance with access information.
        """

        url = f"https://api.unpaywall.org/v2/{self.doi}?email={self.email}"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            oa_status = data.get("oa_status", "unknown") # Not used in current implementation
            is_oa = data.get("is_oa", False) # This is not always accurate (may say OA when no PDF link)
            best_oa_location_url = data.get("best_oa_location", {}).get("url", None) if data.get("best_oa_location") else None # Not used in current implementation
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

    def check_scihub_access(self):
        """
        Checks access to Sci-Hub mirrors for a given DOI and retrieves the PDF URL if available.
        This method iterates over a list of Sci-Hub mirrors, constructs URLs using the DOI,
        and sends HTTP GET requests to check for access. It handles rate limiting by adding
        random delays between requests and uses random user agents to avoid detection.
        If a PDF is found, it updates the instance attributes `on_scihub` and `pdf_urls`.
        """
        mirror_list = ["https://sci-hub.st", "https://sci-hub.ru", "https://sci-hub.se"]
        urls = [f"{mirror}/{self.doi}" for mirror in mirror_list]

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

    def download_pdf(self):
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
            self.create_json_sidecar(download_success=False, pdf_filepath=pdf_path, json_filepath=json_path)
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
                        self.create_json_sidecar(download_success=True, pdf_filepath=pdf_path, json_filepath=json_path, url=pdf_url)
                        return True
            except requests.RequestException as e:
                print(f"Failed to download from {pdf_url} due to: {e}")

        # If no URLs resulted in a successful download
        self.filepath = "unavailable"
        self.create_json_sidecar(download_success=False, pdf_filepath=pdf_path, json_filepath=json_path)
        return False

    def pmid_to_doi(self, pmid):
        """
        Converts a PMID to a DOI using the Entrez API.
        """
        records = entrez_efetch(self.email, pmid)
        try:
            id_list = records['PubmedArticle'][0]['PubmedData']['ArticleIdList']
            for element in id_list:
                if element.attributes.get('IdType') == 'doi':
                    return str(element)
            raise ValueError("No DOI found for this PMID")
        except Exception as e:
            print(f"Error processing PMID {pmid}: {e}")
        return None

    def doi_to_pmid(self):
        """
        Converts a DOI to a PMID using the Entrez API, and if that fails, uses the PMC API.
        """
        records = entrez_efetch(self.email, self.doi)
        try:
            id_list = records['PubmedArticle'][0]['PubmedData']['ArticleIdList']
            for element in id_list:
                if element.attributes.get('IdType') == 'pubmed':
                    if len(str(element)) > 3:
                        self.pmid = str(element)
                        return self            
        except Exception as e:
            print(f"Error converting doi {self.doi} to pmid: {e}")
        try:
            url_base = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?"
            url = f"{url_base}ids={self.doi}&format=json"
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                self.pmid = data.get("records", [{}])[0].get("pmid", None)
                return self
        except Exception as e:
            print(f"Error converting doi {self.doi} to pmid: {e}")
        return self
    
    def find_and_download(self, override_previous_attempt=False):
        """
        Finds and downloads a paper using the DOI or PMID provided.
        """
        if not override_previous_attempt and self._look_for_previous_download():
            return self
        
        self.check_open_access()
        if self.pdf_urls:
            if self.download_pdf():
                return self
        if self.allow_scihub:
            self.pdf_urls = []
            self.check_scihub_access().download_pdf()
        else:
            self.download_pdf()
        return self
    
    def download_from_list(self, url_list):
        """
        Downloads a list of PDF URLs.
        """
        for url in url_list:
            self.pdf_urls.append(url)
        self.download_pdf()
        return self

    def create_json_sidecar(self, download_success, pdf_filepath, json_filepath, url=None):
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
                with fitz.open(file_path) as doc:
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

    def _process_doi(self, doi):
        """
        Encodes a DOI to remove any URL parameters and return a filename-safe encoded DOI.
        """
        clean_doi = doi.split("doi.org/")[-1].split("?")[0]
        return encode_doi(clean_doi)

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

def encode_doi(doi):
    """Encodes a DOI for safe inclusion in a filename using URL encoding."""
    encoded_doi = doi.replace('/', '%2F').replace('-', '%2D').replace('.', '%2E')
    return encoded_doi.strip("'\"")

def decode_doi(encoded_doi):
    """Decodes a previously URL-encoded DOI."""
    decoded_doi = encoded_doi.replace('%2F', '/').replace('%2D', '-').replace('%2E', '.')
    return decoded_doi

def entrez_efetch(email, id):
    """
    Fetches a record from the Entrez API.
    """
    Entrez.email = email
    handle = Entrez.efetch(db="pubmed", id=id, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

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
    parser.add_argument('--allow-scihub', action='store_true', default=True, help='Allow downloading from Sci-Hub if available.')

    args = parser.parse_args()

    retriever = PaperRetriever(
        email=args.email,
        doi=args.doi,
        pmid=args.pmid,
        download_directory=args.dwn_dir,
        filename=args.filename,
        allow_scihub=args.allow_scihub
    )

    retriever.find_and_download()

if __name__ == '__main__':
    main()
