from Bio import Entrez
import requests

def entrez_efetch(email, id):
    """
    Fetches a record from the Entrez API.
    """
    Entrez.email = email
    handle = Entrez.efetch(db="pubmed", id=id, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

def pmid_to_doi(pmid, email):
    """
    Converts a PMID to a DOI using the Entrez API.
    """
    records = entrez_efetch(email, pmid)
    try:
        id_list = records['PubmedArticle'][0]['PubmedData']['ArticleIdList']
        for element in id_list:
            if element.attributes.get('IdType') == 'doi':
                return str(element)
        raise ValueError("No DOI found for this PMID")
    except Exception as e:
        print(f"Error processing PMID {pmid}: {e}")
    return None

def doi_to_pmid(doi, email):
    """
    Converts a DOI to a PMID using the Entrez API, and if that fails, uses the PMC API.
    """
    print(f"[ReferenceRetriever] Converting DOI {doi} to PMID.")
    Entrez.email = email
    try:
        search_handle = Entrez.esearch(db="pubmed", term=doi)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        if search_results['IdList']:
            pmid = search_results['IdList'][0]
            print(f"[ReferenceRetriever] Converted DOI {doi} to PMID: {pmid}")
            return pmid
        else:
            print(f"[ReferenceRetriever] No PMID found for DOI {doi} in Entrez PubMed.")
    except Exception as e:
        print(f"[ReferenceRetriever] Error converting DOI {doi} to PMID via Entrez: {e}")

    try:
        url_base = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?"
        url = f"{url_base}ids={doi}&format=json"
        print(f"[ReferenceRetriever] Requesting PMID via PMC ID Converter from URL: {url}")
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            pmid = data.get("records", [{}])[0].get("pmid", None)
            if pmid:
                print(f"[ReferenceRetriever] Converted DOI {doi} to PMID via PMC: {pmid}")
                return pmid
            else:
                print(f"[ReferenceRetriever] PMC ID Converter did not return a PMID for DOI {doi}.")
        else:
            print(f"[ReferenceRetriever] PMC ID Converter request failed with status code: {response.status_code}")
    except Exception as e:
        print(f"[ReferenceRetriever] Error converting DOI {doi} to PMID via PMC ID Converter: {e}")
    return None

def encode_doi(doi):
    """Encodes a DOI for safe inclusion in a filename using URL encoding."""
    doi = doi.split("doi.org/")[-1].split("?")[0]
    encoded_doi = doi.replace('/', '%2F').replace('-', '%2D').replace('.', '%2E')
    return encoded_doi.strip("'\"")

def decode_doi(encoded_doi):
    """Decodes a previously URL-encoded DOI."""
    decoded_doi = encoded_doi.replace('%2F', '/').replace('%2D', '-').replace('%2E', '.')
    return decoded_doi