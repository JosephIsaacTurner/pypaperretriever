from Bio import Entrez
import requests


def entrez_efetch(email, id):
    """Fetch a record from the Entrez API.

    Args:
        email (str): Email address required by the Entrez API.
        id (str): PubMed identifier to fetch.

    Returns:
        Any: Parsed record returned by the Entrez API.
    """
    Entrez.email = email
    handle = Entrez.efetch(db="pubmed", id=id, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records


def pmid_to_doi(pmid, email):
    """Convert a PMID to a DOI using the Entrez API.

    Args:
        pmid (str): PubMed identifier to convert.
        email (str): Email address required by the Entrez API.

    Returns:
        str | None: DOI if found, otherwise ``None``.

    Raises:
        ValueError: If no DOI is associated with the PMID.
    """
    records = entrez_efetch(email, pmid)
    try:
        id_list = records['PubmedArticle'][0]['PubmedData']['ArticleIdList']
        for element in id_list:
            if element.attributes.get('IdType') == 'doi':
                return str(element)
        raise ValueError("No DOI found for this PMID")
    except Exception as e:  # pragma: no cover - informative print
        print(f"Error processing PMID {pmid}: {e}")
    return None


def doi_to_pmid(doi, email):
    """Convert a DOI to a PMID.

    The function first queries the Entrez API. If that fails, the PMC ID converter
    service is used as a fallback.

    Args:
        doi (str): Digital Object Identifier to convert.
        email (str): Email address required by the Entrez API.

    Returns:
        str | None: PMID if found, otherwise ``None``.
    """
    Entrez.email = email
    try:
        search_handle = Entrez.esearch(db="pubmed", term=doi)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        if search_results['IdList']:
            pmid = search_results['IdList'][0]
            return pmid
    except Exception as e:  # pragma: no cover - informative print
        print(f"[ReferenceRetriever] Error converting DOI {doi} to PMID via Entrez: {e}")

    try:
        url_base = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?"
        url = f"{url_base}ids={doi}&format=json"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            pmid = data.get("records", [{}])[0].get("pmid", None)
            if pmid:
                return pmid
        else:  # pragma: no cover - network behaviour
            print(
                f"[PyPaperRetriever] PMC ID Converter request failed with status code: {response.status_code}"
            )
    except Exception as e:  # pragma: no cover - informative print
        print(f"[PyPaperRetriever] Error converting DOI {doi} to PMID via PMC ID Converter: {e}")
    return None


def encode_doi(doi):
    """Encode a DOI for safe inclusion in file names.

    Args:
        doi (str): DOI to encode. It may include a ``doi.org`` prefix.

    Returns:
        str: URL-encoded DOI suitable for use in file paths.
    """
    doi = doi.split("doi.org/")[-1].split("?")[0]
    encoded_doi = doi.replace('/', '%2F').replace('-', '%2D').replace('.', '%2E')
    return encoded_doi.strip("'\"")


def decode_doi(encoded_doi):
    """Decode a previously encoded DOI.

    Args:
        encoded_doi (str): DOI encoded with :func:`encode_doi`.

    Returns:
        str: The decoded DOI.
    """
    decoded_doi = encoded_doi.replace('%2F', '/').replace('%2D', '-').replace('%2E', '.')
    return decoded_doi