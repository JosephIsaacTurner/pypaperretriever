from unittest.mock import patch, Mock
from pypaperretriever.paper_retriever import PaperRetriever, encode_doi
import json
from pathlib import Path
from io import BytesIO
from Bio import Entrez

# Sample DOI and PMID for testing
TEST_DOI = "10.7759/cureus.76081"
TEST_PMID = "12345678"
TEST_EMAIL = "bob_tester@testemail.com"

# For referencing your actual PDF from tests/data/
TESTS_DIR = Path(__file__).parent  # e.g. /Users/jt041/repos/pypaperretriever/tests
PDF_PATH = TESTS_DIR / "data" / f"doi-{encode_doi(TEST_DOI)}" / f"doi-{encode_doi(TEST_DOI)}.pdf"

def mock_requests_get(url, *args, **kwargs):
    """
    Mocks requests.get for:
    - Unpaywall (returns mock JSON)
    - .pdf (returns bytes from the real PDF in tests/data/)
    - Sci-Hub (returns dummy HTML)
    - Everything else (404).
    """
    if url.startswith("https://api.unpaywall.org"):
        # Mock response for Unpaywall API
        return Mock(
            status_code=200, 
            json=lambda: {
                "is_oa": True,
                "best_oa_location": {"url": "https://example.com/sample.pdf"},
                "oa_locations": [{"url_for_pdf": "https://example.com/sample.pdf"}]
            }
        )
    elif url.endswith(".pdf"):
        # Return the real PDF bytes from tests/data
        mock_resp = Mock(status_code=200)
        if PDF_PATH.exists():
            with open(PDF_PATH, "rb") as pdf_file:
                pdf_bytes = pdf_file.read()
        else:
            print("PDF fixture not found at:", PDF_PATH)
            pdf_bytes = b""  # fallback if the file doesn't exist
        
        # Serve the file's content in one chunk
        mock_resp.iter_content = lambda chunk_size: [pdf_bytes]
        return mock_resp
    
    elif url.startswith("https://sci-hub"):
        # Mock response for Sci-Hub
        return Mock(status_code=200, text="<html></html>")
    else:
        # Default mock response for other URLs
        return Mock(status_code=404, text="<html></html>")

def mock_entrez_efetch_doi(*args, **kwargs):
    """
    This mock should NOT be called when a DOI is provided.
    """
    raise AssertionError("Entrez.efetch should not be called when DOI is provided")

def mock_entrez_efetch_pmid(*args, **kwargs):
    """
    Mocks Entrez.efetch for PMID -> DOI conversion by returning a realistic Entrez response.
    """
    path = TESTS_DIR / "data" / "entrez_record.xml"
    data = open(path, 'r').read()
    data = BytesIO(data.encode())
    return Entrez.read(data)

@patch('pypaperretriever.paper_retriever.requests.get', side_effect=mock_requests_get)
@patch('pypaperretriever.paper_retriever.entrez_efetch', side_effect=mock_entrez_efetch_doi)
def test_fetch_paper_with_doi(mock_efetch, mock_get, tmp_path):
    """
    Test PaperRetriever's ability to fetch and download a paper using a DOI.
    Ensures that Entrez.efetch is NOT called when DOI is provided.
    """
    # Initialize PaperRetriever with test DOI and email
    retriever = PaperRetriever(email=TEST_EMAIL, doi=TEST_DOI, download_directory=str(tmp_path))

    # Invoke the method to find and download the paper
    retriever.find_and_download()

    # Assertions to ensure Entrez.efetch was NOT called
    mock_efetch.assert_not_called()

    # Ensure that requests.get was called at least twice:
    # 1. Unpaywall API call
    # 2. PDF download
    assert mock_get.call_count >= 2

    # Check if PDF was "downloaded" to the tmp_path
    encoded = encode_doi(TEST_DOI)  # e.g. "10.7759%2Fcureus%2E76081"
    expected_subdir = tmp_path / f"doi-{encoded}"
    expected_pdf_path = expected_subdir / f"doi-{encoded}.pdf"
    assert expected_pdf_path.exists(), f"Expected PDF at {expected_pdf_path} not found."

    # Verify the PDF content matches what we stored in tests/data
    with open(PDF_PATH, "rb") as source_file, open(expected_pdf_path, "rb") as dest_file:
        source_content = source_file.read()
        dest_content = dest_file.read()
        assert dest_content == source_content, "Downloaded PDF doesn't match the fixture."

    # Check JSON sidecar
    expected_json_path = expected_pdf_path.with_suffix(".json")
    assert expected_json_path.exists(), f"JSON sidecar at {expected_json_path} not found."
    with open(expected_json_path, 'r') as f:
        json_data = json.load(f)
        assert json_data['doi'] == TEST_DOI, f"Expected DOI '{TEST_DOI}' but got '{json_data['doi']}'"
        assert json_data['download_success'] is True
        assert json_data['pdf_filepath'].endswith(str(expected_pdf_path).replace(str(tmp_path), ""))
        assert json_data['open_access'] is True

@patch('pypaperretriever.paper_retriever.requests.get', side_effect=mock_requests_get)
@patch('pypaperretriever.paper_retriever.entrez_efetch', side_effect=mock_entrez_efetch_pmid)
def test_fetch_paper_with_pmid(mock_efetch, mock_get, tmp_path):
    """
    Test PaperRetriever's ability to fetch and download a paper using a PMID.
    Ensures that Entrez.efetch is called when PMID is provided without a DOI.
    """
    # Initialize PaperRetriever with test PMID and email (DOI is None)
    retriever = PaperRetriever(email=TEST_EMAIL, pmid=TEST_PMID, download_directory=str(tmp_path))

    # Invoke the method to find and download the paper
    retriever.find_and_download()

    # Assertions to ensure Entrez.efetch was called
    mock_efetch.assert_called()

    # Ensure that requests.get was called at least twice:
    # 1. Unpaywall API call
    # 2. PDF download
    assert mock_get.call_count >= 2

    # Check if PDF was "downloaded" to the tmp_path
    expected_subdir = tmp_path / f"pmid-{TEST_PMID}"
    expected_pdf_path = expected_subdir / f"pmid-{TEST_PMID}.pdf"
    assert expected_pdf_path.exists(), f"Expected PDF at {expected_pdf_path} not found."

    # Verify that the downloaded PDF matches our fixture PDF
    with open(PDF_PATH, "rb") as source_file, open(expected_pdf_path, "rb") as dest_file:
        source_content = source_file.read()
        dest_content = dest_file.read()
        assert dest_content == source_content, "Downloaded PDF doesn't match the fixture."

    # Check JSON sidecar
    expected_json_path = expected_pdf_path.with_suffix(".json")
    assert expected_json_path.exists(), f"JSON sidecar at {expected_json_path} not found."
    with open(expected_json_path, 'r') as f:
        json_data = json.load(f)
        assert json_data['doi'] == TEST_DOI, f"Expected DOI '{TEST_DOI}' but got '{json_data['doi']}'"
        assert json_data['download_success'] is True
        assert json_data['pdf_filepath'].endswith(str(expected_pdf_path).replace(str(tmp_path), ""))
        assert json_data['open_access'] is True