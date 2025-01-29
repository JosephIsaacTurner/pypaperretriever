import pytest
import pandas as pd
from unittest.mock import patch
from io import BytesIO
from Bio import Entrez
from pathlib import Path
from pypaperretriever import PubMedSearcher

# Test config
TEST_EMAIL = "bob_tester@testemail.com"
TEST_SEARCH_STRING = "brain lesions"

# Paths to stored XML fixtures
TESTS_DIR = Path(__file__).parent
ESEARCH_XML_PATH = TESTS_DIR / "data" / "entrez_esearch.xml"
EFETCH_XML_PATH = TESTS_DIR / "data" / "entrez_efetch.xml"

def load_mock_xml(filepath):
    """Helper function to load test XML files."""
    with open(filepath, "r", encoding="utf-8") as f:
        return BytesIO(f.read().encode())  # Convert to BytesIO for Entrez.read()

@patch("Bio.Entrez.esearch")
@patch("Bio.Entrez.efetch")
def test_pubmed_search(mock_efetch, mock_esearch):
    """
    Tests PubMedSearcher.search() to ensure:
    - It correctly calls the Entrez API.
    - The DataFrame updates with search results.
    """

    # ✅ Use stored test XML files
    mock_esearch.return_value = load_mock_xml(ESEARCH_XML_PATH)
    mock_efetch.return_value = load_mock_xml(EFETCH_XML_PATH)

    # ✅ Initialize PubMedSearcher
    searcher = PubMedSearcher(search_string=TEST_SEARCH_STRING, email=TEST_EMAIL)

    # ✅ Call the search method
    searcher.search(count=2)

    # ✅ Ensure esearch was called
    mock_esearch.assert_called()

    # ✅ Ensure efetch was called
    mock_efetch.assert_called()

    # ✅ Ensure the DataFrame was updated (even though data is mocked)
    assert not searcher.df.empty, "DataFrame should not be empty after a search."

    # ✅ Check if expected article data is in the DataFrame
    assert "Damaged Relay Station: EEG Neurofeedback Training in Isolated Bilateral Paramedian Thalamic Infarct." in searcher.df["title"].values, "Expected 'Damaged Relay Station: EEG Neurofeedback Training in Isolated Bilateral Paramedian Thalamic Infarct.' in DataFrame."