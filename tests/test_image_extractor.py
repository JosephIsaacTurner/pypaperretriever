import os
import shutil
import json
import pytest
from pathlib import Path
from pypaperretriever import encode_doi
from pypaperretriever import ImageExtractor  

TEST_DOI = "10.7759/cureus.76081"
TESTS_DIR = Path(__file__).parent
PDF_PATH = TESTS_DIR / "data" / f"doi-{encode_doi(TEST_DOI)}" / f"doi-{encode_doi(TEST_DOI)}.pdf"

@pytest.mark.parametrize("pdf_path", [PDF_PATH])
def test_image_extractor(pdf_path, tmp_path):
    """
    Test the ImageExtractor's ability to extract figures from a known PDF.
    This uses the same test PDF file used in other tests. 
    """
    # Dude, ensure that our fixture PDF actually exists
    assert pdf_path.exists(), f"Test PDF not found at: {pdf_path}"

    # Copy the PDF to our pytest tmp_path
    test_pdf_path = tmp_path / "test.pdf"
    shutil.copy(pdf_path, test_pdf_path)

    # Initialize the ImageExtractor with our test PDF
    extractor = ImageExtractor(str(test_pdf_path))

    # Run the extraction
    extractor.extract_images()

    # If the PDF is valid, we expect either the native PDF or image-based logic to run
    if extractor.is_valid_pdf:
        # Check that we have at least one extracted image, or gracefully handle the case of zero
        # (depending on how many images exist in your actual test PDF).
        # For demonstration, we'll just check the final paths list.
        assert len(extractor.img_paths) > 0, "No images were extracted from the PDF."

        # Verify that each extracted image file and its JSON sidecar exist
        for img_file in extractor.img_paths:
            assert os.path.exists(img_file), f"Extracted image not found: {img_file}"
            sidecar_json = Path(img_file).with_suffix(".json")
            assert sidecar_json.exists(), f"No JSON sidecar found for {img_file}"
            
            # Optionally, read and inspect JSON content
            with open(sidecar_json, "r") as f:
                data = json.load(f)
                # Expect some fields like 'img_number'
                assert "img_number" in data, f"JSON sidecar for {img_file} missing 'img_number'."
    else:
        # If it's not a valid PDF, check that no images were extracted
        assert len(extractor.img_paths) == 0, (
            "Images were extracted from an invalid PDF, which shouldn't happen."
        )
