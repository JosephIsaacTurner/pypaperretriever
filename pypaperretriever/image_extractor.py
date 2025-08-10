"""Utilities for extracting figures from PDF documents."""

import os
import fitz as pymupdf
import numpy as np
import cv2
from pdf2image import convert_from_path
import json
from PIL import Image, ImageOps
import io


class ImageExtractor:
    """Extract figures from a PDF file.

    The extractor handles both native PDFs (containing embedded images) and
    scanned PDFs where each page is an image.

    Args:
        pdf_file_path (str): Path to the PDF file to process.

    Attributes:
        filepath (str): Path to the PDF file.
        dir (str): Directory containing the PDF file.
        is_valid_pdf (bool): Whether the file can be opened by PyMuPDF.
        is_native_pdf (bool): ``True`` if the PDF contains embedded text.
        img_paths (list[str]): Paths to extracted image files.
        img_counter (int): Counter used to name extracted images.
        id (str | None): Optional identifier prefix for saved images.
    """

    def __init__(self, pdf_file_path):
        """Initialize the extractor.

        Args:
            pdf_file_path (str): Path to the PDF file to process.
        """
        self.filepath = pdf_file_path
        self.dir = os.path.dirname(pdf_file_path)
        self.is_valid_pdf = False
        self.is_native_pdf = False
        self.img_paths = []
        self.img_counter = 0  # Initialize shared image counter
        self.id = None
        self._determine_if_valid_pdf()  # Sometimes PDFs are corrupted and cannot be opened
        if self.is_valid_pdf:
            self._check_pdf_type()
            self._get_metadata()

    def extract_images(self):
        """Extract images from the PDF.

        The method determines the PDF type and delegates to the appropriate
        extraction routine. Extracted image paths are stored in ``img_paths``.

        Returns:
            ImageExtractor: The current instance with ``img_paths`` populated.
        """
        if not self.is_valid_pdf:
            print("PDF is not valid.")
            return self
        if self.is_native_pdf:
            self.extract_from_native_pdf()
            self.handle_image_based_pdf()
        else:
            self.handle_image_based_pdf()

    def extract_from_native_pdf(self):
        """Extract figures from a native PDF using PyMuPDF.

        Saves each valid image to disk and records its file path.
        """
        try:
            with pymupdf.open(self.filepath) as doc:
                for page_num in range(len(doc)):
                    page = doc.load_page(page_num)
                    for img in page.get_images(full=True):
                        xref = img[0]
                        base_image = doc.extract_image(xref)
                        image_bytes = base_image["image"]
                        color_space = base_image.get("cs")  # Get color space
                        bpc = base_image.get("bpc", 8)     # Get bits per component, default to 8

                        # Use PIL to handle color space conversions
                        image_pil = Image.open(io.BytesIO(image_bytes))

                        # Handle different color spaces
                        if color_space == "DeviceCMYK":
                            image_pil = image_pil.convert("CMYK").convert("RGB")
                        elif color_space == "DeviceGray":
                            image_pil = image_pil.convert("L").convert("RGB")
                        elif color_space == "DeviceRGB":
                            image_pil = image_pil.convert("RGB")
                        else:
                            # Handle other or unknown color spaces if necessary
                            image_pil = image_pil.convert("RGB")

                        # Handle transparency
                        if image_pil.mode in ("RGBA", "LA") or (image_pil.mode == "P" and 'transparency' in image_pil.info):
                            # Create a white background image
                            background = Image.new("RGB", image_pil.size, (255, 255, 255))
                            background.paste(image_pil, mask=image_pil.split()[-1])  # Paste with alpha channel as mask
                            image_pil = background

                        # **Handle BPC Inversion**
                        if bpc == 1 and image_pil.mode == 'L':
                            print(f"Inverting image {self.img_counter} on page {page_num} due to bpc=1 and grayscale mode.")
                            # Invert the image to correct negative appearance
                            image_pil = ImageOps.invert(image_pil)
                            # Optionally, convert back to RGB
                            image_pil = image_pil.convert("RGB")

                        # Convert PIL Image to NumPy array for validation
                        image_np = np.array(image_pil)

                        if self._check_valid_img(image_np):
                            id_prefix = f"id-{self.id}_" if self.id else ""
                            img_filepath = os.path.join(self.dir, "images", f"{id_prefix}img-{self.img_counter}.png")
                            os.makedirs(os.path.dirname(img_filepath), exist_ok=True)
                            
                            # Save the image using PIL to ensure correct color space and handling
                            image_pil.save(img_filepath, "PNG")
                            
                            self._make_json_sidecar(self.img_counter)
                            self.img_paths.append(img_filepath)
                            self.img_counter += 1
        except Exception as e:
            print(f"Error extracting from native PDF: {e}")

    def handle_image_based_pdf(self):
        """Process a scanned PDF.

        Each page is converted to an image and potential figures are extracted
        using :meth:`_crop_boxes_in_image`.
        """
        try:
            pages = convert_from_path(self.filepath, 300)  # DPI set to 300 for good quality

            for page_num, page in enumerate(pages):
                img_filepath = os.path.join(self.dir, "images", f"page_{page_num}.png")
                os.makedirs(os.path.dirname(img_filepath), exist_ok=True)
                page.save(img_filepath, 'PNG')
                self.img_counter = self._crop_boxes_in_image(img_filepath)
                os.remove(img_filepath)

        except Exception as e:
            print(f"Error handling image-based PDF: {e}")

    def _crop_boxes_in_image(self, img_path):
        """Crop figure-like regions from an image.

        Args:
            img_path (str): Path to the page image to analyse.

        Returns:
            int: Updated image counter after cropping operations.
        """
        img = cv2.imread(img_path)
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        _, thresh = cv2.threshold(gray, 150, 255, cv2.THRESH_BINARY_INV)
        kernel = np.ones((5,5),np.uint8)
        closing = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE, kernel)
        contours, _ = cv2.findContours(closing, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        for contour in contours:
            if cv2.contourArea(contour) < 1000:
                continue
            x, y, w, h = cv2.boundingRect(contour)
            if w / float(h) < 0.1 or w / float(h) > 10 or w < 100 or h < 100:
                continue
            cropped_img = img[y:y+h, x:x+w]
            if self._check_valid_img(cropped_img):  # Check if the cropped image is valid
                id_prefix = f"id-{self.id}_" if self.id else ""
                img_filepath = os.path.join(self.dir, "images", f"{id_prefix}img-{self.img_counter}.png")
                cv2.imwrite(img_filepath, cropped_img)
                self._make_json_sidecar(self.img_counter)
                self.img_paths.append(img_filepath)
                self.img_counter += 1

        return self.img_counter

    def _get_metadata(self):
        """Load companion metadata if a JSON sidecar exists."""
        metadata_json = self.filepath.replace(".pdf", ".json")
        if os.path.exists(metadata_json):
            try:
                with open(metadata_json, 'r', encoding='utf-8') as f:
                    metadata = json.load(f)
                    for key, value in metadata.items():
                        setattr(self, key, value)
            except Exception as e:
                print(f"Error reading metadata JSON: {e}")

    def _check_pdf_type(self):
        """Determine whether the PDF is native or image based.

        Checks the first few pages for the presence of text; absence suggests a
        scanned PDF.
        """
        try:
            with pymupdf.open(self.filepath) as doc:
                for page_num in range(min(5, len(doc))):  # Check the first 5 pages
                    page = doc.load_page(page_num)
                    text = page.get_text()
                    if len(text) > 50:  # Assuming more than 50 characters indicates native
                        self.is_native_pdf = True
                        break  # Found substantial text, no need to check further
                else:
                    self.is_native_pdf = False  # No substantial text found in first 5 pages
        except Exception as e:
            print(f"Error checking PDF type: {e}")
            self.is_native_pdf = False

    def _check_valid_img(self, img):
        """Validate whether an image likely represents a figure.

        Args:
            img (numpy.ndarray | PIL.Image.Image): Image to validate.

        Returns:
            bool: ``True`` if the image appears to be a valid figure.
        """
        if isinstance(img, np.ndarray):
            h, w = img.shape[:2]
        else:
            # If img is a PIL Image
            w, h = img.size

        aspect_ratio = w / float(h)
        min_width, min_height = 100, 100  # Minimum dimensions
        min_area = 1000  # Minimum pixel area
        max_aspect_ratio = 10  # Maximum aspect ratio
        if not (w >= min_width and h >= min_height and w * h >= min_area and
                0.1 <= aspect_ratio <= max_aspect_ratio):
            return False  # Image fails basic dimension checks

        # Check for unique pixel intensities and variance
        if isinstance(img, np.ndarray):
            grayscale_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY) if len(img.shape) == 3 else img
        else:
            grayscale_img = img.convert("L")
            grayscale_img = np.array(grayscale_img)
        
        unique_intensities = len(np.unique(grayscale_img))
        variance = np.var(grayscale_img)

        min_unique_intensities = 100  # Expect at least 100 unique intensities for complexity
        min_variance = 800  # Set based on expected variance in data

        if unique_intensities >= min_unique_intensities and variance >= min_variance:
            return True  # Image passes all checks
        return False  # Image is likely not neuroimaging data or is too simple

    def _make_json_sidecar(self, img_counter):
        """Create a JSON sidecar for an extracted image.

        Args:
            img_counter (int): Sequential number of the image being saved.
        """
        id_prefix = f"id-{self.id}_" if self.id else ""
        json_filepath = os.path.join(self.dir, "images", f"{id_prefix}img-{img_counter}.json")
        os.makedirs(os.path.dirname(json_filepath), exist_ok=True)
        json_content = {
            "img_number": img_counter
        }
        for key, value in vars(self).items():
            if key not in ["dir", "is_valid_pdf", "is_native_pdf", "img_paths", "img_counter"]:
                json_content[key] = value

        try:
            with open(json_filepath, "w", encoding='utf-8') as f:
                json.dump(json_content, f, indent=4)
        except Exception as e:
            print(f"Error writing JSON sidecar for image {img_counter}: {e}")

    def _determine_if_valid_pdf(self):
        """Check whether the PDF can be opened by PyMuPDF."""
        try:
            with pymupdf.open(self.filepath) as doc:
                self.is_valid_pdf = True
        except Exception as e:
            print(f"Error opening PDF: {e}")
            self.is_valid_pdf = False
