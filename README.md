## PyPaperRetriever

<img src="logo.png" width="200">

**A python package for retrieving scientific papers from the web.** Inspired by PyPaperBot (https://github.com/ferru97/PyPaperBot) but with improved flexibility and extensibility. Prefers open-access sources but users can opt to use Sci-Hub as a fallback depending on their ethical considerations and local laws. 

### Installation
    
```bash
pip install git+https://github.com/josephisaacturner/pypaperretriever.git
```

### Usage
#### A. Pythonic OOP approach

```python
from pypaperretriever import PyPaperRetriever

doi = "10.1056/NEJMra1706158"
email = "your_email@gmail.com"
download_dir = "pdf_downloads"
allow_scihub = True # If False, will only use open-access sources
filename = "fox_nejfm_2018.pdf" # Optional, defaults to doi-<doi>/doi-<doi>.pdf for interoperability with PyBIDS

retriever = PyPaperRetriever(email=email, doi=doi, download_directory=download_dir, allow_scihub=allow_scihub)
result = retriever.find_and_download()
if result.is_downloaded:
    print("Downloaded to", result.filepath)
```

#### B. Command-line approach

```bash 
python -m pypaperretriever --email your_email@gmail.com --doi 10.1056/NEJMra1706158 --dwn-dir pdf_downloads --allow-scihub --filename fox_nejfm_2018.pdf
```

### Features

- Retrieve papers by DOI or PMID
- Finds all available sources from Unpaywall and Sci-Hub
- Keeps track of sources used for each download by saving a json sidecar file for each download attempt
- Won't download the same paper twice if it's already been downloaded
- Specify your own filepath, or use our default naming convention that follows BIDS standards (https://bids-specification.readthedocs.io/en/stable/)
- Command-line interface for easy use in scripts
- Pythonic class-based for flexibility and extensibility
- Improves upon PyPaperBot by allowing open-access sources and improves web-scraping robustness to find more papers that are overlooked by PyPaperBot