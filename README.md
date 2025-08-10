# PyPaperRetriever

<img src="logo.png" width="200">

[![Docs](https://img.shields.io/badge/Docs-Read%20the%20docs-blue?logo=readthedocs)](https://josephiturner.com/pypaperretriever/)
[![Quickstart](https://img.shields.io/badge/Quickstart-5%20min-brightgreen)](https://josephiturner.com/pypaperretriever/quickstart/)

👉 **[Full Documentation](https://josephiturner.com/pypaperretriever/)**  
👉 **[Start with the Quickstart](https://josephiturner.com/pypaperretriever/quickstart/)**

**A python package for retrieving scientific papers from the web.** Inspired by PyPaperBot (https://github.com/ferru97/PyPaperBot) but with improved flexibility and extensibility. Prefers open-access sources but users can opt to use Sci-Hub as a fallback depending on their ethical considerations and local laws.

### Installation
    
```bash
pip install git+https://github.com/josephisaacturner/pypaperretriever.git
```

### Features

- Download papers using DOI or PubMed ID (PMID)
- Search PubMed programmatically with advanced query options
- Track citation networks (both upstream and downstream) for papers of interest
- Extract images from downloaded PDFs
- Find all available sources from Unpaywall and optional Sci-Hub integration
- Keep track of sources used via JSON sidecar files for each download
- Avoid duplicate downloads with intelligent checking
- BIDS-compatible file naming convention
- Both command-line and Python API interfaces
- Advanced search capabilities with customizable filters
- Citation network analysis tools

### Ethical and legal note on Sci-Hub
Use of Sci-Hub is disabled by default and clearly labeled. Institutions and researchers differ in policy and legal context; PyPaperRetriever exposes an opt-in flag so users can comply with local rules while retaining a complete pipeline for contexts where such access is permitted. The authors of PyPaperRetriever do not endorse or encourage the use of Sci-Hub in violation of local laws or institutional policies. Users are responsible for ensuring compliance with all applicable laws and ethical guidelines when using this tool.

### Usage Examples

For complete examples, see [examples.ipynb](examples.ipynb) in the repository.

#### 1. Download Using DOI

```python
from pypaperretriever import PaperRetriever

retriever = PaperRetriever(
    email="your.email@gmail.com",
    doi="10.7759/cureus.76081",
    download_directory='PDFs'
)
retriever.download()

# Command-line alternative
pypaperretriever --doi 10.7759/cureus.76081 --email your.email@gmail.com --dwn-dir PDFs
```

#### 2. Download Using PubMed ID

```python
from pypaperretriever import PaperRetriever

retriever = PaperRetriever(
    email="your.email@gmail.com",
    pmid="33813262",
    download_directory='PDFs'
)
retriever.download()

# Command-line alternative
pypaperretriever --pmid 33813262 --email your.email@gmail.com --dwn-dir PDFs
```

#### 3. Control Sci-Hub Access

```python
retriever = PaperRetriever(
    email="your.email@gmail.com",
    doi="10.1016/j.revmed.2011.10.009",
    download_directory='PDFs',
    allow_scihub=False  # Set to True to enable Sci-Hub
)
retriever.download()
```

#### 4. Extract Images from PDFs

```python
from pypaperretriever import ImageExtractor

extractor = ImageExtractor('path/to/your/paper.pdf')
extractor.extract_images()
```

#### 5. Search PubMed Programmatically

```python
from pypaperretriever import PubMedSearcher

search_query = """("brain lesions"[MeSH Terms] OR "brain lesion"[Title/Abstract] OR 
                   "cerebral lesion"[Title/Abstract]) AND (case reports[Publication Type])"""

searcher = PubMedSearcher(search_string=search_query, email="your.email@gmail.com")

results = searcher.search(
    count=10,
    order_by='relevance',  # or 'chronological'
    only_open_access=False,
    only_case_reports=False
)

# Download found articles
searcher.download_articles(download_directory='PDFs', allow_scihub=True)

# Extract images from downloaded articles
searcher.extract_images()
```

#### 6. Track Citation Networks

```python
from pypaperretriever import PaperTracker

tracker = PaperTracker(
    email="your.email@gmail.com",
    doi='10.1097/RLU.0000000000001894',
    max_upstream_generations=1,   # Papers referenced by your paper
    max_downstream_generations=1  # Papers that cite your paper
)

results = tracker.track_paper()
```

### Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### License

MIT License

### Citation

If you use PyPaperRetriever in your research, please cite:

[Add citation information here later]