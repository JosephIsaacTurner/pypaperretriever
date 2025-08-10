# Quickstart

This guide highlights common workflows of **PyPaperRetriever**.

## Installation

```bash
pip install git+https://github.com/josephisaacturner/pypaperretriever.git
```

## Download a paper

Download a PDF by DOI:

```python
from pypaperretriever import PaperRetriever

retriever = PaperRetriever(
    email="your.email@gmail.com",
    doi="10.7759/cureus.76081",
    download_directory="PDFs",
    allow_scihub=True # If False, will only use oepn access sources
)
retriever.download()
```

If `allow_scihub=True`, PyPaperRetriever will fall back to Sci-Hub if no open access copy is found. Sci-Hub is not used by default due to legal and ethical considerations.

To download via PubMed ID:

```python
from pypaperretriever import PaperRetriever

retriever = PaperRetriever(
    email="your.email@gmail.com",
    pmid="33813262",
    download_directory="PDFs",
    allow_scihub=True
)
retriever.download()
```

Or by CLI:

```bash
pypaperretriever --doi 10.7759/cureus.76081 --email your.email@gmail.com --dwn-dir PDFs --allow-scihub true
```

## Search PubMed and download results

```python
from pypaperretriever import PubMedSearcher

search_query = "brain lesion case reports"
searcher = PubMedSearcher(search_string=search_query, email="your.email@gmail.com")

results = searcher.search(count=5)
searcher.download_articles(download_directory="PDFs", allow_scihub=True)
```

## Extract images from PDFs

```python
from pypaperretriever import ImageExtractor

extractor = ImageExtractor("PDFs/10.7759_cureus.76081.pdf")
extractor.extract_images()
```

## Track citation networks

```python
from pypaperretriever import PaperTracker

tracker = PaperTracker(
    email="your.email@gmail.com",
    doi="10.1097/RLU.0000000000001894",
    max_upstream_generations=1,
    max_downstream_generations=1,
)
results = tracker.track_paper()
```

## What's next?
Explore the full API reference for more details on available methods and options:
- [API Reference](api/index.md)