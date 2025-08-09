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
    download_directory="PDFs"
)
retriever.download()
```

To download via PubMed ID:

```python
from pypaperretriever import PaperRetriever

retriever = PaperRetriever(
    email="your.email@gmail.com",
    pmid="33813262",
    download_directory="PDFs"
)
retriever.download()
```

## Search PubMed and download results

```python
from pypaperretriever import PubMedSearcher

search_query = "brain lesion case reports"
searcher = PubMedSearcher(search_string=search_query, email="your.email@gmail.com")

results = searcher.search(count=5)
searcher.download_articles(download_directory="PDFs")
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
