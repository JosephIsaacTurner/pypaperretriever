# API Overview

The public API is organized by module (mirrors the package layout):

- `paper_retriever` — single-paper retrieval with OA-first logic and optional Sci-Hub fallback
- `pubmed_searcher` — programmatic PubMed search → DataFrame pipeline
- `image_extractor` — extract figures from PDFs
- `reference_retriever` — fetch references and cited-by
- `paper_tracker` — upstream/downstream citation network
- `utils` — small helpers for power users
