"""Tools to build citation networks around scientific papers."""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import pandas as pd

from .reference_retriever import ReferenceRetriever


class PaperTracker:
    """Track references and citations for a given paper.

    Args:
        email (str): Email address used for API requests.
        max_upstream_generations (int, optional): Depth of reference traversal.
        max_downstream_generations (int, optional): Depth of citation traversal.
        doi (str, optional): DOI of the root paper.
        pmid (str, optional): PMID of the root paper.

    Attributes:
        df (pd.DataFrame): Table describing all tracked papers.
        processed_upstream (set[str]): Papers already expanded upstream.
        processed_downstream (set[str]): Papers already expanded downstream.
    """

    def __init__(
        self,
        email: str,
        max_upstream_generations: int = 1,
        max_downstream_generations: int = 1,
        doi: Optional[str] = None,
        pmid: Optional[str] = None,
    ):
        """
        Initialize the PaperTracker with either a DOI or PMID.

        Args:
            email (str): Email for API authentication
            max_upstream_generations (int, optional): Maximum reference depth to track. Defaults to 1.
            max_downstream_generations (int, optional): Maximum citation depth to track. Defaults to 1.
            doi (str, optional): Digital Object Identifier of the paper. Required if pmid not provided.
            pmid (str, optional): PubMed ID of the paper. Required if doi not provided.

        Raises:
            ValueError: If neither DOI nor PMID is provided
        """
        if not doi and not pmid:
            raise ValueError("Either DOI or PMID must be provided.")
        self.email = email
        self.max_upstream_generations = max_upstream_generations
        self.max_downstream_generations = max_downstream_generations
        self.doi = doi
        self.pmid = pmid
        print(f"[PaperTracker] Initializing with DOI: {doi} and PMID: {pmid}")
        self.root_retriever = ReferenceRetriever(email=email, doi=doi, pmid=pmid)
        self.df_columns = [
            'doi', 'pmid', 'title', 'authors', 'year',
            'upstream_generation', 'downstream_generation',
            'children_identifiers', 'parent_identifiers'
        ]
        self.df = pd.DataFrame(columns=self.df_columns)
        self.processed_upstream = set()
        self.processed_downstream = set()

    def go_upstream(self, doi: Optional[str] = None, pmid: Optional[str] = None) -> List[Dict[str, Any]]:
        """Fetch references for a paper.

        Args:
            doi (str, optional): DOI of the paper.
            pmid (str, optional): PMID of the paper.

        Returns:
            list[dict[str, Any]]: Reference metadata.
        """
        print(f"[PaperTracker] Going upstream for DOI: {doi}, PMID: {pmid}")
        retriever = ReferenceRetriever(email=self.email, doi=doi, pmid=pmid)
        return retriever.fetch_references()

    def go_downstream(self, doi: Optional[str] = None, pmid: Optional[str] = None) -> List[Dict[str, Any]]:
        """Fetch papers that cite the given paper.

        Args:
            doi (str, optional): DOI of the paper.
            pmid (str, optional): PMID of the paper.

        Returns:
            list[dict[str, Any]]: Citing paper metadata.
        """
        print(f"[PaperTracker] Going downstream for DOI: {doi}, PMID: {pmid}")
        retriever = ReferenceRetriever(email=self.email, doi=doi, pmid=pmid)
        return retriever.fetch_cited_by()

    def track_paper(self) -> pd.DataFrame:
        """Build the citation network around the root paper.

        Returns:
            pandas.DataFrame: Table containing all tracked papers and relationships.
        """
        print(f"[PaperTracker] Starting tracking process for DOI: {self.doi}, PMID: {self.pmid}")
        self._track_upstream(self.doi, self.pmid, 0, None)
        self._track_downstream(self.doi, self.pmid, 0, None)
        print("[PaperTracker] Tracking process completed.")
        return self.df

    def _track_upstream(
        self,
        doi: Optional[str],
        pmid: Optional[str],
        generation: int,
        parent_id: Optional[str],
    ) -> None:
        """Recursively follow reference chains.

        Args:
            doi (str | None): DOI of the current paper.
            pmid (str | None): PMID of the current paper.
            generation (int): Current depth in the reference tree.
            parent_id (str | None): Identifier of the citing paper.
        """
        paper_id = doi if doi else pmid
        print(f"[PaperTracker] Tracking upstream - Generation: {generation}, Paper ID: {paper_id}, Parent ID: {parent_id}")

        if generation > self.max_upstream_generations:
            print(f"[PaperTracker] Maximum upstream generations reached: {generation} for Paper ID: {paper_id}")
            return

        if paper_id in self.processed_upstream:
            print(f"[PaperTracker] Paper ID {paper_id} already processed upstream. Skipping.")
            return
        else:
            self.processed_upstream.add(paper_id)

        references = self.go_upstream(doi=doi, pmid=pmid)
        print(f"[PaperTracker] Retrieved {len(references)} references for Paper ID: {paper_id}")

        if not isinstance(references, list):
            print(f"[PaperTracker] Unexpected type for references: {type(references)}. Expected list.")
            return

        paper_data = self._get_paper_metadata(doi, pmid)
        if not paper_data:
            print(f"[PaperTracker] No metadata found for {'DOI: ' + doi if doi else 'PMID: ' + pmid}")
            return
        paper_data.update({
            'upstream_generation': generation,
            'downstream_generation': None,
            'children_identifiers': [],
            'parent_identifiers': [parent_id] if parent_id else []
        })

        print(f"[PaperTracker] Adding Paper ID: {paper_id} to tracker DataFrame.")
        self.df = pd.concat([self.df, pd.DataFrame([paper_data])], ignore_index=True)

        for ref in references:
            if not isinstance(ref, dict):
                print(f"[PaperTracker] Unexpected reference format: {ref}. Skipping.")
                continue
            ref_id = ref.get('doi') or ref.get('pmid')
            if not ref_id:
                print(f"[PaperTracker] Reference missing DOI and PMID: {ref}. Skipping.")
                continue
            print(f"[PaperTracker] Processing reference ID: {ref_id}")
            self._track_upstream(ref.get('doi'), ref.get('pmid'), generation + 1, paper_id)
            idx = self.df[self.df['doi'] == paper_id].index
            if not idx.empty:
                current_children = self.df.at[idx[0], 'children_identifiers']
                if not isinstance(current_children, list):
                    current_children = []
                self.df.at[idx[0], 'children_identifiers'] = current_children + [ref_id]
                print(f"[PaperTracker] Updated children_identifiers for Paper ID: {paper_id}")

    def _track_downstream(
        self,
        doi: Optional[str],
        pmid: Optional[str],
        generation: int,
        parent_id: Optional[str],
    ) -> None:
        """Recursively follow citation chains.

        Args:
            doi (str | None): DOI of the current paper.
            pmid (str | None): PMID of the current paper.
            generation (int): Current depth in the citation tree.
            parent_id (str | None): Identifier of the referenced paper.
        """
        paper_id = doi if doi else pmid
        print(f"[PaperTracker] Tracking downstream - Generation: {generation}, Paper ID: {paper_id}, Parent ID: {parent_id}")

        if generation > self.max_downstream_generations:
            print(f"[PaperTracker] Maximum downstream generations reached: {generation} for Paper ID: {paper_id}")
            return

        if paper_id in self.processed_downstream:
            print(f"[PaperTracker] Paper ID {paper_id} already processed downstream. Skipping.")
            return
        else:
            self.processed_downstream.add(paper_id)

        cited_by = self.go_downstream(doi=doi, pmid=pmid)
        print(f"[PaperTracker] Retrieved {len(cited_by)} citing articles for Paper ID: {paper_id}")

        if not isinstance(cited_by, list):
            print(f"[PaperTracker] Unexpected type for cited_by: {type(cited_by)}. Expected list.")
            return

        paper_data = self._get_paper_metadata(doi, pmid)
        if not paper_data:
            print(f"[PaperTracker] No metadata found for {'DOI: ' + doi if doi else 'PMID: ' + pmid}")
            return
        paper_data.update({
            'upstream_generation': None,
            'downstream_generation': generation,
            'children_identifiers': [],
            'parent_identifiers': [parent_id] if parent_id else []
        })

        print(f"[PaperTracker] Adding Paper ID: {paper_id} to tracker DataFrame.")
        self.df = pd.concat([self.df, pd.DataFrame([paper_data])], ignore_index=True)

        for cite in cited_by:
            if not isinstance(cite, dict):
                print(f"[PaperTracker] Unexpected citing article format: {cite}. Skipping.")
                continue
            cite_id = cite.get('doi') or cite.get('pmid')
            if not cite_id:
                print(f"[PaperTracker] Citing article missing DOI and PMID: {cite}. Skipping.")
                continue
            print(f"[PaperTracker] Processing citing article ID: {cite_id}")
            self._track_downstream(cite.get('doi'), cite.get('pmid'), generation + 1, paper_id)
            idx = self.df[self.df['doi'] == paper_id].index
            if not idx.empty:
                current_children = self.df.at[idx[0], 'children_identifiers']
                if not isinstance(current_children, list):
                    current_children = []
                self.df.at[idx[0], 'children_identifiers'] = current_children + [cite_id]
                print(f"[PaperTracker] Updated children_identifiers for Paper ID: {paper_id}")

    def _get_paper_metadata(self, doi: Optional[str], pmid: Optional[str]) -> Dict[str, Any]:
        """Retrieve metadata for a paper.

        Args:
            doi (str | None): DOI of the paper.
            pmid (str | None): PMID of the paper.

        Returns:
            dict[str, Any]: Metadata including identifiers, title, authors and year.
        """
        print(f"[PaperTracker] Fetching metadata for DOI: {doi}, PMID: {pmid}")
        retriever = ReferenceRetriever(email=self.email, doi=doi, pmid=pmid)
        metadata = retriever.get_paper_metadata()
        if not metadata:
            print(f"[PaperTracker] No metadata retrieved for DOI: {doi}, PMID: {pmid}")
            return {}
        print(f"[PaperTracker] Metadata retrieved for DOI: {doi}, PMID: {pmid}")
        return {
            'doi': metadata.get('doi', ''),
            'pmid': metadata.get('pmid', ''),
            'title': metadata.get('title', ''),
            'authors': metadata.get('authors', ''),
            'year': metadata.get('year', ''),
            'upstream_generation': None,
            'downstream_generation': None,
            'children_identifiers': [],
            'parent_identifiers': []
        }
