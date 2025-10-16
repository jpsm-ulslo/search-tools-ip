Name: Systematic Literature Search & LaTeX Report Generator
Description: |
  This manifest describes the Python script `wally.py` which performs automated literature searches
  on PubMed using the Entrez API, extracts metadata, and generates a reproducible LaTeX report
  with a PRISMA-style flow diagram and an evidence table.

Features: |
  - Search PubMed based on a customizable query, topic, and date range.
  - Extract metadata: title, authors, year, DOI, abstract.
  - Generate a .bib file for citations.
  - Produce a LaTeX report with:
      * PRISMA flow diagram for screening and inclusion/exclusion
      * Synthesized evidence table
      * Clickable DOI links in the table
  - Export results to CSV or JSON for Excel or other analysis.

Installation: |
  1. Clone this repository:
     git clone https://github.com/your-username/repository-name.git
     cd repository-name
  2. Install Python dependencies:
     pip install biopython pandas

Usage: |
  Edit the `wally.py` file to set your search parameters:
    search_topic = "indirect calorimetry"
    query = '"indirect calorimetry"[tiab] AND ("guideline"[tiab] OR "consensus"[tiab]) AND ("energy expenditure"[tiab] OR "resting energy expenditure"[tiab])'
    date_from = "2015/01/01"
    date_to = "2025/12/31"

  Run the script:
    python wally.py

  Outputs:
    - systematic_review.tex → LaTeX report
    - references.bib → BibTeX file
    - records.csv → CSV file compatible with Excel
    - records.json → JSON file

Excel / CSV Considerations: |
  - `records.csv` is UTF-8 encoded, compatible with Excel in any language (including Portuguese from Portugal).
  - Columns include:
      Include, PMID, Title, Authors, Year, DOI, Key

Customization & Reproducibility: |
  - Change `search_topic`, `query`, `date_from`, and `date_to` in `wally.py` to perform searches for different topics.
  - Modify LaTeX templates inside the script to adjust report styling or table formatting.
  - The CSV file allows easy marking of papers for inclusion/exclusion in subsequent runs, which can update the PRISMA flow diagram.

Example: |
  After running `wally.py`, open `records.csv` in Excel, mark papers to include, and rerun the script to update the LaTeX report.

License: MIT License
