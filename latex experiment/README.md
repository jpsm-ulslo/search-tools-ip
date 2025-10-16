Systematic Literature Search & LaTeX Report Generator

Description:
This repository contains the Python script wally.py, which performs automated literature searches on PubMed using the Entrez API, extracts metadata, and generates a reproducible LaTeX report with a PRISMA-style flow diagram and an evidence table.

Features:

- Search PubMed based on a customizable query, topic, and date range.
- Extract metadata: title, authors, year, DOI, abstract.
- Generate a .bib file for citations.
- Produce a LaTeX report with:
  - PRISMA flow diagram for screening and inclusion/exclusion
  - Synthesized evidence table
  - Clickable DOI links in the table
- Export results to CSV or JSON for Excel or other analysis.

Installation:

1. Clone the repository:
   https://github.com/jpsm-ulslo/search-tools-ip/
2. Install Python dependencies:
   pip install biopython pandas

Usage:

1. Edit the wally.py file to set your search parameters:

   search_topic = "indirect calorimetry"
   query = '"indirect calorimetry"[tiab] AND ("guideline"[tiab] OR "consensus"[tiab]) AND ("energy expenditure"[tiab] OR "resting energy expenditure"[tiab])'
   date_from = "2015/01/01"
   date_to = "2025/12/31"

2. Run the script:
   python wally.py

Outputs:

- systematic_review.tex → LaTeX report
- references.bib → BibTeX file
- records.csv → CSV file compatible with Excel
- records.json → JSON file

Excel / CSV Considerations:

- records.csv is UTF-8 encoded, compatible with Excel in any language (including Portuguese from Portugal).
- Columns include:
  - Include → Mark 1 / Yes to include, 0 / No to exclude
  - PMID → PubMed ID
  - Title → Article title
  - Authors → List of authors
  - Year → Year of publication
  - DOI → Digital Object Identifier
  - Key → Unique BibTeX key

> Note: Editing the Include column allows controlling which papers appear in the evidence table and PRISMA flow diagram. The script reads this column to update the LaTeX report.

Customization & Reproducibility:

- Change search_topic, query, date_from, and date_to in wally.py to perform searches on different topics.
- Modify LaTeX templates inside the script to adjust report styling or table formatting.
- The CSV file allows easy marking of papers for inclusion/exclusion in subsequent runs, which updates the PRISMA flow diagram.

Example:

1. Run wally.py.
2. Open records.csv in Excel.
3. Mark papers to include in the Include column.
4. Rerun the script to regenerate the LaTeX report reflecting inclusion decisions.

License:

MIT License
