Systematic Literature Search & LaTeX Report Generator
Description

This repository contains the Python script wally.py which automates literature searches on PubMed using the Entrez API. It extracts metadata, generates a reproducible LaTeX report with a PRISMA-style flow diagram, and produces a synthesized evidence table with clickable DOI links.

Features

Search PubMed based on a customizable query, topic, and date range.

Extract metadata: Title, Authors, Year, DOI, Abstract.

Generate a .bib file for LaTeX citations.

Produce a LaTeX report with:

PRISMA flow diagram for screening and inclusion/exclusion.

Synthesized evidence table.

Clickable DOI links in the table.

Export results to CSV and JSON for Excel or further analysis.

Include/Exclude workflow via CSV marking.

Installation

Clone the repository:
git clone https://github.com/your-username/repository-name.git

cd repository-name

Install Python dependencies:
pip install biopython pandas

Usage

Edit wally.py to set your search parameters:
search_topic = "indirect calorimetry"
query = '"indirect calorimetry"[tiab] AND ("guideline"[tiab] OR "consensus"[tiab]) AND ("energy expenditure"[tiab] OR "resting energy expenditure"[tiab])'
date_from = "2015/01/01"
date_to = "2025/12/31"

Run the script:
python wally.py

Outputs:

systematic_review.tex → LaTeX report

references.bib → BibTeX file

records.csv → CSV file compatible with Excel

records.json → JSON file

Excel / CSV Considerations

records.csv is UTF-8 encoded, compatible with Excel in any language (including Portuguese from Portugal).

Columns include: Include, PMID, Title, Authors, Year, DOI, Key.

The Include column can be updated to 0 (exclude) or 1 (include). Re-running wally.py will update the LaTeX report and PRISMA flow diagram automatically.

Customization & Reproducibility

Change search_topic, query, date_from, and date_to in wally.py for different topics.

Modify LaTeX templates inside the script to adjust report styling or table formatting.

CSV marking allows flexible decision-making for which papers to include in the final synthesis.

JSON export provides a reproducible machine-readable record of all results.

Example Workflow

Run wally.py to fetch PubMed results.

Open records.csv in Excel.

Mark papers to include (1) or exclude (0) in the Include column.

Save the CSV and re-run wally.py to update:

The PRISMA flow diagram in the LaTeX report.

The synthesized evidence table.

Clickable DOI links for included papers.

License

MIT License
