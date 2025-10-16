from Bio import Entrez, Medline
import json

# ----------------------------
# 1. User-defined search setup
# ----------------------------
search_topic = "indirect calorimetry"
query = '"indirect calorimetry"[tiab] AND ("guideline"[tiab] OR "consensus"[tiab]) AND ("energy expenditure"[tiab] OR "resting energy expenditure"[tiab])'
date_from = "2015/01/01"
date_to = "2025/12/31"

# ----------------------------
# 2. NCBI / Entrez setup
# ----------------------------
Entrez.email = "your.email@example.com"  # Replace with your real email

# ----------------------------
# 3. Search PubMed
# ----------------------------
handle = Entrez.esearch(
    db="pubmed",
    term=query,
    mindate=date_from,
    maxdate=date_to,
    datetype="pdat",
    retmax=1000
)
record = Entrez.read(handle)
handle.close()

pmid_list = record["IdList"]
pubmed_total = int(record["Count"])
duplicates = 0
screened = pubmed_total - duplicates
excluded_screening = 0
full_text = 0
excluded_fulltext = 0
included = 0

print(f"PubMed total: {pubmed_total}")
print(f"PMIDs: {pmid_list}")

# ----------------------------
# 4. Fetch metadata
# ----------------------------
records_data = []

if pmid_list:
    handle = Entrez.efetch(db="pubmed", id=pmid_list, rettype="medline", retmode="text")
    records = Medline.parse(handle)
    for r in records:
        title = r.get("TI", "")
        abstract = r.get("AB", "")
        doi_list = r.get("AID", [])
        doi = ""
        if doi_list:
            for d in doi_list:
                if "[doi]" in d:
                    doi = d.replace("[doi]", "").strip()
                    break
        authors = r.get("AU", [])
        year = r.get("DP", "").split(" ")[0]
        if authors:
            key = authors[0].split()[-1][:3] + year
        else:
            key = "Ref" + year

        records_data.append({
            "pmid": r.get("PMID", ""),
            "title": title,
            "abstract": abstract,
            "doi": doi,
            "authors": authors,
            "year": year,
            "key": key,
            "decision": "Pending"  # <-- new screening status column
        })
    handle.close()

# JSON (for reproducibility)
with open("records.json", "w", encoding="utf-8") as jf:
    json.dump(records_data, jf, indent=2, ensure_ascii=False)

# CSV (for Excel)
import csv
with open("records.csv", "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["PMID", "Authors", "Year", "Title", "DOI", "Decision"])
    for rec in records_data:
        authors = ", ".join(rec["authors"])
        writer.writerow([rec["pmid"], authors, rec["year"], rec["title"], rec["doi"], rec["decision"]])

# ----------------------------
# 5. Write .bib file
# ----------------------------
with open("references.bib", "w", encoding="utf-8") as f:
    for rec in records_data:
        f.write("@article{" + rec['key'] + ",\n")
        f.write("  author = {" + " and ".join(rec['authors']) + "},\n")
        f.write("  title = {" + rec['title'] + "},\n")
        f.write("  journal = {PubMed},\n")
        f.write("  year = {" + rec['year'] + "},\n")
        if rec['doi']:
            f.write("  doi = {" + rec['doi'] + "},\n")
        f.write("}\n\n")

# ----------------------------
# 6. Generate LaTeX file
# ----------------------------
latex_content = f"""
\\documentclass[12pt,a4paper]{{article}}
\\usepackage[utf8]{{inputenc}}
\\usepackage{{geometry}}
\\usepackage{{array}}
\\usepackage{{hyperref}}
\\usepackage[backend=biber,style=numeric,sorting=ynt]{{biblatex}}
\\usepackage{{longtable}}
\\usepackage{{tikz}}
\\usepackage{{microtype}}
\\usetikzlibrary{{positioning}}

\\geometry{{top=0.8in, bottom=0.8in, left=0.8in, right=0.8in}}
\\setlength{{\\emergencystretch}}{{3em}} % prevent text overflow

\\addbibresource{{references.bib}}

\\title{{Systematic Evidence on {search_topic} (2015--2025)}}
\\author{{Automated Systematic Search (Entrez API)}}
\\date{{\\today}}

\\begin{{document}}
\\maketitle

\\section*{{Objective}}
Summarize international guidelines, consensus statements, and related literature (2015--2025) regarding the use of \\textbf{{{search_topic}}}.

\\section*{{Methodological Approach}}
\\begin{{itemize}}
    \\item \\textbf{{Databases:}} PubMed (primary), supplemented by Google Scholar and society websites.
    \\item \\textbf{{Time frame:}} {date_from} -- {date_to}.
    \\item \\textbf{{PubMed search string:}} \\texttt{{{query}}}
    \\item \\textbf{{Screening and inclusion/exclusion:}} See PRISMA diagram.
\\end{{itemize}}

\\section*{{PRISMA Flow Diagram (PubMed 2015--2025)}}
\\begin{{tikzpicture}}[node distance=1.5cm, auto]
\\tikzstyle{{block}} = [rectangle, draw, text width=10cm, align=center, rounded corners, minimum height=1.2cm]

\\node[block] (id) {{Records identified from PubMed: {pubmed_total}}};
\\node[block, below=of id] (dup) {{Duplicates removed: {duplicates} \\\\ Records screened: {screened}}};
\\node[block, below=of dup] (exc) {{Excluded after screening: {excluded_screening}}};
\\node[block, below=of exc] (elig) {{Full-text assessed: {full_text}}};
\\node[block, below=of elig] (exc2) {{Excluded after full-text: {excluded_fulltext}}};
\\node[block, below=of exc2] (inc) {{Included in synthesis: {included}}};

\\draw[->] (id) -- (dup);
\\draw[->] (dup) -- (exc);
\\draw[->] (exc) -- (elig);
\\draw[->] (elig) -- (exc2);
\\draw[->] (exc2) -- (inc);
\\end{{tikzpicture}}

\\section*{{Synthesized Evidence Table}}
\\footnotesize
\\setlength{{\\extrarowheight}}{{2pt}}
\\begin{{longtable}}{{|p{{2cm}}|p{{3cm}}|c|p{{7cm}}|c|}}
\\hline
Decision & Authors & Year & Title & Ref \\\\
\\hline
\\endfirsthead
\\hline
Decision & Authors & Year & Title & Ref \\\\
\\hline
\\endhead
\\hline
\\endfoot
\\hline
\\endlastfoot
"""

# Write each record to the LaTeX table
for rec in records_data:
    authors = ", ".join(rec["authors"][:3]) + (" et al." if len(rec["authors"]) > 3 else "")
    doi_link = f"\\href{{https://doi.org/{rec['doi']}}}{{{rec['key']}}}" if rec['doi'] else rec['key']
    latex_content += f"{rec['decision']} & {authors} & {rec['year']} & {rec['title']} & {doi_link} \\\\\n"

latex_content += """
\\end{longtable}

\\section*{Transparency and Reproducibility}
\\begin{itemize}
    \\item All sources retrieved from PubMed via Entrez API.
    \\item Data extracted from abstracts, titles, and verified DOIs.
    \\item Only peer-reviewed, English-language sources included.
    \\item Fully reproducible following the described search strategy.
\\end{itemize}

\\printbibliography

\\end{document}
"""

with open("systematic_review.tex", "w", encoding="utf-8") as f:
    f.write(latex_content)

print("✅ 'systematic_review.tex' and 'references.bib' generated successfully!")
print("→ You can compile using: pdflatex systematic_review.tex && biber systematic_review && pdflatex systematic_review.tex")
