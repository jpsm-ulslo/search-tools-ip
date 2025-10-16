from Bio import Entrez, Medline
import csv
import json
import os
from datetime import datetime

# ===============================================================
# 1. USER-DEFINED SEARCH SETUP (EASY TO CUSTOMIZE)
# ===============================================================
search_topic = "indirect calorimetry"
query = '"indirect calorimetry"[tiab] AND ("guideline"[tiab] OR "consensus"[tiab]) AND ("energy expenditure"[tiab] OR "resting energy expenditure"[tiab])'
date_from = "2015/01/01"
date_to = "2025/12/31"
Entrez.email = "your.email@example.com"  # REQUIRED for NCBI access

# ===============================================================
# 2. HELPER: CHECK IF CSV EXISTS (avoid overwriting manual edits)
# ===============================================================
csv_filename = "records.csv"
json_filename = "records.json"
bib_filename = "references.bib"
tex_filename = "systematic_review.tex"

# ===============================================================
# 3. IF CSV EXISTS → READ DECISIONS AND SKIP NEW FETCH
# ===============================================================
records_data = []

if os.path.exists(csv_filename):
    print(f"📂 Existing '{csv_filename}' detected — reading manual decisions...")
    with open(csv_filename, encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter=';')
        records_data = [row for row in reader]

    # Count PRISMA elements
    total_records = len(records_data)
    included = sum(1 for r in records_data if r.get("Decision", "").lower() == "include")
    excluded = sum(1 for r in records_data if r.get("Decision", "").lower() == "exclude")
    screened = total_records
    duplicates = 0
    full_text = included + excluded
    excluded_fulltext = excluded
    print(f"✅ Decisions loaded: {included} included / {excluded} excluded")

else:
    # ===========================================================
    # 4. FETCH FROM PUBMED (FIRST RUN ONLY)
    # ===========================================================
    print("🔎 Fetching records from PubMed...")
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
    total_records = int(record["Count"])
    duplicates = 0
    screened = total_records - duplicates
    excluded = 0
    full_text = 0
    excluded_fulltext = 0
    included = 0

    print(f"📄 PubMed total: {total_records}")

    if pmid_list:
        handle = Entrez.efetch(db="pubmed", id=pmid_list, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        for r in records:
            title = r.get("TI", "")
            doi = ""
            for d in r.get("AID", []):
                if "[doi]" in d:
                    doi = d.replace("[doi]", "").strip()
                    break
            authors = r.get("AU", [])
            year = r.get("DP", "").split(" ")[0]
            key = (authors[0].split()[-1][:3] + year) if authors else ("Ref" + year)

            records_data.append({
                "PMID": r.get("PMID", ""),
                "Authors": "; ".join(authors),
                "Year": year,
                "Title": title,
                "DOI": doi,
                "Key": key,
                "Decision": ""  # blank for later inclusion/exclusion marking
            })
        handle.close()

        # -------------------------------------------------------
        # Save results for manual review
        # -------------------------------------------------------
        with open(csv_filename, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(records_data[0].keys()), delimiter=';')
            writer.writeheader()
            writer.writerows(records_data)
        with open(json_filename, "w", encoding="utf-8") as f:
            json.dump(records_data, f, ensure_ascii=False, indent=2)

        print(f"💾 Saved '{csv_filename}' and '{json_filename}' — ready for Excel review.")

# ===============================================================
# 5. WRITE .BIB FILE (for LaTeX)
# ===============================================================
with open(bib_filename, "w", encoding="utf-8") as f:
    for rec in records_data:
        f.write("@article{" + rec['Key'] + ",\n")
        f.write("  author = {" + rec['Authors'].replace(";", " and") + "},\n")
        f.write("  title = {" + rec['Title'] + "},\n")
        f.write("  journal = {PubMed},\n")
        f.write("  year = {" + rec['Year'] + "},\n")
        if rec['DOI']:
            f.write("  doi = {" + rec['DOI'] + "},\n")
        f.write("}\n\n")

# ===============================================================
# 6. GENERATE LATEX FILE (with PRISMA + clickable DOIs)
# ===============================================================
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

\\title{{Systematic Evidence on {search_topic} ({date_from}--{date_to})}}
\\author{{Automated Systematic Search (Entrez API)}}
\\date{{\\today}}

\\begin{{document}}
\\maketitle

\\section*{{Objective}}
Summarize international guidelines, consensus statements, and related literature ({date_from}--{date_to}) on \\textbf{{{search_topic}}}.

\\section*{{Methodological Approach}}
\\begin{{itemize}}
    \\item \\textbf{{Database:}} PubMed (primary)
    \\item \\textbf{{Time frame:}} {date_from} -- {date_to}
    \\item \\textbf{{Search string:}} \\texttt{{{query}}}
    \\item \\textbf{{Screening process:}} See PRISMA diagram below
\\end{{itemize}}

\\section*{{PRISMA Flow Diagram}}
\\begin{{tikzpicture}}[node distance=1.5cm, auto]
\\tikzstyle{{block}} = [rectangle, draw, text width=10cm, align=center, rounded corners, minimum height=1.2cm]

\\node[block] (id) {{Records identified from PubMed: {total_records}}};
\\node[block, below=of id] (dup) {{Duplicates removed: {duplicates} \\\\ Records screened: {screened}}};
\\node[block, below=of dup] (exc) {{Excluded after screening: {excluded}}};
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
\\begin{{longtable}}{{|p{{1.8cm}}|p{{3cm}}|c|p{{7cm}}|c|}}
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

for rec in records_data:
    authors = rec["Authors"].split(";")[0:3]
    author_text = ", ".join(a.strip() for a in authors) + (" et al." if len(authors) >= 3 else "")
    doi_link = f"\\href{{https://doi.org/{rec['DOI']}}}{{{rec['Key']}}}" if rec['DOI'] else rec['Key']
    decision = rec.get("Decision", "")
    latex_content += f"{decision} & {author_text} & {rec['Year']} & {rec['Title']} & {doi_link} \\\\\n"

latex_content += """
\\end{longtable}

\\section*{Transparency and Reproducibility}
\\begin{itemize}
    \\item All data retrieved via Entrez API (NCBI)
    \\item Manual screening performed using exported CSV
    \\item BibTeX and LaTeX generated automatically
    \\item PRISMA diagram updates dynamically from 'Decision' column
\\end{itemize}

\\printbibliography

\\end{document}
"""

with open(tex_filename, "w", encoding="utf-8") as f:
    f.write(latex_content)

print(f"📘 Generated '{tex_filename}', '{bib_filename}', and PRISMA summary.")
