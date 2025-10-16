from Bio import Entrez, Medline
import pandas as pd
import os

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
Entrez.email = "your.email@example.com"  # Replace with your actual email

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
        doi_list = r.get("AID", [])
        doi = ""
        if doi_list:
            for d in doi_list:
                if "[doi]" in d:
                    doi = d.replace("[doi]", "").strip()
                    break
        authors = r.get("AU", [])
        # Use et al. truncation (first 3 authors only)
        if len(authors) > 3:
            formatted_authors = ", ".join(authors[:3]) + " et al."
        else:
            formatted_authors = ", ".join(authors)

        year = r.get("DP", "").split(" ")[0]
        key = (authors[0].split()[-1][:3] + year) if authors else "Ref" + year
        records_data.append({
            "PMID": str(r.get("PMID", "")),
            "Title": title,
            "Authors": formatted_authors,
            "Year": year,
            "DOI": doi,
            "Key": key
        })
    handle.close()

# ----------------------------
# 5. Write BibTeX
# ----------------------------
with open("references.bib", "w", encoding="utf-8") as f:
    for rec in records_data:
        f.write("@article{{{}}},\n".format(rec["Key"]))
        f.write("  author = {{{}}},\n".format(rec["Authors"]))
        f.write("  title = {{{}}},\n".format(rec["Title"]))
        f.write("  journal = {PubMed},\n")
        f.write("  year = {{{}}},\n".format(rec["Year"]))
        if rec["DOI"]:
            f.write("  doi = {{{}}},\n".format(rec["DOI"]))
        f.write("}\n\n")

# ----------------------------
# 6. Save CSV/JSON for Excel with Include column preserved
# ----------------------------
records_df_new = pd.DataFrame(records_data)
csv_path = "records.csv"
json_path = "records.json"

if os.path.exists(csv_path):
    # Read previous version to preserve Include column
    records_df_old = pd.read_csv(csv_path, sep=";", encoding="utf-8", dtype={"PMID": str})
    records_df_new["PMID"] = records_df_new["PMID"].astype(str)

    if "Include" in records_df_old.columns:
        records_df = pd.merge(
            records_df_new,
            records_df_old[["PMID", "Include"]],
            on="PMID",
            how="left"
        )
        records_df["Include"] = records_df["Include"].fillna(0).astype(int)
    else:
        records_df = records_df_new.copy()
        records_df.insert(0, "Include", 0)
else:
    records_df = records_df_new.copy()
    records_df.insert(0, "Include", 0)

records_df.to_csv(csv_path, index=False, sep=";", encoding="utf-8")
records_df.to_json(json_path, orient="records", force_ascii=False)

# ----------------------------
# 7. PRISMA counts
# ----------------------------
included = records_df["Include"].sum()
excluded_fulltext = len(records_df) - included

# ----------------------------
# 8. Generate LaTeX report
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
\\usepackage{{ragged2e}}
\\usetikzlibrary{{positioning}}

\\geometry{{top=0.8in, bottom=0.8in, left=0.7in, right=0.7in}}
\\addbibresource{{references.bib}}

\\title{{Systematic Evidence on {search_topic} (2015--2025)}}
\\author{{Automated Systematic Search (Entrez API)}}
\\date{{\\today}}

\\begin{{document}}
\\maketitle

\\section*{{Objective}}
Summarize international guidelines, consensus statements, and related literature (2015--2025) regarding \\textbf{{{search_topic}}}.

\\section*{{Methodological Approach}}
\\begin{{itemize}}
    \\item \\textbf{{Database:}} PubMed (primary)
    \\item \\textbf{{Time frame:}} {date_from} -- {date_to}
    \\item \\textbf{{Query:}} \\texttt{{{query}}}
    \\item \\textbf{{Screening:}} See PRISMA diagram
\\end{{itemize}}

\\section*{{PRISMA Flow Diagram}}
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
\\renewcommand{{\\arraystretch}}{{1.3}}
\\footnotesize
\\setlength{{\\extrarowheight}}{{2pt}}
\\begin{{longtable}}{{|p{{1cm}}|p{{3.5cm}}|p{{1cm}}|p{{8cm}}|p{{1.5cm}}|}}
\\hline
\\textbf{{Inc}} & \\textbf{{Authors}} & \\textbf{{Year}} & \\textbf{{Title}} & \\textbf{{Ref}} \\\\
\\hline
\\endfirsthead
\\hline
\\textbf{{Inc}} & \\textbf{{Authors}} & \\textbf{{Year}} & \\textbf{{Title}} & \\textbf{{Ref}} \\\\
\\hline
\\endhead
\\hline
\\endfoot
\\hline
\\endlastfoot
"""

for _, rec in records_df.iterrows():
    authors = rec["Authors"].replace("&", "\\&")
    doi_link = f"\\href{{https://doi.org/{rec['DOI']}}}{{{rec['Key']}}}" if rec["DOI"] else rec["Key"]
    latex_content += f"{rec['Include']} & \\RaggedRight {authors} & {rec['Year']} & \\RaggedRight {rec['Title']} & {doi_link} \\\\\n\\hline\n"

latex_content += """
\\end{longtable}

\\section*{Transparency and Reproducibility}
\\begin{itemize}
    \\item Author lists truncated to first 3 names followed by “et al.” for clarity.
    \\item UTF-8 semicolon-delimited CSV compatible with Excel (Portuguese locale).
    \\item Inclusion decisions persist across runs.
    \\item PRISMA diagram updates automatically.
\\end{itemize}

\\printbibliography
\\end{document}
"""

with open("systematic_review.tex", "w", encoding="utf-8") as f:
    f.write(latex_content)

print("✅ LaTeX file 'systematic_review.tex' and BibTeX 'references.bib' generated successfully!")
print("✅ CSV and JSON updated (records.csv, records.json).")
