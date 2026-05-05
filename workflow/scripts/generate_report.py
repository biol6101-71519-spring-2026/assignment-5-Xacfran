import os
import csv
import datetime
from pathlib import Path

# Output
out_path   = Path(snakemake.output[0])
report_dir = out_path.parent
report_dir.mkdir(parents=True, exist_ok=True)


# Utilities

def rel_pdf(pdf_path: str) -> str:
    """Relative URL from the report HTML file to a PDF result file."""
    return os.path.relpath(pdf_path, report_dir).replace("\\", "/")


def read_text(path: str) -> str:
    try:
        return Path(path).read_text(encoding="utf-8").strip()
    except Exception:
        return "(file not available)"


def safe_float(val: str) -> str:
    """Format a string as a float with 2 dp if numeric, else return as-is."""
    try:
        return f"{float(val):.2f}"
    except (ValueError, TypeError):
        return str(val) if val else "—"


#  AMAS alignment statistics

def parse_amas_stats(path: str) -> tuple[list, list]:
    """
    AMAS.py summary produces a TSV with a header row.
    Returns (headers, rows) where each row is a list of strings.
    Multiple gene files are concatenated; we keep the header only once.
    """
    headers: list[str] = []
    rows: list[list[str]] = []
    seen_header = False
    try:
        with open(path, newline="", encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                # Normalise first field in case two files were cat'd without a
                # trailing newline: e.g. "0Alignment_name" → strip leading digits
                first = parts[0].lstrip("0123456789")
                if first.startswith("Alignment_name"):
                    parts[0] = first          # clean up the field
                    if not seen_header:
                        headers = parts
                        seen_header = True
                    continue
                rows.append(parts)
    except FileNotFoundError:
        pass
    return headers, rows


#  model summary

def parse_model_summary(path: str) -> tuple[list, list]:
    """
    Parses the tabular section of best_models_summary.txt.
    """
    headers: list[str] = []
    rows: list[list[str]] = []
    try:
        with open(path, encoding="utf-8") as fh:
            in_table = False
            for line in fh:
                line = line.rstrip()
                if line.startswith("Gene") and "Best Model" in line:
                    headers = line.split()
                    in_table = True
                    continue
                if in_table and set(line.strip()) <= {"-", " ", ""}:
                    continue
                if in_table and line.startswith(("Total", "Unique", "=")):
                    break
                if in_table and line:
                    rows.append(line.split())
    except FileNotFoundError:
        pass
    return headers, rows


#  Parse RF distance table

def parse_rf_table(path: str) -> tuple[list, list]:
    """
    Parses robinson_foulds_distances.txt into (headers, rows).
    Expects lines like:  TreeA   TreeB   RF
    after the separator.
    """
    headers: list[str] = []
    rows: list[list[str]] = []
    try:
        with open(path, encoding="utf-8") as fh:
            in_table = False
            for line in fh:
                line = line.rstrip()
                if line.startswith("Tree A") or line.startswith("Tree A"):
                    headers = [x.strip() for x in line.split() if x.strip()]
                    in_table = True
                    continue
                if in_table and set(line.strip()) <= {"-", " ", ""}:
                    continue
                if in_table and line:
                    parts = line.split()
                    if len(parts) >= 3:
                        rows.append(parts)
    except FileNotFoundError:
        pass
    return headers, rows


#  HTML builders

def html_table(headers: list, rows: list, css_id: str = "") -> str:
    if not rows:
        return '<p class="no-data">No data available.</p>'
    id_attr = f' id="{css_id}"' if css_id else ""
    cells = "".join(f"<th>{h}</th>" for h in headers)
    head  = f"<thead><tr>{cells}</tr></thead>"
    body_rows = []
    for row in rows:
        tds = "".join(f"<td>{v}</td>" for v in row)
        body_rows.append(f"<tr>{tds}</tr>")
    body = "<tbody>" + "".join(body_rows) + "</tbody>"
    return f'<div class="table-wrapper"><table{id_attr} class="data-table">{head}{body}</table></div>'


def embed_pdf(pdf_path: str, title: str, height: int = 520) -> str:
    rel = rel_pdf(pdf_path)
    return f"""
    <div class="pdf-block">
      <p class="pdf-title">{title}</p>
      <object data="{rel}" type="application/pdf"
              width="100%" height="{height}px">
        <p>Your browser cannot display PDFs inline.
           <a href="{rel}" target="_blank">Open {title}</a></p>
      </object>
    </div>"""


#  Gather inputs

stats_headers, stats_rows = parse_amas_stats(snakemake.input.summary)
model_headers, model_rows = parse_model_summary(snakemake.input.models)
rf_headers,    rf_rows    = parse_rf_table(snakemake.input.rf)

topo_text        = read_text(snakemake.input.topo)
concordance_text = read_text(snakemake.input.concordance)

# Per-gene PDF list
gene_pdfs  = list(snakemake.input.gene_pdfs) \
             if not isinstance(snakemake.input.gene_pdfs, str) \
             else [snakemake.input.gene_pdfs]
gene_names = [Path(p).stem.replace("_tree", "") for p in gene_pdfs]

gene_tree_embeds = "\n".join(
    embed_pdf(p, f"{g} — gene tree (ML consensus + bootstrap)")
    for g, p in zip(gene_names, gene_pdfs)
)

now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

#  Assemble HTML

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Phylogenetic Analysis Report</title>
  <style>
    /*  Reset & base  */
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: 'Segoe UI', Arial, sans-serif;
      font-size: 15px; line-height: 1.6;
      background: #f0f4f8; color: #2d3748;
    }}
    a {{ color: #2b6cb0; }}

    /*  Layout  */
    .page-header {{
      background: linear-gradient(135deg, #1a365d 0%, #2b6cb0 100%);
      color: #fff; padding: 36px 48px 28px;
    }}
    .page-header h1 {{ font-size: 2em; font-weight: 700; letter-spacing: -0.5px; }}
    .page-header .meta {{ margin-top: 8px; opacity: 0.82; font-size: 0.92em; }}
    .page-header .tags {{ margin-top: 12px; }}
    .tag {{
      display: inline-block; background: rgba(255,255,255,0.18);
      border-radius: 12px; padding: 2px 10px; font-size: 0.82em;
      margin-right: 6px; margin-top: 4px;
    }}

    .container {{ max-width: 1100px; margin: 32px auto; padding: 0 20px 60px; }}

    /*  TOC  */
    .toc {{
      background: #fff; border-radius: 8px; padding: 20px 28px;
      margin-bottom: 28px; box-shadow: 0 1px 4px rgba(0,0,0,0.07);
    }}
    .toc h2 {{ font-size: 1em; text-transform: uppercase; letter-spacing: 0.08em;
               color: #718096; margin-bottom: 10px; }}
    .toc ol {{ padding-left: 18px; }}
    .toc li {{ margin: 4px 0; }}
    .toc a {{ text-decoration: none; color: #2b6cb0; }}
    .toc a:hover {{ text-decoration: underline; }}

    /*  Sections  */
    section {{
      background: #fff; border-radius: 8px; padding: 28px 32px;
      margin-bottom: 24px; box-shadow: 0 1px 4px rgba(0,0,0,0.07);
    }}
    section h2 {{
      color: #1a365d; font-size: 1.3em; font-weight: 700;
      border-bottom: 2px solid #bee3f8; padding-bottom: 8px;
      margin-bottom: 18px;
    }}
    section h3 {{
      color: #2b6cb0; font-size: 1.05em; font-weight: 600;
      margin: 24px 0 10px;
    }}
    p {{ margin-bottom: 10px; }}

    /*  Tables  */
    .table-wrapper {{
      overflow-x: auto;
      -webkit-overflow-scrolling: touch;
      margin: 12px 0 20px;
    }}
    .data-table {{
      width: 100%; min-width: 600px; border-collapse: collapse;
      font-size: 0.86em;
    }}
    .data-table th {{
      background: #2b6cb0; color: #fff;
      padding: 9px 14px; text-align: left;
      font-weight: 600; font-size: 0.9em;
    }}
    .data-table td {{
      padding: 8px 14px; border-bottom: 1px solid #e2e8f0;
      vertical-align: top;
    }}
    .data-table tr:nth-child(even) td {{ background: #f7fafc; }}
    .data-table tr:hover td {{ background: #ebf4ff; }}
    .no-data {{ color: #a0aec0; font-style: italic; margin: 8px 0; }}

    /*  Pre / code blocks  */
    pre {{
      background: #f7fafc; border: 1px solid #e2e8f0; border-radius: 6px;
      padding: 16px 20px; font-size: 0.82em; font-family: 'Consolas',
      'Courier New', monospace; overflow-x: auto; white-space: pre-wrap;
      word-break: break-word; color: #2d3748; margin: 10px 0 20px;
    }}

    /*  PDF embeds  */
    .pdf-block {{ margin: 16px 0 28px; }}
    .pdf-title {{
      font-weight: 600; color: #2d3748; font-size: 0.93em;
      margin-bottom: 6px;
    }}
    object {{
      display: block; border: 1px solid #e2e8f0;
      border-radius: 6px; background: #f7fafc;
    }}

    /*  Note / callout  */
    .note {{
      background: #fffbeb; border-left: 4px solid #f6ad55;
      padding: 12px 16px; border-radius: 0 6px 6px 0;
      font-size: 0.9em; margin: 14px 0;
    }}
    .info {{
      background: #ebf4ff; border-left: 4px solid #4299e1;
      padding: 12px 16px; border-radius: 0 6px 6px 0;
      font-size: 0.9em; margin: 14px 0;
    }}
  </style>
</head>
<body>

<div class="page-header">
  <h1>Phylogenetic Analysis Report</h1>
  <div class="meta">Generated: {now}</div>
  <div class="tags">
    {"".join(f'<span class="tag">{g}</span>' for g in gene_names)}
    <span class="tag">IQ-TREE</span>
    <span class="tag">ASTRAL</span>
    <span class="tag">trimAl</span>
    <span class="tag">MAFFT</span>
  </div>
</div>

<div class="container">

  <!--  Table of Contents  -->
  <nav class="toc">
    <h2>Contents</h2>
    <ol>
      <li><a href="#alignments">Alignment Statistics</a></li>
      <li><a href="#models">Substitution Model Selection</a></li>
      <li><a href="#gene-trees">Individual Gene Trees (Coalescent approach)</a></li>
      <li><a href="#coalescent">Coalescent Species Tree (ASTRAL)</a></li>
      <li><a href="#supermatrix">Supermatrix Tree (IQ-TREE partitioned)</a></li>
      <li><a href="#comparison">Tree Comparison &amp; Robinson-Foulds Distances</a></li>
      <li><a href="#concordance">Topological Concordance Analysis</a></li>
      <li><a href="#discussion">Discussion &amp; Interpretation</a></li>
    </ol>
  </nav>

  <!--  1. Alignment Statistics  -->
  <section id="alignments">
    <h2>1. Alignment Statistics</h2>
    <p>Raw sequences were aligned independently per gene using
       <strong>MAFFT</strong> (--auto) and trimmed with
       <strong>trimAl</strong> (-gappyout). Statistics were computed on the
       raw alignments using <strong>AMAS</strong>.</p>
    {html_table(stats_headers, stats_rows)}
    <div class="note">
      <strong>Note:</strong> Proportion of missing data and parsimony-informative
      sites are key indicators of alignment quality. High missing data (&gt;30%)
      or very low parsimony-informative sites may reduce phylogenetic resolution.
    </div>
  </section>

  <!--  2. Model Selection  -->
  <section id="models">
    <h2>2. Substitution Model Selection</h2>
    <p>Best-fit substitution models were identified per gene using
       <strong>IQ-TREE ModelFinder</strong> (-m MF), which tests a large set
       of nucleotide models and selects the best by BIC.</p>
    {html_table(model_headers, model_rows)}
    <div class="info">
      Different genes may favour different substitution models, reflecting
      differences in evolutionary rates and base composition. Per-gene models
      were used for the individual gene-tree analyses.
    </div>
  </section>

  <!--  3. Individual Gene Trees  -->
  <section id="gene-trees">
    <h2>3. Individual Gene Trees</h2>
    <p>Maximum-likelihood trees were inferred independently for each gene
       using <strong>IQ-TREE</strong> with 100 non-parametric bootstrap
       replicates. The consensus tree with bootstrap support is shown
       (values &ge;50 displayed on nodes).</p>
    {gene_tree_embeds}
  </section>

  <!--  4. Coalescent Species Tree  -->
  <section id="coalescent">
    <h2>4. Coalescent Species Tree (ASTRAL)</h2>
    <p>Individual ML gene trees were combined into a coalescent-based species
       tree using <strong>ASTRAL-III</strong>. Node support is expressed as
       local posterior probability (LPP); values reflect quartet support for
       each bipartition.</p>
    <div class="info">
      The coalescent model explicitly accounts for incomplete lineage sorting
      (ILS). It is preferred when the studied taxa diverged rapidly or when
      effective population sizes were large relative to divergence times.
    </div>
    {embed_pdf(snakemake.input.coalescent_pdf,
               "Coalescent species tree — ASTRAL", height=540)}
  </section>

  <!--  5. Supermatrix Tree  -->
  <section id="supermatrix">
    <h2>5. Supermatrix Tree (Concatenated, IQ-TREE)</h2>
    <p>Trimmed alignments were concatenated into a supermatrix using
       <strong>AMAS</strong>, and a partitioned ML tree was inferred with
       <strong>IQ-TREE</strong> (-m TEST, 100 bootstrap replicates).
       Each gene partition was allowed its own substitution model.</p>
    <div class="note">
      The concatenation approach implicitly assumes a single topology across
      all loci (no ILS). It may produce high bootstrap support even for
      incorrect topologies when there is strong ILS.
    </div>
    {embed_pdf(snakemake.input.concat_pdf,
               "Supermatrix ML tree — IQ-TREE partitioned", height=540)}
  </section>

  <!--  6. Tree Comparison  -->
  <section id="comparison">
    <h2>6. Tree Comparison &amp; Robinson-Foulds Distances</h2>
    <p>Robinson-Foulds (RF) distances (symmetric difference) between all
       tree pairs were calculated using <strong>dendropy</strong>.
       A distance of 0 means identical topologies; larger values indicate
       more discordance.</p>

    <h3>RF Distance Table</h3>
    {html_table(rf_headers, rf_rows)}

    <h3>Topological Summary</h3>
    <pre>{topo_text}</pre>

    <h3>Method Comparison Figure</h3>
    {embed_pdf(snakemake.input.comparison_pdf,
               "Coalescent vs. Supermatrix — side-by-side comparison",
               height=560)}
  </section>

  <!--  7. Concordance  -->
  <section id="concordance">
    <h2>7. Topological Concordance Analysis</h2>
    <p>The proportion of bipartitions in each gene tree that are shared with
       the coalescent and concatenated species trees was calculated to assess
       gene-tree discordance.</p>
    <pre>{concordance_text}</pre>
    <div class="note">
      Low concordance between gene trees and the species tree can indicate
      ILS, hybridisation, or artefacts (alignment errors, long-branch
      attraction). Compare concordance scores with bootstrap support to
      evaluate the robustness of each node.
    </div>
  </section>

  <!--  8. Discussion  -->
  <section id="discussion">
    <h2>8. Discussion &amp; Interpretation</h2>

    <h3>Coalescent vs. Supermatrix</h3>
    <p>Two complementary approaches were applied to infer the species
       phylogeny:</p>
    <ul style="margin: 8px 0 12px 22px; line-height: 1.8;">
      <li><strong>Coalescent (ASTRAL):</strong> statistically consistent
          under the multi-species coalescent; robust to ILS; sensitive to
          gene-tree estimation errors when loci are short.</li>
      <li><strong>Supermatrix (concatenation):</strong> leverages all
          phylogenetic signal jointly; can be misleading when ILS is
          pervasive, producing inflated bootstrap support for incorrect
          relationships.</li>
    </ul>
    <p>Examine the RF distances and concordance table above: if the two
       methods recover the same topology (RF = 0) with high support,
       confidence in the result is elevated. Topological conflict warrants
       closer investigation of gene trees for signatures of ILS or
       horizontal gene transfer.</p>

    <h3>Bootstrap &amp; Local Posterior Probability</h3>
    <p>Non-parametric bootstrap values (IQ-TREE, 100 replicates) and local
       posterior probabilities (ASTRAL) measure support from different
       perspectives. BS values &ge;70 and LPP &ge;0.95 are commonly
       considered well-supported.</p>

    <h3>Recommendations</h3>
    <div class="info">
      Refer to the concordance analysis to identify which nodes are
      consistently recovered across individual gene trees — these are the
      most reliable parts of the species tree regardless of method.
      Discordant nodes should be interpreted cautiously and may require
      additional loci.
    </div>
  </section>

</div><!-- /container -->
</body>
</html>
"""

out_path.write_text(html, encoding="utf-8")
print(f"Report written to: {out_path}")
