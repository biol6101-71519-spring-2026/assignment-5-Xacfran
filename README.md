[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/NNshSCAo)
# Assignment 5: Phylogenetic Analysis Pipeline

## Overview
In this assignment, you will build a **Snakemake workflow** for phylogenetic analysis. You will start with a set of homologous gene sequences, perform multiple sequence alignment, trim and filter the alignments, and construct phylogenetic trees using two different approaches: the **species coalescent method** and the **supermatrix (concatenation) method**. You will assess tree support through bootstrap analysis, compare the results from both approaches, and visualize your findings. This assignment combines computational phylogenetics with reproducible workflow development and critical evaluation of phylogenetic methods.

## Learning Objectives
By completing this assignment, you will:
- Perform multiple sequence alignment on homologous sequences
- Trim and filter alignments
- Construct phylogenetic trees using maximum likelihood methods
- Compare species coalescent vs. supermatrix approaches to phylogeny
- Understand the assumptions and limitations of each phylogenetic method
- Evaluate tree reliability through bootstrap and coalescent support values
- Select appropriate evolutionary models for phylogenetic inference
- Assess concordance and discordance among gene trees
- Visualize and annotate phylogenetic trees
- Build reproducible phylogenetic workflows with Snakemake

### Sequence alignment tools:
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) (multiple alignment)
- [trimAl](http://trimal.cgenomics.org/) (alignment trimming)
- [BMGE](https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-10-210) (alignment filtering)

### Phylogenetic inference tools:
- [IQ-TREE](http://www.iqtree.org/) (maximum likelihood + model selection)
- [RAxML-NG](https://github.com/amkozlov/raxml-ng) (maximum likelihood)
- [ASTRAL](https://github.com/smirarab/ASTRAL) (coalescent-based species tree estimation)

### Concatenation tools:
- [AMAS](https://github.com/marekborowiec/AMAS) (alignment manipulation and statistics)
- [catfasta2phyml](https://github.com/nylander/catfasta2phyml) (concatenate FASTA alignments)

### Tree analysis and visualization:
- [ETE3](http://etetoolkit.org/) (Python tree manipulation)
- R with packages: ape, phytools, ggtree, treeio, phangorn
- [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) (optional, for manual viewing)

**Important:** Use **per-rule conda environments** or **Apptainer containers** (stored in `workflow/envs/`) to ensure reproducibility and proper dependency management. You must create separate environment YAML files for each tool/rule.

Example base installation:
```bash
# Create a base environment with Snakemake if you don't have it already
mamba create -n snakemake -c bioconda -c conda-forge snakemake

# Create individual tool environment files in workflow/envs/
# Example: workflow/envs/mafft.yaml, workflow/envs/iqtree.yaml, etc.

# Alternatively, use Apptainer/Singularity containers
# Snakemake can automatically pull containers from conda packages
# Run with: snakemake --use-conda --use-apptainer
```

## Data
For this assignment, you need to provide your own homologous gene sequences for phylogenetic analysis. Download or prepare sequences from public databases or your own research.

**Data sources:**
- **NCBI**: Download gene sequences using Entrez Direct or the web interface
- **Ensembl**: Obtain gene sequences from model organisms
- **OrthoDB**: Download orthologous gene sets
- **OrthoFinder output**: Single-copy orthologous genes from genome comparisons

**Requirements:**
- Each FASTA file should contain sequences from the same gene/locus across different taxa
- Use 3-5 genes for a multi-locus analysis
- Start with 10-20 taxa (species/strains)
- Each gene alignment should be 500-2000 bp

**FASTA format example:**
```
>Species_A
ATCGATCG...
>Species_B
ATCGATGG...
```

## Assignment Tasks

### Task 1: Multiple Sequence Alignment
Implement Snakemake rules to:
1. Align each gene/locus independently using **MAFFT**
2. Compute alignment statistics using **AMAS**
3. Trim/filter poorly aligned regions using **trimAl**
4. Generate a report on alignment quality for each gene

#### Expected outputs:
- `results/alignments/raw/{gene}.aln` (raw alignments)
- `results/alignments/trimmed/{gene}.trimmed.aln` (filtered alignments)
- `results/alignment_stats/{gene}_stats.txt`
- `results/alignment_stats/summary_report.txt`

### Task 2: Evolutionary Model Selection
Implement Snakemake rules to:
1. Test substitution models for each gene using **IQ-TREE** or **ModelTest-NG**
2. Identify the best-fit model
3. Summarize model selection results across all genes

#### Expected outputs:
- `results/model_selection/{gene}.model`
- `results/model_selection/best_models_summary.txt`

**Note:** Different genes may have different evolutionary models.

### Task 3: Phylogenetic Tree Inference
Build phylogenetic trees using two approaches to compare methodological assumptions and results.

#### Task 3a: Species Coalescent Approach
Implement Snakemake rules to:
1. Construct ML trees for each gene independently using **IQ-TREE**
2. Perform bootstrap analysis (100 real bootstrap replicates)
3. Generate trees with bootstrap support values for each gene
4. Estimate a species tree from individual gene trees using a coalescent-based method (e.g., **ASTRAL**)

##### Expected outputs:
- `results/trees/individual/{gene}.treefile` (best ML tree per gene)
- `results/trees/individual/{gene}.contree` (consensus tree with bootstrap support)
- `results/trees/individual/{gene}.log` (log file with likelihood scores)
- `results/trees/coalescent/species_tree.treefile` (coalescent-based species tree)

#### Task 3b: Supermatrix Approach
Implement Snakemake rules to:
1. Concatenate alignments from multiple genes into a supermatrix
2. Create a partition file defining gene boundaries
3. Construct a phylogenetic tree from the concatenated alignment using **IQ-TREE**
4. Perform bootstrap analysis on the concatenated alignment

##### Expected outputs:
- `results/concatenated/concatenated.aln` (supermatrix alignment)
- `results/concatenated/partition_info.txt` (partition boundaries for each gene)
- `results/trees/concatenated/concatenated.treefile` (ML tree from supermatrix; name may vary based on tool)
- `results/trees/concatenated/concatenated.support` (ML tree with bootstrap support; name may vary based on tool)

### Task 4: Tree Comparison and Visualization
Implement Snakemake rules to:
1. Compare the species coalescent tree (Task 3a) with the supermatrix tree (Task 3b)
2. Compare individual gene trees with both species trees
3. Calculate Robinson-Foulds distances between all tree pairs
4. Assess topological concordance across genes and methods
5. Visualize trees with support values using **R (ggtree)** or **Python (ETE3)**
6. Create comparison plots showing both approaches side-by-side

#### Expected outputs:
- `results/visualizations/{gene}_tree.pdf` or `.png` (individual gene trees)
- `results/visualizations/coalescent_tree.pdf` (species coalescent tree)
- `results/visualizations/concatenated_tree.pdf` (supermatrix tree)
- `results/visualizations/method_comparison.pdf` (both approaches side-by-side)
- `results/tree_comparison/robinson_foulds_distances.txt`
- `results/tree_comparison/topological_summary.txt`
- `results/tree_comparison/concordance_analysis.txt`

**Visualization should include:**
- Branch lengths
- Bootstrap/support values
- Species labels
- Scale bar
- Optional: tip labels colored by taxonomy/function

### Task 5: Phylogenetic Analysis Report
Create analysis scripts in `workflow/scripts/` and generate a report that includes:
1. Comparison of species coalescent vs. supermatrix approaches
2. Assessment of node support (bootstrap/coalescent support values)
3. Discussion of any discordance between gene trees
4. Evaluation of which approach is more appropriate for your dataset
5. Final species tree recommendation with justification

**You must write your own scripts for:**
- Tree comparison (calculating Robinson-Foulds distances)
- Coalescent tree inference (e.g., using ASTRAL)
- Tree visualization (using ggtree, ETE3, or similar)
- Statistical summaries and concordance analysis
- Report generation

#### Expected output:
- `results/report/phylogenetic_analysis_report.pdf` or `.html`

## Workflow Structure
Your final directory structure should look like:
```
.
├── README.md
├── data/                 # Place your input FASTA files here
├── workflow/
│   ├── Snakefile
│   ├── scripts/          # Create your own analysis scripts here
│   │   ├── [your scripts for tree comparison, visualization, etc.]
│   └── envs/             # Create conda environment files for each tool
│       ├── [your environment YAML files]
├── config/
│   └── config.yaml
└── results/              # Generated by workflow
    ├── alignments/
    │   ├── raw/
    │   ├── trimmed/
    │   └── concatenated/     # Task 3b: concatenated alignment
    ├── alignment_stats/
    ├── model_selection/
    ├── trees/
    │   ├── individual/       # Task 3a: per-gene trees
    │   ├── coalescent/       # Task 3a: species coalescent tree
    │   └── concatenated/     # Task 3b: supermatrix tree
    ├── visualizations/
    ├── tree_comparison/
    └── report/
```

**Note:** 
- Place your input sequence FASTA files in the `data/` directory and specify their location in `config/config.yaml`
- You must create all environment YAML files and analysis scripts yourself
- **AI Usage**: You may use AI assistants (e.g., ChatGPT, GitHub Copilot) to help write analysis scripts (in `workflow/scripts/`), but the Snakemake workflow itself (`Snakefile`) must be your own work

## Running Your Workflow
Once your workflow is complete, run it with:
```bash
# Dry run to check workflow
snakemake -n

# Run workflow with 4 cores using conda
snakemake --cores 4 --use-conda

# Run workflow with Apptainer containers (alternative to conda)
snakemake --cores 4 --use-conda --use-apptainer

# Generate workflow diagram
snakemake --dag | dot -Tpng > workflow_dag.png
```

## Tips for Success
- **Start small**: Test your workflow on a subset of genes before running all
- **Validate alignments**: Visually inspect a few alignments to ensure quality
- **Monitor resources**: Bootstrap analyses can be computationally intensive
- **Understand the methods**: Read about the assumptions of coalescent vs. concatenation approaches
- **Check for discordance**: Examine individual gene trees for signals of incomplete lineage sorting
- **Create conda environments**: Write YAML files for each tool specifying channels and dependencies (or use Apptainer)
- **Write helper scripts**: Create scripts for tree comparison, visualization, and statistics
- **Document decisions**: Note why you chose specific parameters or models
- **Version control**: Use Git to track changes to your workflow
- **Read documentation**: Consult tool manuals and Snakemake docs for proper syntax

## Additional Resources
- [IQ-TREE tutorial](http://www.iqtree.org/doc/)
- [ASTRAL tutorial](https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md)
- [Species tree inference: a guide to methods and their application](https://www.nature.com/articles/s41559-017-0126) (review paper)
- [ggtree documentation](https://yulab-smu.top/treedata-book/)
- [ETE3 tutorial](http://etetoolkit.org/docs/latest/tutorial/index.html)
- [Snakemake phylogenetics examples](https://github.com/snakemake-workflows)
- [MAFFT manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html)
- [Concatenation vs. coalescence debate](https://doi.org/10.1146/annurev-ecolsys-120213-091627)

## Submission
Submit the following:
1. Complete workflow directory (including `workflow/Snakefile`, `config/config.yaml`, scripts, and environment files)
2. A brief report (`report.md` or `.pdf`) with:
   - Description of your dataset
   - Summary of key results from both phylogenetic approaches
   - Comparison of species coalescent vs. supermatrix methods
   - Discussion of concordance/discordance among gene trees
   - Evaluation of which method is more appropriate for your data and why
   - Any challenges encountered and how you solved them
3. Example output files (at least one tree visualization from each approach)

**Do not include** raw data files or large result files in your submission (provide download links if needed).