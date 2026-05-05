# Inputs from Snakemake
#   gene_trees    – IQ-TREE best ML trees (.treefile) per gene
#   gene_contrees – IQ-TREE bootstrap consensus trees (.contree) per gene
#   coalescent    – ASTRAL species tree (.treefile)
#   concatenated  – IQ-TREE supermatrix best ML tree (.treefile)
#
# Outputs from Snakemake:
#   gene_pdfs      – one PDF per gene (consensus + bootstrap)
#   coalescent_pdf – ASTRAL species tree
#   concat_pdf     – supermatrix ML tree
#   comparison_pdf – side-by-side panel (coalescent vs supermatrix)

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(phangorn)    # midpoint() for tree rooting
  library(patchwork)   # layout side-by-side
})

#  Helpers -----------------------------------

# Extract a display-ready label from an ASTRAL node annotation string.
# ASTRAL-III writes either a plain float ("0.95") or a bracketed string
# "[q1=0.33;q2=0.33;q3=0.34;f1=...;pp1=0.95;...]".
# Returns a formatted string or "" if not parseable.
clean_astral_label <- function(x) {
  if (is.na(x) || nchar(trimws(x)) == 0) return("")
  # Try bracketed format: extract pp1 (local posterior probability)
  if (grepl("\\[", x)) {
    m <- regmatches(x, regexpr("pp1=([0-9.]+)", x))
    if (length(m) > 0 && nchar(m) > 0) {
      val <- as.numeric(sub("pp1=", "", m))
      return(if (is.na(val)) "" else sprintf("%.2f", val))
    }
    return("")
  }
  # Plain numeric
  val <- suppressWarnings(as.numeric(x))
  if (!is.na(val)) return(sprintf("%.2f", val))
  return("")
}

# Read a Newick tree file safely; return NULL on failure.
safe_read_tree <- function(path) {
  tryCatch(read.tree(path), error = function(e) {
    warning("Could not read tree: ", path, "\n  ", e$message)
    NULL
  })
}

# Root a tree at midpoint; silently return as-is if rooting fails
safe_midpoint <- function(tr) {
  tryCatch(phangorn::midpoint(tr), error = function(e) tr)
}

# Root by outgroup if the taxon is present in the tree; fall back to midpoint.
safe_root <- function(tr, outgroup) {
  if (is.null(tr)) return(tr)
  if (!is.na(outgroup) && nchar(trimws(outgroup)) > 0 &&
      outgroup %in% tr$tip.label) {
    tryCatch(
      ape::root(tr, outgroup = outgroup, resolve.root = TRUE),
      error = function(e) {
        warning("Outgroup rooting failed for '", outgroup, "', falling back to midpoint.")
        safe_midpoint(tr)
      }
    )
  } else {
    warning("Outgroup '", outgroup, "' not found in tree tips — using midpoint rooting.")
    safe_midpoint(tr)
  }
}

# Core ggtree plot for a bootstrap-annotated consensus tree.
#   tr         – ape phylo object 
#   title      – plot title string
#   subtitle   – plot subtitle string
#   node_color – colour for bootstrap labels
#   bs_cutoff  – minimum bootstrap to display (default 50)
plot_bs_tree <- function(tr, title, outgroup = NULL,
                         subtitle = "Bootstrap values ≥ 50 shown",
                         node_color = "firebrick3", bs_cutoff = 50) {
  if (is.null(tr)) {
    return(ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Tree unavailable",
               size = 6, color = "grey50") +
      theme_void())
  }

  tr <- safe_root(tr, outgroup)

  # Blank out bootstrap labels below cutoff
  bs <- suppressWarnings(as.numeric(tr$node.label))
  tr$node.label <- ifelse(!is.na(bs) & bs >= bs_cutoff,
                           as.character(round(bs)), "")

  ggtree(tr, branch.length = "branch.length", color = "grey25", linewidth = 0.45) +
    geom_tiplab(size = 2.7, hjust = -0.06, color = "black") +
    geom_nodelab(
      aes(label = label),
      size = 2.1, hjust = 1.15, vjust = -0.35,
      color = node_color, fontface = "bold"
    ) +
    geom_treescale(
      x = 0, y = -0.6,
      fontsize = 2.5, linesize = 0.5, offset = 0.25,
      color = "grey40"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.38))) +
    coord_cartesian(clip = "off") +
    theme_tree2() +
    theme(
      axis.text.x    = element_text(size = 7, color = "grey50"),
      axis.title.x   = element_text(size = 8, color = "grey40"),
      plot.title     = element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle  = element_text(size = 8,  color = "grey50", hjust = 0),
      plot.margin    = margin(6, 10, 6, 6)
    ) +
    labs(title = title, subtitle = subtitle, x = "Substitutions per site")
}

# Variant for ASTRAL coalescent trees
plot_coalescent_tree <- function(tr, title, outgroup = NULL) {
  if (is.null(tr)) {
    return(ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Tree unavailable",
               size = 6, color = "grey50") +
      theme_void())
  }

  tr <- safe_root(tr, outgroup)

  if (!is.null(tr$node.label)) {
    tr$node.label <- sapply(tr$node.label, clean_astral_label)
  }

  ggtree(tr, branch.length = "branch.length", color = "steelblue4", linewidth = 0.45) +
    geom_tiplab(size = 2.7, hjust = -0.06, color = "black") +
    geom_nodelab(
      aes(label = label),
      size = 2.1, hjust = 1.15, vjust = -0.35,
      color = "steelblue4", fontface = "bold"
    ) +
    geom_treescale(
      x = 0, y = -0.6,
      fontsize = 2.5, linesize = 0.5, offset = 0.25,
      color = "grey40"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.38))) +
    coord_cartesian(clip = "off") +
    theme_tree2() +
    theme(
      axis.text.x    = element_text(size = 7, color = "grey50"),
      axis.title.x   = element_text(size = 8, color = "grey40"),
      plot.title     = element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle  = element_text(size = 8,  color = "grey50", hjust = 0),
      plot.margin    = margin(6, 10, 6, 6)
    ) +
    labs(
      title    = title,
      subtitle = "Node labels: local posterior probability (ASTRAL)",
      x        = "Coalescent units"
    )
}

# Snakemake bindings
gene_contrees  <- unlist(snakemake@input[["gene_contrees"]])
gene_treefiles <- unlist(snakemake@input[["gene_trees"]])
coalescent_in  <- snakemake@input[["coalescent"]]
concat_in      <- snakemake@input[["concatenated"]]

gene_pdfs      <- unlist(snakemake@output[["gene_pdfs"]])
coalescent_pdf <- snakemake@output[["coalescent_pdf"]]
concat_pdf     <- snakemake@output[["concat_pdf"]]
comparison_pdf <- snakemake@output[["comparison_pdf"]]

outgroup <- snakemake@params[["outgroup"]]

gene_names <- sub("\\.contree$", "", basename(gene_contrees))

#  Per-gene consensus trees
for (i in seq_along(gene_contrees)) {
  tr <- safe_read_tree(gene_contrees[i])
  p  <- plot_bs_tree(
    tr,
    title    = paste0(gene_names[i], "  —  ML consensus tree"),
    outgroup = outgroup,
    subtitle = "Bootstrap values ≥ 50 shown on internal nodes (100 replicates, IQ-TREE)"
  )
  ggsave(gene_pdfs[i], plot = p, width = 10, height = 8, device = "pdf")
  message("Written: ", gene_pdfs[i])
}

#  ASTRAL coalescent species tree 
tr_coal  <- safe_read_tree(coalescent_in)
p_coal   <- plot_coalescent_tree(tr_coal, "Coalescent species tree  —  ASTRAL", outgroup = outgroup)
ggsave(coalescent_pdf, plot = p_coal, width = 10, height = 8, device = "pdf")
message("Written: ", coalescent_pdf)

# Concatenated supermatrix ML tree
tr_concat <- safe_read_tree(concat_in)
p_concat  <- plot_bs_tree(
  tr_concat,
  title    = "Supermatrix tree  —  IQ-TREE partitioned ML",
  outgroup = outgroup,
  subtitle = "Bootstrap values ≥ 50 shown on internal nodes (100 replicates)"
)
ggsave(concat_pdf, plot = p_concat, width = 10, height = 8, device = "pdf")
message("Written: ", concat_pdf)

#  comparison panel 
p_coal_cmp   <- plot_coalescent_tree(tr_coal,   "Coalescent (ASTRAL)",   outgroup = outgroup)
p_concat_cmp <- plot_bs_tree(tr_concat,
                  title    = "Supermatrix (IQ-TREE)",
                  outgroup = outgroup,
                  subtitle = "Bootstrap ≥ 50 shown")

comparison <- (p_coal_cmp | p_concat_cmp) +
  plot_annotation(
    title    = "Species-tree comparison: coalescent vs. supermatrix",
    subtitle = paste0("Genes analysed: ", paste(gene_names, collapse = ", ")),
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5,
                                   margin = margin(b = 4)),
      plot.subtitle = element_text(size = 10, color = "grey45", hjust = 0.5,
                                   margin = margin(b = 12))
    )
  )

ggsave(comparison_pdf, plot = comparison, width = 18, height = 9, device = "pdf")
message("Written: ", comparison_pdf)

message("visualize_trees.R complete.")
