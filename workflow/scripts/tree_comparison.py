import os
import itertools
import dendropy
from dendropy.calculate import treecompare

# Load trees

def load_tree(path, taxon_namespace):
    return dendropy.Tree.get(
        path=path,
        schema="newick",
        taxon_namespace=taxon_namespace,
        preserve_underscores=True,
    )


tns = dendropy.TaxonNamespace()

coalescent_tree = load_tree(snakemake.input.coalescent, tns)
concatenated_tree = load_tree(snakemake.input.concatenated, tns)
gene_trees = [load_tree(p, tns) for p in snakemake.input.gene_trees]
gene_names = [
    os.path.basename(p).replace(".treefile", "")
    for p in snakemake.input.gene_trees
]

# Collect all trees for pairwise comparisons
all_trees = gene_trees + [coalescent_tree, concatenated_tree]
all_labels = gene_names + ["coalescent_species_tree", "concatenated_supermatrix"]

# Unroot all trees for RF distance
for t in all_trees:
    t.is_rooted = False
    t.collapse_basal_bifurcation(set_as_unrooted_tree=True)

# RF distances

rf_records = []
for (i, label_i), (j, label_j) in itertools.combinations(
    enumerate(all_labels), 2
):
    rf = treecompare.symmetric_difference(all_trees[i], all_trees[j])
    rf_records.append((label_i, label_j, rf))

with open(snakemake.output.rf, "w") as out:
    out.write("Robinson-Foulds (symmetric difference) distances\n")
    out.write("=" * 60 + "\n")
    out.write(f"{'Tree A':<35} {'Tree B':<35} {'RF'}\n")
    out.write("-" * 60 + "\n")
    for a, b, rf in rf_records:
        out.write(f"{a:<35} {b:<35} {rf}\n")

# Summary

gene_vs_coal = [
    treecompare.symmetric_difference(gt, coalescent_tree) for gt in gene_trees
]
gene_vs_concat = [
    treecompare.symmetric_difference(gt, concatenated_tree) for gt in gene_trees
]
coal_vs_concat = treecompare.symmetric_difference(coalescent_tree, concatenated_tree)

with open(snakemake.output.topo, "w") as out:
    out.write("Topological Summary\n")
    out.write("=" * 60 + "\n\n")

    out.write("RF distances: individual gene trees vs. coalescent species tree\n")
    out.write("-" * 50 + "\n")
    for name, rf in zip(gene_names, gene_vs_coal):
        out.write(f"  {name}: {rf}\n")
    if gene_vs_coal:
        out.write(f"  Mean: {sum(gene_vs_coal)/len(gene_vs_coal):.2f}\n\n")

    out.write("RF distances: individual gene trees vs. concatenated tree\n")
    out.write("-" * 50 + "\n")
    for name, rf in zip(gene_names, gene_vs_concat):
        out.write(f"  {name}: {rf}\n")
    if gene_vs_concat:
        out.write(f"  Mean: {sum(gene_vs_concat)/len(gene_vs_concat):.2f}\n\n")

    out.write("RF distance: coalescent vs. concatenated species tree\n")
    out.write("-" * 50 + "\n")
    out.write(f"  RF = {coal_vs_concat}\n")

# Concordance analysis

def concordant_bipartitions(ref_tree, query_trees, labels):
    ref_bip = ref_tree.encode_bipartitions()
    ref_set = {b.split_bitmask for b in ref_tree.bipartition_encoding}
    results = []
    for qt, lbl in zip(query_trees, labels):
        qt.encode_bipartitions()
        q_set = {b.split_bitmask for b in qt.bipartition_encoding}
        shared = len(ref_set & q_set)
        total = len(ref_set)
        results.append((lbl, shared, total))
    return results


coal_concordance = concordant_bipartitions(coalescent_tree, gene_trees, gene_names)
concat_concordance = concordant_bipartitions(concatenated_tree, gene_trees, gene_names)

with open(snakemake.output.concordance, "w") as out:
    out.write("Bipartition Concordance Analysis\n")
    out.write("=" * 60 + "\n\n")

    out.write("Gene tree concordance with coalescent species tree\n")
    out.write("-" * 50 + "\n")
    out.write(f"{'Gene':<20} {'Shared bipartitions':<22} {'Reference total'}\n")
    for lbl, shared, total in coal_concordance:
        pct = 100 * shared / total if total else 0
        out.write(f"  {lbl:<18} {shared}/{total} ({pct:.1f}%)\n")

    out.write("\nGene tree concordance with concatenated tree\n")
    out.write("-" * 50 + "\n")
    out.write(f"{'Gene':<20} {'Shared bipartitions':<22} {'Reference total'}\n")
    for lbl, shared, total in concat_concordance:
        pct = 100 * shared / total if total else 0
        out.write(f"  {lbl:<18} {shared}/{total} ({pct:.1f}%)\n")

    out.write(
        "\nNote: bipartition counts exclude the root/trivial bipartitions.\n"
    )
