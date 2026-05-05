# Phylogenetic analysis report

## Dataset

I used three genes from 24 mammalian taxa from a bat and whales project I was recently working on: RAG1 and RAG2 (nuclear) and CYTB (mitochondrial). The taxa spans cetaceans, a diverse chiropteran sample pf nine bat species across 5 families, and 4 outgroup placental mammals: oryx (*Oryx dammah*), groundhog (*Marmota monax*), human (*Homo sapiens*), and blue-eyed black lemur (*Eulemur flavifrons*). RAG2 sequences were missing for three taxa relative to RAG1 and CYTB.

Sequences were downloaded from NCBI and placed as unaligned FASTA files in `data/`.

## Alignment quality

After MAFFT alignment and trimAl trimming, the three loci differed immensely in quality:

| Gene | Taxa | Aligned length (bp) | Missing data (%) | Parsimony-informative sites |
|------|------|--------------------|--------------------|----------------------------|
| RAG1 | 24 | 74,831 | 92.5 | 672 (0.9%) |
| RAG2 | 21 | 8,144 | 59.3 | 364 (4.5%) |
| CYTB | 24 | 1,141 | 9.1 | 507 (44.4%) |

RAG1's alignment has 92.5% missing data with only 672 parsimony-informative sites out of ~75,000 positions. This might be because the sequences in GenBank for different taxa cover different exon regions, so large regions are represented by only a few taxa. CYTB is the best marker: short, nearly complete across all taxa, and 44% of sites are parsimony-informative. RAG2 is between the two.

## Substitution models

ModelFinder in IQ-TREE `-m MF`, selected a different model for each gene, as expected:

| Gene | Best model | lnL | AIC | BIC |
|------|-----------|-----|-----|-----|
| RAG1 | TPM3u+R2 | −12,762.4 | 25,622.8 | 25,912.4 |
| RAG2 | TPM3u+F+G4 | −16,227.7 | 32,545.4 | 32,848.6 |
| CYTB | TIM2+F+R3 | −11,392.4 | 22,894.8 | 23,171.8 |

CYTB got the most parameter-rich model whereas RAG1 and RAG2 resolved to TPM3u-based models with different rate treatments: +R2 and +G4, respectively.

## Gene trees

Individual Maximum likelihood (ML) trees were inferred with IQ-TREE with 100 bootstrap replicates, per-gene model. The RAG1 tree and the CYTB tree recovered the major mammalian groupings mostly with low to medium bootstrap support (> 50 & < 80), and major super families with high support (100). The long branch of *Myotis myotis* in the RAG1 tree is a major concern due to long branch attraction. The RAG2 tree resolved better the cetacean clade and also deeper nodes, despite the 59% missing data and lower phylogenetic signal per site relative to CYTB. This discordance is not surprising and is one of the reasons why resolving the Laurasiatherian tree is so difficult.

## Species tree with ASTRAL

The three ML gene trees were used as input for ASTRAL-III. ASTRAL uses quartet frequencies rather than branch lengths, so the missing-data problem in RAG1 is less damaging than it would be in a likelihood framework.

The coalescent tree places the mysticetes (gray whale + blue whale) sister to a clade containing bats, primates, artiodactyls, and rodents, with odontocetes occupying several positions elsewhere in the tree. Node support (local posterior probability in this case) is uniformly 0.67 across most internal nodes, which reflects limited signal from only three gene trees. ASTRAL cannot distinguish between the three possible resolutions of any given quartet with high confidence due to the low number of loci. A larger gene set could improve support, yet again, this tree has been historically hard to resolve even with hundreds of loci, so the low support is not surprising.

## Species tree: supermatrix approach (IQ-TREE, partitioned)

The three trimmed alignments were concatenated with AMAS into a 84,116 bp supermatrix. IQ-TREE inferred a partitioned ML tree (`-m TEST`, per-partition models, 100 bootstrap replicates) treating each gene as an independent partition.

## Tree comparison

Robinson-Foulds distances between all tree pairs:

| Tree A | Tree B | RF |
|--------|--------|----|
| RAG1 | RAG2 | 34 |
| RAG1 | CYTB | 18 |
| RAG1 | coalescent | 0 |
| RAG1 | concatenated | 10 |
| RAG2 | CYTB | 28 |
| RAG2 | coalescent | 34 |
| RAG2 | concatenated | 26 |
| CYTB | coalescent | 18 |
| CYTB | concatenated | 10 |
| coalescent | concatenated | 10 |

The most striking and concerning result is that RAG1's tree is topologically identical to the ASTRAL coalescent tree (RF = 0), while RAG2 disagrees with both species trees substantially (RF = 34 and 26). The concatenated tree sits between the two, and is different from the coalescent tree by RF = 10.

## Concordance analysis

Bipartitions in each gene tree that are shared with the two species trees:

| Gene | Concordance with coalescent | Concordance with concatenated |
|------|----------------------------|-------------------------------|
| RAG1 | 46/46 (100%) | 41/46 (89%) |
| RAG2 | 26/46 (57%) | 30/46 (65%) |
| CYTB | 37/46 (80%) | 41/46 (89%) |

RAG1 is perfectly concordant with the coalescent tree, meaning, every bipartition in the coalescent topology appears in the RAG1 gene tree. This explains the RF = 0 result above. RAG2 is the outlier since it disagrees with both species trees at 35–43% of bipartitions. CYTB agrees with the concatenated tree as strongly as RAG1 does (89%), but more poorly with the coalescent tree (80%).

The high RAG2 discordance deserves scrutiny. Incomplete lineage sorting is one explanation, but the 59% missing data makes it difficult to rule out alignment artifacts. The three Rhinolophid + Pteropodid bat relationships were most affected; those clades were also the ones with the most fragmented RAG2 coverage in GenBank.

## Method comparison

The coalescent and concatenated trees agree on most major groupings (RF = 10), and were only different at five bipartitions. Their mean gene-tree RF distances are similar, 17.3 for the coalescent tree and 15.3 for the concatenated tree. On average, neither method is closer to the gene trees.

For this dataset, the coalescent approach is conceptually preferable because the taxon set includes both rapid radiations (dolphins) and deep divergences (cetaceans vs. bats), a combination where incomplete lineage sorting is expected at shallow nodes. The uniform ASTRAL support (all 0.67) is a practical problem. However with three genes, ASTRAL lacks the resolution to distinguish high quality quartet topologies. More loci might be needed before the coalescent tree can be trusted over the concatenated one at most nodes.

The concatenated tree has better apparent support, but this conflates signal from mitochondrial and nuclear genomes. CYTB evolves under maternal inheritance, so forcing it into a single partition model can exaggerate support for the wrong topology when nuclear and mitochondrial histories disagree.

I would use the concatenated ML tree for downstream analyses as a working hypothesis while treating any node that disagrees with the coalescent tree.

## Challenges and solutions

Once again, most challenges were associated with software and data quality rather than the methods themselves. Although my knowledge on ASTRAL is a little rusty, I had no problem in setting up the parameters and data to run it.

Major problems included **R conda environment conflicts**. This is a known problem in our cluster. So mixing the `defaults` channel with `conda-forge` and `bioconda` made the solver unable to resolve `libdeflate` versions. Removing `- defaults` from all environment YAML files resolved every conflict. Also, `r-phytools` depends on `r-rgl` which had no working build for linux-64 on the cluster.

Also had problems with **trimAl `-automated1` crashed** on sequences with lowercase ambiguous bases (`n`), and solved it by switching it to `-gappyout` to resolve this. Also to get the best fit models in IQTREE, I had to fix two issues with the IQ-TREE report parsing.

Finally, I couldn't get the **Partition file** correctly since AMAS writes partitions as `p1_RAG1 = 1-2724` but IQ-TREE's RAxML partition format requires `DNA, RAG1 = 1-2724`. I included a one-line `sed` call inside the `concatenated_tree` rule to convert the format before IQ-TREE ran.

