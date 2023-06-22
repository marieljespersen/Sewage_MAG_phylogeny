# Sewage_MAG_phylogeny

This rep contains scriipts used to make the analysis in the paper: Jespersen et al. "Global phylogeny of metagenome-assembled genomes from sewage suggests that local selection shapes geographical bacterial phylogenetic clustering".

The scripts are arranged in a similar way to the methods section in the paper.



# Data preparing

## genome binning:
```
genome_binning/snake.gs2.vamb.py
```
### checkM:
```
genome_binning/checkm/filter_binsize.sh
genome_binning/checkm/checkm.pbs
```
### dRep:
```
compare: genome_binning/drep/dRep.pbs
dereplicate: genome_binning/drep/dRep_drep_prep.pbs
```
### GTDB_Tk:
```
genome_binning/gtdb_tk/gtdbtk.pbs
```
### Mash: 
```
genome_binning/mash/mash.pbs
```
### CoverM:
```
genome_binning/coverm/snake.coverm.py
genome_binning/coverm/extract_rpkm_ab.R
genome_binning/coverm/expected_coverage.R
```

## Phylogeny:
### ASTRAL:
```
phylogeny/astral/sonicparanoid_to_nucltree.nofilter.pbs
phylogeny/astral/sonic2fasta.all.py
phylogeny/astral/sonicparanoid_prep/Ntogap.sh
phylogeny/astral/phylip.concatenate.py
phylogeny/astral/astral.pbs
phylogeny/astral/clean_tree.sh
phylogeny/astral/unroot_tree.R
phylogeny/astral/gene_iqtrees.pbs
```
### FastTree:
```
phylogeny/fasttree/fasttree.pbs
```

## Functional annotation:
```
functional_annotation/interproscan.sh
functional_annotation/cluster_GO.pbs
functional_annotation/cluster_GO.R
functional_annotation/get_geneids.R
functional_annotation/combine_perm_GO.R
```

## Statistical testing
### PERMANOVA: 
```
statistical_testing/permanova/permanova.genetrees.gs2.R
```
### Gene variance: 
```
statistical_testing/gene_variance/get_gene_div.sh
statistical_testing/gene_variance/seqdiv.py
```
### Gene lengths:
```
statistical_testing/gene_lengths/get_gene_lengths.sh
```
### Regional entropy: 
```
statistical_testing/regional_entropy/get.reg_entropy.R
```


## dN/dS:
### codeml
```
/home/projects/cge/people/maloj/src/dnds/check_start_codon.py
/home/projects/cge/people/maloj/src/dnds/remove_stop_codons.py
/home/projects/cge/data/projects/5001/Binning_vamb/clusters/C3/codeml/codeml.ctl
/home/projects/cge/people/maloj/src/dnds/codeml/prepare_codeml.sh
```
### CSI phylogeny: 
```
SNP pos analyzer: /home/projects/cge/people/maloj/src/C14.dnds.pbs
```

# Plotting 

## Fig 1:
```
map: plotting/worldmap.R
```

## Fig 2:
```
A: plotting/phyla_region_plot.R
B: plotting/plot_oxy.R
C: plotting/PCA.R
D: plotting/diversity.R
```

## Fig 3:
```
A: plotting/plot_trees.R
B: plotting/R2_genomes_pvals_plot.R
C: plotting/ANOVA_R2_allclusters.R
D: plotting/R2_genomes_pvals_plot.R
```

## Supplementary:
```
Fig S1: plotting/plot_metadata.R 
Fig S2: plotting/plot_trees.R
Fig S3: plotting/heatmap_ANI.R
Fig S4: plotting/plot_trees.R
Fig S5: plotting/plot_trees.R
Fig S6: plotting/plot_trees.R
Fig S7: plotting/ANOVA_R2.R
Fig S8: plotting/variance_organelle.R
Fig S9: plotting/variance_organelle.R
Fig S10: plotting/dnds.plot.R
Fig S11: plotting/negative_r2.R
Table S2: plotting/negative_r2.R
```






