---
description: Identify genes that are significantly over or
  under-expressed between conditions in specific cell populations.
subtitle:  Scanpy Toolkit
title:  Differential gene expression
---
<div>

> **Note**
>
> Code chunks run Python commands unless it starts with `%%bash`, in
> which case, those chunks run shell commands.

</div>

In this tutorial we will cover differential gene expression, which
comprises an extensive range of topics and methods. In single cell,
differential expresison can have multiple functionalities such as
identifying marker genes for cell populations, as well as identifying
differentially regulated genes across conditions (healthy vs control).
We will also cover controlling batch effect in your test.

Differential expression is performed with the function rank_genes_group.
The default method to compute differential expression is the
t-test_overestim_var. Other implemented methods are: logreg, t-test and
wilcoxon.

By default, the .raw attribute of AnnData is used in case it has been
initialized, it can be changed by setting use_raw=False.

The clustering with resolution 0.6 seems to give a reasonable number of
clusters, so we will use that clustering for all DE tests.

First, let's import libraries and fetch the clustered data from the
previous lab.

Read in the clustered data object.




    AnnData object with n_obs × n_vars = 7332 × 3984
        obs: 'type', 'sample', 'batch', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb', 'pct_counts_hb', 'percent_mt2', 'n_counts', 'n_genes', 'percent_chrY', 'XIST-counts', 'S_score', 'G2M_score', 'phase', 'doublet_scores', 'predicted_doublets', 'doublet_info', 'leiden', 'leiden_0.4', 'leiden_0.6', 'leiden_1.0', 'leiden_1.4', 'kmeans5', 'kmeans10', 'kmeans15', 'hclust_5', 'hclust_10', 'hclust_15'
        var: 'gene_ids', 'feature_types', 'genome', 'mt', 'ribo', 'hb', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts', 'n_cells', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'mean', 'std'
        uns: 'dendrogram_leiden_0.6', 'doublet_info_colors', 'hclust_10_colors', 'hclust_15_colors', 'hclust_5_colors', 'hvg', 'kmeans10_colors', 'kmeans15_colors', 'kmeans5_colors', 'leiden', 'leiden_0.4', 'leiden_0.4_colors', 'leiden_0.6', 'leiden_0.6_colors', 'leiden_1.0', 'leiden_1.0_colors', 'leiden_1.4', 'leiden_1.4_colors', 'log1p', 'neighbors', 'pca', 'phase_colors', 'sample_colors', 'tsne', 'umap'
        obsm: 'Scanorama', 'X_pca', 'X_pca_combat', 'X_pca_harmony', 'X_tsne', 'X_tsne_bbknn', 'X_tsne_combat', 'X_tsne_harmony', 'X_tsne_scanorama', 'X_tsne_uncorr', 'X_umap', 'X_umap_bbknn', 'X_umap_combat', 'X_umap_harmony', 'X_umap_scanorama', 'X_umap_uncorr'
        varm: 'PCs'
        obsp: 'connectivities', 'distances'



Check what you have in the different matrices.

    (7332, 3984)
    <class 'anndata._core.raw.Raw'>
    [[-0.17232972 -0.54321549 -0.31634794 -0.90417343 -0.16910255 -0.05583347
      -0.25571551 -0.20106716 -0.07050795 -0.35735729]
     [ 5.54532801 -0.54321549 -0.31634794 -0.90417343 -0.16910255 -0.05583347
      -0.25571551 -0.20106716 -0.07050795 -0.35735729]
     [-0.17232972 -0.54321549 -0.31634794  0.17987218 -0.16910255 -0.05583347
      -0.25571551 -0.20106716 -0.07050795 -0.35735729]
     [-0.17232972 -0.54321549 -0.31634794 -0.90417343 -0.16910255 -0.05583347
      -0.25571551  6.01645508 -0.07050795 -0.35735729]
     [-0.17232972 -0.54321549 -0.31634794 -0.90417343 -0.16910255 -0.05583347
      -0.25571551 -0.20106716 -0.07050795 -0.35735729]
     [-0.17232972 -0.54321549 -0.31634794 -0.90417343 -0.16910255 -0.05583347
      -0.25571551 -0.20106716 -0.07050795 -0.35735729]
     [-0.17232972 -0.54321549 -0.31634794 -0.90417343 -0.16910255 -0.05583347
      -0.25571551 -0.20106716 -0.07050795 -0.35735729]
     [-0.17232972 -0.54321549 -0.31634794 -0.90417343 -0.16910255 -0.05583347
      -0.25571551 -0.20106716 -0.07050795 -0.35735729]
     [-0.17232972  1.93648021 -0.31634794 -0.90417343 -0.16910255 -0.05583347
      -0.25571551 -0.20106716 -0.07050795 -0.35735729]
     [-0.17232972 -0.54321549 -0.31634794 -0.90417343 -0.16910255 -0.05583347
      -0.25571551 -0.20106716 -0.07050795 -0.35735729]]


As you can see, the X matrix contains all genes and the data looks
logtransformed.

For DGE analysis we would like to run with all genes, on normalized
values, so if you did subset the `adata.X` for variable genes you would
have to revert back to the raw matrix with
`adata = adata.raw.to_adata()`. In case you have raw counts in the
matrix you also have to renormalize and logtransform.

Now lets look at the clustering of the object we loaded in the umap. We
will use leiden_0.6 clustering in this exercise. If you recall from the
previous exercise, we set the default umap to the umap created with
Harmony.


    
![png](scanpy_05_dge_files/scanpy_05_dge_8_0.png)
    


## T-test

    ranking genes
        finished (0:00:00)



    
![png](scanpy_05_dge_files/scanpy_05_dge_10_1.png)
    





    dict_keys(['dendrogram_leiden_0.6', 'doublet_info_colors', 'hclust_10_colors', 'hclust_15_colors', 'hclust_5_colors', 'hvg', 'kmeans10_colors', 'kmeans15_colors', 'kmeans5_colors', 'leiden', 'leiden_0.4', 'leiden_0.4_colors', 'leiden_0.6', 'leiden_0.6_colors', 'leiden_1.0', 'leiden_1.0_colors', 'leiden_1.4', 'leiden_1.4_colors', 'log1p', 'neighbors', 'pca', 'phase_colors', 'sample_colors', 'tsne', 'umap', 't-test'])



## T-test overestimated_variance

    ranking genes
        finished (0:00:00)



    
![png](scanpy_05_dge_files/scanpy_05_dge_12_1.png)
    


## Wilcoxon rank-sum

The result of a Wilcoxon rank-sum (Mann-Whitney-U) test is very similar.
We recommend using the latter in publications, see e.g., Sonison &
Robinson (2018). You might also consider much more powerful differential
testing packages like MAST, limma, DESeq2 and, for python, the recent
diffxpy.

    ranking genes
        finished (0:00:03)



    
![png](scanpy_05_dge_files/scanpy_05_dge_14_1.png)
    


## Logistic regression test

As an alternative, let us rank genes using logistic regression. For
instance, this has been suggested by Natranos et al. (2018). The
essential difference is that here, we use a multi-variate appraoch
whereas conventional differential tests are uni-variate. Clark et
al. (2014) has more details.

    ranking genes
        finished (0:00:07)



    
![png](scanpy_05_dge_files/scanpy_05_dge_16_1.png)
    


## Compare genes

Take all significant DE genes for cluster0 with each test and compare
the overlap.


    
![png](scanpy_05_dge_files/scanpy_05_dge_18_0.png)
    


As you can see, the Wilcoxon test and the T-test with overestimated
variance gives very similar result. Also the regular T-test has good
overlap.

## Visualization

There are several ways to visualize the expression of top DE genes. Here
we will plot top 5 genes per cluster from Wilcoxon test as heatmap,
dotplot, violin plots or a matrix with average expression.


    
![png](scanpy_05_dge_files/scanpy_05_dge_20_0.png)
    



    
![png](scanpy_05_dge_files/scanpy_05_dge_20_1.png)
    



    
![png](scanpy_05_dge_files/scanpy_05_dge_20_2.png)
    



    
![png](scanpy_05_dge_files/scanpy_05_dge_20_3.png)
    


## Compare specific clusters

We can also do pairwise comparisons of individual clusters on one vs
many clusters. For instance, clusters 1 & 2 have very similar expression
profiles.

    ranking genes
        finished (0:00:00)



    
![png](scanpy_05_dge_files/scanpy_05_dge_22_1.png)
    


Plot as violins for those two groups, or across all the clusters.


    
![png](scanpy_05_dge_files/scanpy_05_dge_24_0.png)
    



    
![png](scanpy_05_dge_files/scanpy_05_dge_24_1.png)
    


## DGE across conditions

The second way of computing differential expression is to answer which
genes are differentially expressed within a cluster. For example, in our
case we have libraries comming from patients and controls and we would
like to know which genes are influenced the most in a particular cell
type. For this end, we will first subset our data for the desired cell
cluster, then change the cell identities to the variable of comparison
(which now in our case is the **type**, e.g. Covid/Ctrl).

    ranking genes
        finished (0:00:00)



    
![png](scanpy_05_dge_files/scanpy_05_dge_26_1.png)
    



    
![png](scanpy_05_dge_files/scanpy_05_dge_27_0.png)
    



    
![png](scanpy_05_dge_files/scanpy_05_dge_27_1.png)
    


We can also plot these genes across all clusters, but split by "type",
to check if the genes are also up/downregulated in other celltypes.




    <seaborn.axisgrid.FacetGrid at 0x14f483a10>




    
![png](scanpy_05_dge_files/scanpy_05_dge_29_1.png)
    


As you can see, we have many sex chromosome related genes among the top
DE genes. And if you remember from the QC lab, we have inbalanced sex
distribution among our subjects, so this is probably not related to
covid at all.

### Remove sex chromosome genes

To remove some of the bias due to inbalanced sex in the subjects we can
remove the sex chromosome related genes.

    134
    3984
    3850


Rerun differential expression.

    ranking genes
        finished (0:00:00)



    
![png](scanpy_05_dge_files/scanpy_05_dge_33_1.png)
    


Now at least we do not have the sex chromosome genes as DE but still,
some of the differences between patient and control could still be
related to sex.

### Patient batch effects

When we are testing for Covid vs Control we are running a DGE test for 4
vs 4 individuals. That will be very sensitive to sample differences
unless we find a way to control for it. So first, lets check how the top
DGEs are expressed in that cluster, across the individuals:


    
![png](scanpy_05_dge_files/scanpy_05_dge_35_0.png)
    



    
![png](scanpy_05_dge_files/scanpy_05_dge_35_1.png)
    


As you can see, many of the genes detected as DGE in Covid are unique to
one or 2 patients.

We can also plot the top Covid and top Ctrl genes as a dotplot:


    
![png](scanpy_05_dge_files/scanpy_05_dge_37_0.png)
    


Clearly many of the top Covid genes are only high in the covid_17
sample, and not a general feature of covid patients.

This is also the patient with the highest number of cells in this
cluster:




    sample
    covid_17    181
    ctrl_5      153
    covid_1      96
    ctrl_13      68
    ctrl_14      61
    ctrl_19      42
    covid_16     41
    covid_15     32
    Name: count, dtype: int64



### Subsample

So one obvious thing to consider is an equal amount of cells per
individual so that the DGE results are not dominated by a single sample.

So we will downsample to an equal number of cells per sample, in this
case 34 cells per sample as it is the lowest number among all samples




    sample
    covid_1     37
    covid_16    37
    covid_17    37
    ctrl_5      37
    ctrl_13     37
    ctrl_14     37
    ctrl_19     37
    covid_15    32
    Name: count, dtype: int64



    ranking genes
        finished (0:00:00)



    
![png](scanpy_05_dge_files/scanpy_05_dge_42_1.png)
    



    
![png](scanpy_05_dge_files/scanpy_05_dge_43_0.png)
    


It looks much better now. But if we look per subject you can see that we
still have some genes that are dominated by a single patient. Still, it
is often a good idea to control the number of cells from each sample
when doing differential expression.

There are many different ways to try and resolve the issue of patient
batch effects, however most of them require R packages. These can be run
via rpy2 as is demonstraded in this compendium:
https://www.sc-best-practices.org/conditions/differential_gene_expression.html

However, we have not included it here as of now. So please have a look
at the patient batch effect section in the seurat DGE tutorial where we
run EdgeR on pseudobulk and MAST with random effect:
https://nbisweden.github.io/workshop-scRNAseq/labs/seurat/seurat_05_dge.html

## Gene Set Analysis (GSA)

### Hypergeometric enrichment test

Having a defined list of differentially expressed genes, you can now
look for their combined function using hypergeometric test.

    ['ARCHS4_Cell-lines', 'ARCHS4_IDG_Coexp', 'ARCHS4_Kinases_Coexp', 'ARCHS4_TFs_Coexp', 'ARCHS4_Tissues', 'Achilles_fitness_decrease', 'Achilles_fitness_increase', 'Aging_Perturbations_from_GEO_down', 'Aging_Perturbations_from_GEO_up', 'Allen_Brain_Atlas_10x_scRNA_2021', 'Allen_Brain_Atlas_down', 'Allen_Brain_Atlas_up', 'Azimuth_2023', 'Azimuth_Cell_Types_2021', 'BioCarta_2013', 'BioCarta_2015', 'BioCarta_2016', 'BioPlanet_2019', 'BioPlex_2017', 'CCLE_Proteomics_2020', 'CORUM', 'COVID-19_Related_Gene_Sets', 'COVID-19_Related_Gene_Sets_2021', 'Cancer_Cell_Line_Encyclopedia', 'CellMarker_2024', 'CellMarker_Augmented_2021', 'ChEA_2013', 'ChEA_2015', 'ChEA_2016', 'ChEA_2022', 'Chromosome_Location', 'Chromosome_Location_hg19', 'ClinVar_2019', 'DGIdb_Drug_Targets_2024', 'DSigDB', 'Data_Acquisition_Method_Most_Popular_Genes', 'DepMap_CRISPR_GeneDependency_CellLines_2023', 'DepMap_WG_CRISPR_Screens_Broad_CellLines_2019', 'DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019', 'Descartes_Cell_Types_and_Tissue_2021', 'Diabetes_Perturbations_GEO_2022', 'DisGeNET', 'Disease_Perturbations_from_GEO_down', 'Disease_Perturbations_from_GEO_up', 'Disease_Signatures_from_GEO_down_2014', 'Disease_Signatures_from_GEO_up_2014', 'DrugMatrix', 'Drug_Perturbations_from_GEO_2014', 'Drug_Perturbations_from_GEO_down', 'Drug_Perturbations_from_GEO_up', 'ENCODE_Histone_Modifications_2013', 'ENCODE_Histone_Modifications_2015', 'ENCODE_TF_ChIP-seq_2014', 'ENCODE_TF_ChIP-seq_2015', 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X', 'ESCAPE', 'Elsevier_Pathway_Collection', 'Enrichr_Libraries_Most_Popular_Genes', 'Enrichr_Submissions_TF-Gene_Coocurrence', 'Enrichr_Users_Contributed_Lists_2020', 'Epigenomics_Roadmap_HM_ChIP-seq', 'FANTOM6_lncRNA_KD_DEGs', 'GO_Biological_Process_2013', 'GO_Biological_Process_2015', 'GO_Biological_Process_2017', 'GO_Biological_Process_2017b', 'GO_Biological_Process_2018', 'GO_Biological_Process_2021', 'GO_Biological_Process_2023', 'GO_Biological_Process_2025', 'GO_Cellular_Component_2013', 'GO_Cellular_Component_2015', 'GO_Cellular_Component_2017', 'GO_Cellular_Component_2017b', 'GO_Cellular_Component_2018', 'GO_Cellular_Component_2021', 'GO_Cellular_Component_2023', 'GO_Cellular_Component_2025', 'GO_Molecular_Function_2013', 'GO_Molecular_Function_2015', 'GO_Molecular_Function_2017', 'GO_Molecular_Function_2017b', 'GO_Molecular_Function_2018', 'GO_Molecular_Function_2021', 'GO_Molecular_Function_2023', 'GO_Molecular_Function_2025', 'GTEx_Aging_Signatures_2021', 'GTEx_Tissue_Expression_Down', 'GTEx_Tissue_Expression_Up', 'GTEx_Tissues_V8_2023', 'GWAS_Catalog_2019', 'GWAS_Catalog_2023', 'GeDiPNet_2023', 'GeneSigDB', 'Gene_Perturbations_from_GEO_down', 'Gene_Perturbations_from_GEO_up', 'Genes_Associated_with_NIH_Grants', 'Genome_Browser_PWMs', 'GlyGen_Glycosylated_Proteins_2022', 'HDSigDB_Human_2021', 'HDSigDB_Mouse_2021', 'HMDB_Metabolites', 'HMS_LINCS_KinomeScan', 'HomoloGene', 'HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression', 'HuBMAP_ASCTplusB_augmented_2022', 'HumanCyc_2015', 'HumanCyc_2016', 'Human_Gene_Atlas', 'Human_Phenotype_Ontology', 'IDG_Drug_Targets_2022', 'InterPro_Domains_2019', 'Jensen_COMPARTMENTS', 'Jensen_DISEASES', 'Jensen_DISEASES_Curated_2025', 'Jensen_DISEASES_Experimental_2025', 'Jensen_TISSUES', 'KEA_2013', 'KEA_2015', 'KEGG_2013', 'KEGG_2015', 'KEGG_2016', 'KEGG_2019_Human', 'KEGG_2019_Mouse', 'KEGG_2021_Human', 'KOMP2_Mouse_Phenotypes_2022', 'Kinase_Perturbations_from_GEO_down', 'Kinase_Perturbations_from_GEO_up', 'L1000_Kinase_and_GPCR_Perturbations_down', 'L1000_Kinase_and_GPCR_Perturbations_up', 'LINCS_L1000_CRISPR_KO_Consensus_Sigs', 'LINCS_L1000_Chem_Pert_Consensus_Sigs', 'LINCS_L1000_Chem_Pert_down', 'LINCS_L1000_Chem_Pert_up', 'LINCS_L1000_Ligand_Perturbations_down', 'LINCS_L1000_Ligand_Perturbations_up', 'Ligand_Perturbations_from_GEO_down', 'Ligand_Perturbations_from_GEO_up', 'MAGMA_Drugs_and_Diseases', 'MAGNET_2023', 'MCF7_Perturbations_from_GEO_down', 'MCF7_Perturbations_from_GEO_up', 'MGI_Mammalian_Phenotype_2013', 'MGI_Mammalian_Phenotype_2017', 'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotype_Level_4', 'MGI_Mammalian_Phenotype_Level_4_2019', 'MGI_Mammalian_Phenotype_Level_4_2021', 'MGI_Mammalian_Phenotype_Level_4_2024', 'MSigDB_Computational', 'MSigDB_Hallmark_2020', 'MSigDB_Oncogenic_Signatures', 'Metabolomics_Workbench_Metabolites_2022', 'Microbe_Perturbations_from_GEO_down', 'Microbe_Perturbations_from_GEO_up', 'MoTrPAC_2023', 'Mouse_Gene_Atlas', 'NCI-60_Cancer_Cell_Lines', 'NCI-Nature_2016', 'NIBR_DRUGseq_2025_down', 'NIBR_DRUGseq_2025_up', 'NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions', 'NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions', 'NIH_Funded_PIs_2017_Human_AutoRIF', 'NIH_Funded_PIs_2017_Human_GeneRIF', 'NURSA_Human_Endogenous_Complexome', 'OMIM_Disease', 'OMIM_Expanded', 'Old_CMAP_down', 'Old_CMAP_up', 'Orphanet_Augmented_2021', 'PFOCR_Pathways', 'PFOCR_Pathways_2023', 'PPI_Hub_Proteins', 'PanglaoDB_Augmented_2021', 'Panther_2015', 'Panther_2016', 'PerturbAtlas', 'Pfam_Domains_2019', 'Pfam_InterPro_Domains', 'PheWeb_2019', 'PhenGenI_Association_2021', 'Phosphatase_Substrates_from_DEPOD', 'ProteomicsDB_2020', 'Proteomics_Drug_Atlas_2023', 'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO', 'RNAseq_Automatic_GEO_Signatures_Human_Down', 'RNAseq_Automatic_GEO_Signatures_Human_Up', 'RNAseq_Automatic_GEO_Signatures_Mouse_Down', 'RNAseq_Automatic_GEO_Signatures_Mouse_Up', 'Rare_Diseases_AutoRIF_ARCHS4_Predictions', 'Rare_Diseases_AutoRIF_Gene_Lists', 'Rare_Diseases_GeneRIF_ARCHS4_Predictions', 'Rare_Diseases_GeneRIF_Gene_Lists', 'Reactome_2013', 'Reactome_2015', 'Reactome_2016', 'Reactome_2022', 'Reactome_Pathways_2024', 'Rummagene_kinases', 'Rummagene_signatures', 'Rummagene_transcription_factors', 'SILAC_Phosphoproteomics', 'SubCell_BarCode', 'SynGO_2022', 'SynGO_2024', 'SysMyo_Muscle_Gene_Sets', 'TF-LOF_Expression_from_GEO', 'TF_Perturbations_Followed_by_Expression', 'TG_GATES_2020', 'TRANSFAC_and_JASPAR_PWMs', 'TRRUST_Transcription_Factors_2019', 'Table_Mining_of_CRISPR_Studies', 'Tabula_Muris', 'Tabula_Sapiens', 'TargetScan_microRNA', 'TargetScan_microRNA_2017', 'The_Kinase_Library_2023', 'The_Kinase_Library_2024', 'Tissue_Protein_Expression_from_Human_Proteome_Map', 'Tissue_Protein_Expression_from_ProteomicsDB', 'Transcription_Factor_PPIs', 'UK_Biobank_GWAS_v1', 'Virus-Host_PPI_P-HIPSTer_2020', 'VirusMINT', 'Virus_Perturbations_from_GEO_down', 'Virus_Perturbations_from_GEO_up', 'WikiPathway_2021_Human', 'WikiPathway_2023_Human', 'WikiPathways_2013', 'WikiPathways_2015', 'WikiPathways_2016', 'WikiPathways_2019_Human', 'WikiPathways_2019_Mouse', 'WikiPathways_2024_Human', 'WikiPathways_2024_Mouse', 'dbGaP', 'huMAP', 'lncHUB_lncRNA_Co-Expression', 'miRTarBase_2017']


Get the significant DEGs for the Covid patients.

    6





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Gene_set</th>
      <th>Term</th>
      <th>Overlap</th>
      <th>P-value</th>
      <th>Adjusted P-value</th>
      <th>Old P-value</th>
      <th>Old Adjusted P-value</th>
      <th>Odds Ratio</th>
      <th>Combined Score</th>
      <th>Genes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>GO_Biological_Process_2018</td>
      <td>negative regulation of transforming growth fac...</td>
      <td>1/6</td>
      <td>0.001799</td>
      <td>0.028515</td>
      <td>0</td>
      <td>0</td>
      <td>799.560000</td>
      <td>5053.709148</td>
      <td>HSP90AB1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>GO_Biological_Process_2018</td>
      <td>negative regulation of phospholipid metabolic ...</td>
      <td>1/6</td>
      <td>0.001799</td>
      <td>0.028515</td>
      <td>0</td>
      <td>0</td>
      <td>799.560000</td>
      <td>5053.709148</td>
      <td>PIK3IP1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>GO_Biological_Process_2018</td>
      <td>leukocyte aggregation (GO:0070486)</td>
      <td>1/7</td>
      <td>0.002098</td>
      <td>0.028515</td>
      <td>0</td>
      <td>0</td>
      <td>666.266667</td>
      <td>4108.590826</td>
      <td>S100A9</td>
    </tr>
    <tr>
      <th>3</th>
      <td>GO_Biological_Process_2018</td>
      <td>positive regulation of protein localization to...</td>
      <td>1/8</td>
      <td>0.002398</td>
      <td>0.028515</td>
      <td>0</td>
      <td>0</td>
      <td>571.057143</td>
      <td>3445.289983</td>
      <td>HSP90AB1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>GO_Biological_Process_2018</td>
      <td>virion attachment to host cell (GO:0019062)</td>
      <td>1/9</td>
      <td>0.002697</td>
      <td>0.028515</td>
      <td>0</td>
      <td>0</td>
      <td>499.650000</td>
      <td>2955.689739</td>
      <td>HSP90AB1</td>
    </tr>
  </tbody>
</table>
</div>



Some databases of interest:\
`GO_Biological_Process_2017b``KEGG_2019_Human``KEGG_2019_Mouse``WikiPathways_2019_Human``WikiPathways_2019_Mouse`\
You visualize your results using a simple barplot, for example:




    <Axes: title={'center': 'GO_Biological_Process_2018'}, xlabel='$- \\log_{10}$ (Adjusted P-value)'>




    
![png](scanpy_05_dge_files/scanpy_05_dge_50_1.png)
    


## Gene Set Enrichment Analysis (GSEA)

Besides the enrichment using hypergeometric test, we can also perform
gene set enrichment analysis (GSEA), which scores ranked genes list
(usually based on fold changes) and computes permutation test to check
if a particular gene set is more present in the Up-regulated genes,
among the DOWN_regulated genes or not differentially regulated.

We need a table with all DEGs and their log foldchanges. However, many
lowly expressed genes will have high foldchanges and just contribue
noise, so also filter for expression in enough cells.




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>names</th>
      <th>logfoldchanges</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>151</th>
      <td>PF4</td>
      <td>28.017332</td>
    </tr>
    <tr>
      <th>210</th>
      <td>MTRNR2L1</td>
      <td>27.951775</td>
    </tr>
    <tr>
      <th>316</th>
      <td>CXCL8</td>
      <td>27.453955</td>
    </tr>
    <tr>
      <th>380</th>
      <td>G0S2</td>
      <td>27.165596</td>
    </tr>
    <tr>
      <th>682</th>
      <td>CCL3</td>
      <td>26.821484</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>17544</th>
      <td>CSTA</td>
      <td>-25.910336</td>
    </tr>
    <tr>
      <th>17546</th>
      <td>MIAT</td>
      <td>-25.982347</td>
    </tr>
    <tr>
      <th>17538</th>
      <td>TTC34</td>
      <td>-26.053516</td>
    </tr>
    <tr>
      <th>17901</th>
      <td>CEP41</td>
      <td>-26.089569</td>
    </tr>
    <tr>
      <th>17902</th>
      <td>REM2</td>
      <td>-26.238188</td>
    </tr>
  </tbody>
</table>
<p>3850 rows × 2 columns</p>
</div>



Once our list of genes are sorted, we can proceed with the enrichment
itself. We can use the package to get gene set from the Molecular
Signature Database (MSigDB) and select KEGG pathways as an example.

    ['ARCHS4_Cell-lines', 'ARCHS4_IDG_Coexp', 'ARCHS4_Kinases_Coexp', 'ARCHS4_TFs_Coexp', 'ARCHS4_Tissues', 'Achilles_fitness_decrease', 'Achilles_fitness_increase', 'Aging_Perturbations_from_GEO_down', 'Aging_Perturbations_from_GEO_up', 'Allen_Brain_Atlas_10x_scRNA_2021', 'Allen_Brain_Atlas_down', 'Allen_Brain_Atlas_up', 'Azimuth_2023', 'Azimuth_Cell_Types_2021', 'BioCarta_2013', 'BioCarta_2015', 'BioCarta_2016', 'BioPlanet_2019', 'BioPlex_2017', 'CCLE_Proteomics_2020', 'CORUM', 'COVID-19_Related_Gene_Sets', 'COVID-19_Related_Gene_Sets_2021', 'Cancer_Cell_Line_Encyclopedia', 'CellMarker_2024', 'CellMarker_Augmented_2021', 'ChEA_2013', 'ChEA_2015', 'ChEA_2016', 'ChEA_2022', 'Chromosome_Location', 'Chromosome_Location_hg19', 'ClinVar_2019', 'DGIdb_Drug_Targets_2024', 'DSigDB', 'Data_Acquisition_Method_Most_Popular_Genes', 'DepMap_CRISPR_GeneDependency_CellLines_2023', 'DepMap_WG_CRISPR_Screens_Broad_CellLines_2019', 'DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019', 'Descartes_Cell_Types_and_Tissue_2021', 'Diabetes_Perturbations_GEO_2022', 'DisGeNET', 'Disease_Perturbations_from_GEO_down', 'Disease_Perturbations_from_GEO_up', 'Disease_Signatures_from_GEO_down_2014', 'Disease_Signatures_from_GEO_up_2014', 'DrugMatrix', 'Drug_Perturbations_from_GEO_2014', 'Drug_Perturbations_from_GEO_down', 'Drug_Perturbations_from_GEO_up', 'ENCODE_Histone_Modifications_2013', 'ENCODE_Histone_Modifications_2015', 'ENCODE_TF_ChIP-seq_2014', 'ENCODE_TF_ChIP-seq_2015', 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X', 'ESCAPE', 'Elsevier_Pathway_Collection', 'Enrichr_Libraries_Most_Popular_Genes', 'Enrichr_Submissions_TF-Gene_Coocurrence', 'Enrichr_Users_Contributed_Lists_2020', 'Epigenomics_Roadmap_HM_ChIP-seq', 'FANTOM6_lncRNA_KD_DEGs', 'GO_Biological_Process_2013', 'GO_Biological_Process_2015', 'GO_Biological_Process_2017', 'GO_Biological_Process_2017b', 'GO_Biological_Process_2018', 'GO_Biological_Process_2021', 'GO_Biological_Process_2023', 'GO_Biological_Process_2025', 'GO_Cellular_Component_2013', 'GO_Cellular_Component_2015', 'GO_Cellular_Component_2017', 'GO_Cellular_Component_2017b', 'GO_Cellular_Component_2018', 'GO_Cellular_Component_2021', 'GO_Cellular_Component_2023', 'GO_Cellular_Component_2025', 'GO_Molecular_Function_2013', 'GO_Molecular_Function_2015', 'GO_Molecular_Function_2017', 'GO_Molecular_Function_2017b', 'GO_Molecular_Function_2018', 'GO_Molecular_Function_2021', 'GO_Molecular_Function_2023', 'GO_Molecular_Function_2025', 'GTEx_Aging_Signatures_2021', 'GTEx_Tissue_Expression_Down', 'GTEx_Tissue_Expression_Up', 'GTEx_Tissues_V8_2023', 'GWAS_Catalog_2019', 'GWAS_Catalog_2023', 'GeDiPNet_2023', 'GeneSigDB', 'Gene_Perturbations_from_GEO_down', 'Gene_Perturbations_from_GEO_up', 'Genes_Associated_with_NIH_Grants', 'Genome_Browser_PWMs', 'GlyGen_Glycosylated_Proteins_2022', 'HDSigDB_Human_2021', 'HDSigDB_Mouse_2021', 'HMDB_Metabolites', 'HMS_LINCS_KinomeScan', 'HomoloGene', 'HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression', 'HuBMAP_ASCTplusB_augmented_2022', 'HumanCyc_2015', 'HumanCyc_2016', 'Human_Gene_Atlas', 'Human_Phenotype_Ontology', 'IDG_Drug_Targets_2022', 'InterPro_Domains_2019', 'Jensen_COMPARTMENTS', 'Jensen_DISEASES', 'Jensen_DISEASES_Curated_2025', 'Jensen_DISEASES_Experimental_2025', 'Jensen_TISSUES', 'KEA_2013', 'KEA_2015', 'KEGG_2013', 'KEGG_2015', 'KEGG_2016', 'KEGG_2019_Human', 'KEGG_2019_Mouse', 'KEGG_2021_Human', 'KOMP2_Mouse_Phenotypes_2022', 'Kinase_Perturbations_from_GEO_down', 'Kinase_Perturbations_from_GEO_up', 'L1000_Kinase_and_GPCR_Perturbations_down', 'L1000_Kinase_and_GPCR_Perturbations_up', 'LINCS_L1000_CRISPR_KO_Consensus_Sigs', 'LINCS_L1000_Chem_Pert_Consensus_Sigs', 'LINCS_L1000_Chem_Pert_down', 'LINCS_L1000_Chem_Pert_up', 'LINCS_L1000_Ligand_Perturbations_down', 'LINCS_L1000_Ligand_Perturbations_up', 'Ligand_Perturbations_from_GEO_down', 'Ligand_Perturbations_from_GEO_up', 'MAGMA_Drugs_and_Diseases', 'MAGNET_2023', 'MCF7_Perturbations_from_GEO_down', 'MCF7_Perturbations_from_GEO_up', 'MGI_Mammalian_Phenotype_2013', 'MGI_Mammalian_Phenotype_2017', 'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotype_Level_4', 'MGI_Mammalian_Phenotype_Level_4_2019', 'MGI_Mammalian_Phenotype_Level_4_2021', 'MGI_Mammalian_Phenotype_Level_4_2024', 'MSigDB_Computational', 'MSigDB_Hallmark_2020', 'MSigDB_Oncogenic_Signatures', 'Metabolomics_Workbench_Metabolites_2022', 'Microbe_Perturbations_from_GEO_down', 'Microbe_Perturbations_from_GEO_up', 'MoTrPAC_2023', 'Mouse_Gene_Atlas', 'NCI-60_Cancer_Cell_Lines', 'NCI-Nature_2016', 'NIBR_DRUGseq_2025_down', 'NIBR_DRUGseq_2025_up', 'NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions', 'NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions', 'NIH_Funded_PIs_2017_Human_AutoRIF', 'NIH_Funded_PIs_2017_Human_GeneRIF', 'NURSA_Human_Endogenous_Complexome', 'OMIM_Disease', 'OMIM_Expanded', 'Old_CMAP_down', 'Old_CMAP_up', 'Orphanet_Augmented_2021', 'PFOCR_Pathways', 'PFOCR_Pathways_2023', 'PPI_Hub_Proteins', 'PanglaoDB_Augmented_2021', 'Panther_2015', 'Panther_2016', 'PerturbAtlas', 'Pfam_Domains_2019', 'Pfam_InterPro_Domains', 'PheWeb_2019', 'PhenGenI_Association_2021', 'Phosphatase_Substrates_from_DEPOD', 'ProteomicsDB_2020', 'Proteomics_Drug_Atlas_2023', 'RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO', 'RNAseq_Automatic_GEO_Signatures_Human_Down', 'RNAseq_Automatic_GEO_Signatures_Human_Up', 'RNAseq_Automatic_GEO_Signatures_Mouse_Down', 'RNAseq_Automatic_GEO_Signatures_Mouse_Up', 'Rare_Diseases_AutoRIF_ARCHS4_Predictions', 'Rare_Diseases_AutoRIF_Gene_Lists', 'Rare_Diseases_GeneRIF_ARCHS4_Predictions', 'Rare_Diseases_GeneRIF_Gene_Lists', 'Reactome_2013', 'Reactome_2015', 'Reactome_2016', 'Reactome_2022', 'Reactome_Pathways_2024', 'Rummagene_kinases', 'Rummagene_signatures', 'Rummagene_transcription_factors', 'SILAC_Phosphoproteomics', 'SubCell_BarCode', 'SynGO_2022', 'SynGO_2024', 'SysMyo_Muscle_Gene_Sets', 'TF-LOF_Expression_from_GEO', 'TF_Perturbations_Followed_by_Expression', 'TG_GATES_2020', 'TRANSFAC_and_JASPAR_PWMs', 'TRRUST_Transcription_Factors_2019', 'Table_Mining_of_CRISPR_Studies', 'Tabula_Muris', 'Tabula_Sapiens', 'TargetScan_microRNA', 'TargetScan_microRNA_2017', 'The_Kinase_Library_2023', 'The_Kinase_Library_2024', 'Tissue_Protein_Expression_from_Human_Proteome_Map', 'Tissue_Protein_Expression_from_ProteomicsDB', 'Transcription_Factor_PPIs', 'UK_Biobank_GWAS_v1', 'Virus-Host_PPI_P-HIPSTer_2020', 'VirusMINT', 'Virus_Perturbations_from_GEO_down', 'Virus_Perturbations_from_GEO_up', 'WikiPathway_2021_Human', 'WikiPathway_2023_Human', 'WikiPathways_2013', 'WikiPathways_2015', 'WikiPathways_2016', 'WikiPathways_2019_Human', 'WikiPathways_2019_Mouse', 'WikiPathways_2024_Human', 'WikiPathways_2024_Mouse', 'dbGaP', 'huMAP', 'lncHUB_lncRNA_Co-Expression', 'miRTarBase_2017']


Next, we will run GSEA. This will result in a table containing
information for several pathways. We can then sort and filter those
pathways to visualize only the top ones. You can select/filter them by
either `p-value` or normalized enrichment score (`NES`).

    2025-05-23 17:26:55,432 [WARNING] Duplicated values found in preranked stats: 13.74% of genes
    The order of those genes will be arbitrary, which may produce unexpected results.





    0                       Sphingolipid signaling pathway
    1                    T cell receptor signaling pathway
    2             C-type lectin receptor signaling pathway
    3    PD-L1 expression and PD-1 checkpoint pathway i...
    4                                  Cellular senescence
    5                                  Cholinergic synapse
    6                     Th1 and Th2 cell differentiation
    7                              IL-17 signaling pathway
    8                                        Leishmaniasis
    9       Growth hormone synthesis, secretion and action
    Name: Term, dtype: object






    [<Axes: xlabel='Gene Rank', ylabel='Ranked metric'>,
     <Axes: >,
     <Axes: >,
     <Axes: ylabel='Enrichment Score'>]




    
![png](scanpy_05_dge_files/scanpy_05_dge_57_1.png)
    


<div>

> **Discuss**
>
> Which KEGG pathways are upregulated in this cluster? Which KEGG
> pathways are dowregulated in this cluster? Change the pathway source
> to another gene set (e.g. CP:WIKIPATHWAYS or CP:REACTOME or
> CP:BIOCARTA or GO:BP) and check the if you get similar results?

</div>

Finally, let's save the integrated data for further analysis.

## Session info

```{=html}
<details>
```
```{=html}
<summary>
```
Click here
```{=html}
</summary>
```




<table class=table>
            <thead style="position: sticky; top: 0; background-color: var(--jp-layout-color0, var(--vscode-editor-background, white));">
        <tr><th>Package</th><th>Version</th></tr>
    </thead>
    <tbody>
        <tr><td><strong>pandas</strong></td><td>2.2.3</td></tr>
        <tr><td><strong>numpy</strong></td><td>2.2.5</td></tr>
        <tr><td><strong>scanpy</strong></td><td>1.11.1</td></tr>
        <tr><td><strong>gseapy</strong></td><td>1.1.8</td></tr>
        <tr><td><strong>matplotlib</strong></td><td>3.10.3</td></tr>
        <tr><td><strong>anndata</strong></td><td>0.11.4</td></tr>
        <tr><td><strong>matplotlib-venn</strong></td><td>1.1.2</td></tr>
        <tr><td><strong>seaborn</strong></td><td>0.13.2</td></tr>
    </tbody>
    <thead style="position: sticky; top: 0; background-color: var(--jp-layout-color0, var(--vscode-editor-background, white));">
        <tr><th>Component</th><th>Info</th></tr>
    </thead>
    <tbody>
        <tr><td>Python</td><td>3.12.9 | packaged by Anaconda, Inc. | (main, Feb  6 2025, 12:55:12) [Clang 14.0.6 ]</td></tr>
        <tr><td>OS</td><td>macOS-15.5-arm64-arm-64bit</td></tr>
        <tr><td>CPU</td><td>12 logical CPU cores, arm</td></tr>
        <tr><td>GPU</td><td>No GPU found</td></tr>
        <tr><td>Updated</td><td>2025-05-23 09:27</td></tr>
    </tbody>
        </table>

        <details>
        <summary>Dependencies</summary>
                <div style="max-height: min(500px, 80vh); overflow-y: auto;">
    <table class=table>
            <thead style="position: sticky; top: 0; background-color: var(--jp-layout-color0, var(--vscode-editor-background, white));">
    <tr><th>Dependency</th><th>Version</th></tr>
</thead>
<tbody>
    <tr><td>parso</td><td>0.8.4</td></tr>
    <tr><td>attrs</td><td>24.3.0</td></tr>
    <tr><td>requests</td><td>2.32.3</td></tr>
    <tr><td>pytz</td><td>2025.2</td></tr>
    <tr><td>urllib3</td><td>2.3.0</td></tr>
    <tr><td>scipy</td><td>1.15.3</td></tr>
    <tr><td>h5py</td><td>3.13.0</td></tr>
    <tr><td>threadpoolctl</td><td>3.6.0</td></tr>
    <tr><td>scikit-learn</td><td>1.5.2</td></tr>
    <tr><td>pyparsing</td><td>3.2.3</td></tr>
    <tr><td>future</td><td>1.0.0</td></tr>
    <tr><td>ipython</td><td>9.1.0</td></tr>
    <tr><td>defusedxml</td><td>0.7.1</td></tr>
    <tr><td>jupyter_client</td><td>8.6.3</td></tr>
    <tr><td>pure-eval</td><td>0.2.2</td></tr>
    <tr><td>matplotlib-inline</td><td>0.1.6</td></tr>
    <tr><td>cattrs</td><td>24.1.3</td></tr>
    <tr><td>typing_extensions</td><td>4.12.2</td></tr>
    <tr><td>jedi</td><td>0.19.2</td></tr>
    <tr><td>idna</td><td>3.7</td></tr>
    <tr><td>traitlets</td><td>5.14.3</td></tr>
    <tr><td>url-normalize</td><td>2.2.1</td></tr>
    <tr><td>six</td><td>1.17.0</td></tr>
    <tr><td>pillow</td><td>11.2.1</td></tr>
    <tr><td>Pygments</td><td>2.19.1</td></tr>
    <tr><td>debugpy</td><td>1.8.11</td></tr>
    <tr><td>patsy</td><td>1.0.1</td></tr>
    <tr><td>pycparser</td><td>2.21</td></tr>
    <tr><td>packaging</td><td>24.2</td></tr>
    <tr><td>Cython</td><td>3.1.0</td></tr>
    <tr><td>Brotli</td><td>1.0.9</td></tr>
    <tr><td>statsmodels</td><td>0.14.4</td></tr>
    <tr><td>certifi</td><td>2025.4.26 (2025.04.26)</td></tr>
    <tr><td>cycler</td><td>0.12.1</td></tr>
    <tr><td>natsort</td><td>8.4.0</td></tr>
    <tr><td>igraph</td><td>0.11.8</td></tr>
    <tr><td>python-dateutil</td><td>2.9.0.post0</td></tr>
    <tr><td>ipykernel</td><td>6.29.5</td></tr>
    <tr><td>executing</td><td>0.8.3</td></tr>
    <tr><td>numba</td><td>0.61.2</td></tr>
    <tr><td>decorator</td><td>5.1.1</td></tr>
    <tr><td>jupyter_core</td><td>5.7.2</td></tr>
    <tr><td>leidenalg</td><td>0.10.2</td></tr>
    <tr><td>cffi</td><td>1.17.1</td></tr>
    <tr><td>joblib</td><td>1.5.0</td></tr>
    <tr><td>charset-normalizer</td><td>3.3.2</td></tr>
    <tr><td>platformdirs</td><td>4.3.7</td></tr>
    <tr><td>texttable</td><td>1.7.0</td></tr>
    <tr><td>llvmlite</td><td>0.44.0</td></tr>
    <tr><td>psutil</td><td>5.9.0</td></tr>
    <tr><td>session-info2</td><td>0.1.2</td></tr>
    <tr><td>prompt-toolkit</td><td>3.0.43</td></tr>
    <tr><td>appnope</td><td>0.1.3</td></tr>
    <tr><td>kiwisolver</td><td>1.4.8</td></tr>
    <tr><td>pybiomart</td><td>0.2.0</td></tr>
    <tr><td>PySocks</td><td>1.7.1</td></tr>
    <tr><td>PyYAML</td><td>6.0.2</td></tr>
    <tr><td>setuptools</td><td>78.1.1</td></tr>
    <tr><td>wcwidth</td><td>0.2.5</td></tr>
    <tr><td>legacy-api-wrap</td><td>1.4.1</td></tr>
    <tr><td>tornado</td><td>6.4.2</td></tr>
    <tr><td>pyzmq</td><td>26.2.0</td></tr>
    <tr><td>asttokens</td><td>3.0.0</td></tr>
    <tr><td>stack-data</td><td>0.2.0</td></tr>
    <tr><td>requests-cache</td><td>1.2.1</td></tr>
    <tr><td>comm</td><td>0.2.1</td></tr>
</tbody>
    </table>
</div>
    </details>
        <details>
            <summary>Copyable Markdown</summary>
            <pre>| Package         | Version |
| --------------- | ------- |
| pandas          | 2.2.3   |
| numpy           | 2.2.5   |
| scanpy          | 1.11.1  |
| gseapy          | 1.1.8   |
| matplotlib      | 3.10.3  |
| anndata         | 0.11.4  |
| matplotlib-venn | 1.1.2   |
| seaborn         | 0.13.2  |

| Dependency         | Version                |
| ------------------ | ---------------------- |
| parso              | 0.8.4                  |
| attrs              | 24.3.0                 |
| requests           | 2.32.3                 |
| pytz               | 2025.2                 |
| urllib3            | 2.3.0                  |
| scipy              | 1.15.3                 |
| h5py               | 3.13.0                 |
| threadpoolctl      | 3.6.0                  |
| scikit-learn       | 1.5.2                  |
| pyparsing          | 3.2.3                  |
| future             | 1.0.0                  |
| ipython            | 9.1.0                  |
| defusedxml         | 0.7.1                  |
| jupyter_client     | 8.6.3                  |
| pure-eval          | 0.2.2                  |
| matplotlib-inline  | 0.1.6                  |
| cattrs             | 24.1.3                 |
| typing_extensions  | 4.12.2                 |
| jedi               | 0.19.2                 |
| idna               | 3.7                    |
| traitlets          | 5.14.3                 |
| url-normalize      | 2.2.1                  |
| six                | 1.17.0                 |
| pillow             | 11.2.1                 |
| Pygments           | 2.19.1                 |
| debugpy            | 1.8.11                 |
| patsy              | 1.0.1                  |
| pycparser          | 2.21                   |
| packaging          | 24.2                   |
| Cython             | 3.1.0                  |
| Brotli             | 1.0.9                  |
| statsmodels        | 0.14.4                 |
| certifi            | 2025.4.26 (2025.04.26) |
| cycler             | 0.12.1                 |
| natsort            | 8.4.0                  |
| igraph             | 0.11.8                 |
| python-dateutil    | 2.9.0.post0            |
| ipykernel          | 6.29.5                 |
| executing          | 0.8.3                  |
| numba              | 0.61.2                 |
| decorator          | 5.1.1                  |
| jupyter_core       | 5.7.2                  |
| leidenalg          | 0.10.2                 |
| cffi               | 1.17.1                 |
| joblib             | 1.5.0                  |
| charset-normalizer | 3.3.2                  |
| platformdirs       | 4.3.7                  |
| texttable          | 1.7.0                  |
| llvmlite           | 0.44.0                 |
| psutil             | 5.9.0                  |
| session-info2      | 0.1.2                  |
| prompt-toolkit     | 3.0.43                 |
| appnope            | 0.1.3                  |
| kiwisolver         | 1.4.8                  |
| pybiomart          | 0.2.0                  |
| PySocks            | 1.7.1                  |
| PyYAML             | 6.0.2                  |
| setuptools         | 78.1.1                 |
| wcwidth            | 0.2.5                  |
| legacy-api-wrap    | 1.4.1                  |
| tornado            | 6.4.2                  |
| pyzmq              | 26.2.0                 |
| asttokens          | 3.0.0                  |
| stack-data         | 0.2.0                  |
| requests-cache     | 1.2.1                  |
| comm               | 0.2.1                  |

| Component | Info                                                                                |
| --------- | ----------------------------------------------------------------------------------- |
| Python    | 3.12.9 | packaged by Anaconda, Inc. | (main, Feb  6 2025, 12:55:12) [Clang 14.0.6 ] |
| OS        | macOS-15.5-arm64-arm-64bit                                                          |
| CPU       | 12 logical CPU cores, arm                                                           |
| GPU       | No GPU found                                                                        |
| Updated   | 2025-05-23 09:27                                                                    |</pre>
        </details>



```{=html}
</details>
```

## archive

    [NbConvertApp] Converting notebook scanpy_05_dge.ipynb to markdown
    [NbConvertApp] Support files will be in scanpy_05_dge_files/
    [NbConvertApp] Making directory archive/scanpy_05_dge_files
    [NbConvertApp] Writing 15106 bytes to archive/scanpy_05_dge.md

