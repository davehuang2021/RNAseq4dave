# Pseudobulk Enrichment Analysis

When cell identity clusters are well-defined, it can be advantageous to perform
analyses at the pseudobulk level rather than at the single-cell level.
Pseudobulking involves aggregating counts across cells of the same type within
each sample, effectively creating sample-level gene expression profiles per cell
type. This approach helps mitigate the effects of technical noise and dropouts
common in single-cell data, enabling the detection of lowly expressed genes that
might otherwise be missed.

Moreover, conducting differential expression analysis (DEA) at the pseudobulk
level, treating each biological sample as the unit of observation, is
statistically more robust.
Unlike single-cell DEA, which assumes cells are independent (an assumption
that is violated when cells originate from the same individual), sample-level
pseudobulk analysis avoids inflation of p-values by reducing the number of
observations and by correctly modeling biological replication {cite:p}`psbulk`.

The resulting gene-level statistics from pseudobulk DEA can then be used as
input for downstream enrichment analyses.

In this notebook, we demonstrate how to use `decoupler` to infer transcription factor
(TF) and pathway enrichment scores from a multi-sample scRNA-seq human dataset.

The dataset includes 5k peripheral blood mononuclear cells (PBMCs) from healthy and
COVID-19 infected patients {cite:p}`covid5k`. It publicly available at the Single Cell
Expression Atlas ([E-MTAB-9221](https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-9221/)).


- [ref](https://decoupler.readthedocs.io/en/latest/notebooks/scell/rna_psbk.html)
> 数据是从以上链接下载。如果发现不能访问，可以尝试去以下链接寻找

- [backup ref link (decoupler)](https://decoupler.readthedocs.io/)


## Loading Packages

## Loading The Dataset




    AnnData object with n_obs × n_vars = 4903 × 14120
        obs: 'individual', 'sex', 'disease', 'celltype'
        obsm: 'X_umap'



The obtained {class}`anndata.AnnData` consist of raw integer transcript counts for ~5k cells
with measurements for ~15k genes.

The cell metadata stored in {attr}`anndata.AnnData.obs` can be inspected.




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
      <th>individual</th>
      <th>sex</th>
      <th>disease</th>
      <th>celltype</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>SAMEA6979322-AGGGTTTAGGGCGAAG</th>
      <td>SARS-CoV2 pos Severe #3</td>
      <td>not available</td>
      <td>COVID-19</td>
      <td>B cell</td>
    </tr>
    <tr>
      <th>SAMEA6979322-CCGTTCAGTTGCTCGG</th>
      <td>SARS-CoV2 pos Severe #3</td>
      <td>not available</td>
      <td>COVID-19</td>
      <td>B cell</td>
    </tr>
    <tr>
      <th>SAMEA6979322-AACGGGAGTCCTGAAT</th>
      <td>SARS-CoV2 pos Severe #3</td>
      <td>not available</td>
      <td>COVID-19</td>
      <td>T cell</td>
    </tr>
    <tr>
      <th>SAMEA6979322-AGGAATATCACGGTCG</th>
      <td>SARS-CoV2 pos Severe #3</td>
      <td>not available</td>
      <td>COVID-19</td>
      <td>T cell</td>
    </tr>
    <tr>
      <th>SAMEA6979322-AGACAAACAACAGCTT</th>
      <td>SARS-CoV2 pos Severe #3</td>
      <td>not available</td>
      <td>COVID-19</td>
      <td>T cell</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>SAMEA6979315-CACCGTTGTTACTCAG</th>
      <td>Control #3</td>
      <td>female</td>
      <td>normal</td>
      <td>T cell</td>
    </tr>
    <tr>
      <th>SAMEA6979315-CATCCACGTCTGTCAA</th>
      <td>Control #3</td>
      <td>female</td>
      <td>normal</td>
      <td>T cell</td>
    </tr>
    <tr>
      <th>SAMEA6979315-CATGCAAAGAAATTGC</th>
      <td>Control #3</td>
      <td>female</td>
      <td>normal</td>
      <td>T cell</td>
    </tr>
    <tr>
      <th>SAMEA6979315-ACCTGTCTCGAAGCCC</th>
      <td>Control #3</td>
      <td>female</td>
      <td>normal</td>
      <td>T cell</td>
    </tr>
    <tr>
      <th>SAMEA6979315-GTCGTAAGTACAACGG</th>
      <td>Control #3</td>
      <td>female</td>
      <td>normal</td>
      <td>T cell</td>
    </tr>
  </tbody>
</table>
<p>4903 rows × 4 columns</p>
</div>



And visualized.


    
![png](rna_psbk_files/rna_psbk_9_0.png)
    


## Pseudobulking

The pseudo-bulk approach involves the following steps:
1. Subsetting the cell type of interest
2. Extracting their raw integer counts
3. Summing their counts per gene into a single profile if they pass quality control

Then, DEA can be performed if there are at least two biological replicates per condition (more replicates are recommended).

Pseudobulking can easily be performed using the function {func}`decoupler.pp.pseudobulk`.
In this example, the counts are just summed, though other modes such as the mean or 
any custom aggregation function are available.
For more information, refer to the `mode` argument.

A profile has been generated for each sample and cell type.
Quality control metrics can be then visualized.


    
![png](rna_psbk_files/rna_psbk_13_0.png)
    


Low-quality samples can be filtered based on two criteria: the number
of cells (`adata.obs.psbulk_cells`) and the total count sum (`adata.obs.psbulk_counts`).

In this dataset, the plots show that some generated profiles contain
fewer than 10 cells and 1000 counts.

These can be excluded using the function {func}`decoupler.pp.filter_samples`.

Thresholds are always dataset-dependent and arbitrary, but a common guideline
is to retain samples with at least 10 cells and 1000 total counts.

After removing low-quality samples, the remaining profiles count can be visualized.


    
![png](rna_psbk_files/rna_psbk_17_0.png)
    


### Variability Exploration
With pseudobulk profiles generated for each cell type and sample,
variability across them can now be explored.

This involves some basic preprocessing followed by principal component analysis (PCA).

After computing the PCs, associations or correlations between each inferred PC and the variables
in the metadata can be tested, depending on whether the variables are categorical or continuous.

This type of analysis is applicable to any dimensionality reduction method, such as factors
derived from non-negative matrix factorization.

The importance of each principal component (based on its explained variance ratio) and its associations
with metadata variables can then be visualized.


    
![png](rna_psbk_files/rna_psbk_23_0.png)
    





    
![png](rna_psbk_files/rna_psbk_23_1.png)
    



In this dataset, PC1 appears to explain the largest proportion of variance and is associated
with the metadata variables individual, disease, and cell type.

Metadata variables associated with PCs that capture a substantial amount of variance are
important and should be accounted for as relevant covariates in downstream differential
expression analysis when possible.

The principal components can also be directly visualized, colored by these metadata variables.


    
![png](rna_psbk_files/rna_psbk_25_0.png)
    


A clear separation can be observed between control and infected patients.

### Feature selection

In addition to filtering low-quality samples, lowly or noisily expressed genes can also be filtered
prior to downstream analyses such as differential expression analysis. This step should be performed
at the cell type level, as different cell types may express distinct sets of genes.

In this vignette, we will focus on investigating the effects of COVID-19 on T cells but this should be done for every cell type of interest.

The first step
is to select the relevant samples.




    AnnData object with n_obs × n_vars = 7 × 14120
        obs: 'individual', 'celltype', 'sex', 'disease', 'psbulk_cells', 'psbulk_counts'
        var: 'mean', 'std'
        uns: 'log1p', 'pca', 'rank_obsm', 'individual_colors', 'disease_colors', 'sex_colors'
        obsm: 'X_pca'
        varm: 'PCs'
        layers: 'psbulk_props', 'counts', 'X'



Two strategies are used to filter genes:

1) `decoupler.pp.filter_by_expr`: Retains genes with a minimum total number of reads across all samples (`min_total_count`)
   and a minimum number of counts in a given number of samples (`min_count`). This approach was introduced in
   [edgeR](https://rdrr.io/bioc/edgeR/man/filterByExpr.html){cite:p}`edger`.
2) `decoupler.pp.filter_by_prop`: Retains genes that are expressed in at least a specified proportion of cells (`min_prop`)
   across a minimum number of samples (`min_smpls`).

The number of retained genes can be visualized, and the filtering parameters can be adjusted interactively.


    
![png](rna_psbk_files/rna_psbk_30_0.png)
    



    
![png](rna_psbk_files/rna_psbk_30_1.png)
    


The top plot displays gene frequencies based on the `filter_by_expr` metrics,
while the bottom plot corresponds to `filter_by_prop`.

Dashed lines indicate the current threshold values.
In the top plot, only genes in the upper-right quadrant are retained, in the
bottom plot, only those to the right of the vertical line are kept.

Although filtering thresholds are arbitrary, a common heuristic is to look
for bimodal distributions and set thresholds that separate low-quality
genes from the rest.

In this example, the default parameters retain a
substantial number of genes while removing potentially noisy ones.

Once the threshold parameters are set, the actual gene filtering can be performed by simply changing `pl` to `pp`.




    AnnData object with n_obs × n_vars = 7 × 1481
        obs: 'individual', 'celltype', 'sex', 'disease', 'psbulk_cells', 'psbulk_counts'
        var: 'mean', 'std'
        uns: 'log1p', 'pca', 'rank_obsm', 'individual_colors', 'disease_colors', 'sex_colors'
        obsm: 'X_pca'
        varm: 'PCs'
        layers: 'psbulk_props', 'counts', 'X'



### Differential Expression Analysis

Differential expression analysis (DEA) can be used to identify genes that differ most between disease and control conditions.
Given the observed association of gene expression with both sex and disease status, the experimental design will include
these covariates to account for their effects.

We will use the Python implementation of the `DESeq2` framework {cite:p}`pydeseq`, though other tools like `limma`
{cite:p}`limma` or `edgeR` {cite:p}`edger` could also be used.
For a deeper understanding of how `pyDESeq2` works, refer to its
[official documentation](https://pydeseq2.readthedocs.io/en/latest/).

Note that even more complex experimental designs can be accommodated by adding additional factors to the `design_factors` argument.

    Using None as control genes, passed at DeseqDataSet initialization


    /var/folders/dq/m9jbf2_x2_v6cnd1c0nt1rc40000gn/T/ipykernel_20803/362220359.py:7: DeprecationWarning: design_factors is deprecated and will soon be removed.Please consider providing a formulaic formula using the design argumentinstead.
      dds = DeseqDataSet(
    Fitting size factors...
    ... done in 0.00 seconds.
    
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/sklearn/linear_model/_base.py:290: RuntimeWarning: divide by zero encountered in matmul
      return X @ coef_.T + self.intercept_
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/sklearn/linear_model/_base.py:290: RuntimeWarning: overflow encountered in matmul
      return X @ coef_.T + self.intercept_
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/sklearn/linear_model/_base.py:290: RuntimeWarning: invalid value encountered in matmul
      return X @ coef_.T + self.intercept_
    Fitting dispersions...
    ... done in 0.08 seconds.
    
    Fitting dispersion trend curve...
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/pydeseq2/default_inference.py:209: RuntimeWarning: divide by zero encountered in matmul
      mu = covariates_fit @ coeffs
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/pydeseq2/default_inference.py:209: RuntimeWarning: overflow encountered in matmul
      mu = covariates_fit @ coeffs
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/pydeseq2/default_inference.py:209: RuntimeWarning: invalid value encountered in matmul
      mu = covariates_fit @ coeffs
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/pydeseq2/default_inference.py:213: RuntimeWarning: divide by zero encountered in matmul
      mu = covariates_fit @ coeffs
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/pydeseq2/default_inference.py:213: RuntimeWarning: overflow encountered in matmul
      mu = covariates_fit @ coeffs
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/pydeseq2/default_inference.py:213: RuntimeWarning: invalid value encountered in matmul
      mu = covariates_fit @ coeffs
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/pydeseq2/default_inference.py:230: RuntimeWarning: divide by zero encountered in matmul
      return coeffs, covariates_fit @ coeffs, res.success
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/pydeseq2/default_inference.py:230: RuntimeWarning: overflow encountered in matmul
      return coeffs, covariates_fit @ coeffs, res.success
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/pydeseq2/default_inference.py:230: RuntimeWarning: invalid value encountered in matmul
      return coeffs, covariates_fit @ coeffs, res.success
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/pydeseq2/dds.py:820: UserWarning: The dispersion trend curve fitting did not converge. Switching to a mean-based dispersion trend.
      self._fit_parametric_dispersion_trend(vst)
    ... done in 0.01 seconds.
    
    /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/pydeseq2/dds.py:548: UserWarning: As the residual degrees of freedom is less than 3, the distribution of log dispersions is especially asymmetric and likely to be poorly estimated by the MAD.
      self.fit_dispersion_prior()
    Fitting MAP dispersions...
    ... done in 0.07 seconds.
    
    Fitting LFCs...


    Log2 fold change & Wald test p-value: disease COVID-19 vs normal
                baseMean  log2FoldChange     lfcSE      stat    pvalue      padj
    TRAC       81.875820       -0.930538  0.542822 -1.714261  0.086481  0.993036
    UCP2       68.862556       -0.548205  0.422626 -1.297142  0.194582  0.993036
    CMBL       43.585692        0.572860  0.460769  1.243269  0.213769  0.993036
    PJA2       46.747220        0.125016  0.468412  0.266893  0.789552  0.993036
    RNF125     59.000601        0.960666  0.458653  2.094537  0.036212  0.993036
    ...              ...             ...       ...       ...       ...       ...
    LYRM7      70.932088        0.057458  0.478781  0.120008  0.904476  0.996294
    PNN       116.380620        0.213590  0.417045  0.512152  0.608545  0.993036
    CDC42EP3   52.644710        0.256493  0.589972  0.434755  0.663740  0.993036
    CHCHD2    167.927239       -0.497018  0.449454 -1.105825  0.268802  0.993036
    RPS6KA5   107.058310        0.379362  0.458318  0.827725  0.407826  0.993036
    
    [1481 rows x 6 columns]


    ... done in 0.05 seconds.
    
    Calculating cook's distance...
    ... done in 0.00 seconds.
    
    Replacing 0 outlier genes.
    
    Running Wald tests...
    ... done in 0.05 seconds.
    





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
      <th>baseMean</th>
      <th>log2FoldChange</th>
      <th>lfcSE</th>
      <th>stat</th>
      <th>pvalue</th>
      <th>padj</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TRAC</th>
      <td>81.875820</td>
      <td>-0.930538</td>
      <td>0.542822</td>
      <td>-1.714261</td>
      <td>0.086481</td>
      <td>0.993036</td>
    </tr>
    <tr>
      <th>UCP2</th>
      <td>68.862556</td>
      <td>-0.548205</td>
      <td>0.422626</td>
      <td>-1.297142</td>
      <td>0.194582</td>
      <td>0.993036</td>
    </tr>
    <tr>
      <th>CMBL</th>
      <td>43.585692</td>
      <td>0.572860</td>
      <td>0.460769</td>
      <td>1.243269</td>
      <td>0.213769</td>
      <td>0.993036</td>
    </tr>
    <tr>
      <th>PJA2</th>
      <td>46.747220</td>
      <td>0.125016</td>
      <td>0.468412</td>
      <td>0.266893</td>
      <td>0.789552</td>
      <td>0.993036</td>
    </tr>
    <tr>
      <th>RNF125</th>
      <td>59.000601</td>
      <td>0.960666</td>
      <td>0.458653</td>
      <td>2.094537</td>
      <td>0.036212</td>
      <td>0.993036</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>LYRM7</th>
      <td>70.932088</td>
      <td>0.057458</td>
      <td>0.478781</td>
      <td>0.120008</td>
      <td>0.904476</td>
      <td>0.996294</td>
    </tr>
    <tr>
      <th>PNN</th>
      <td>116.380620</td>
      <td>0.213590</td>
      <td>0.417045</td>
      <td>0.512152</td>
      <td>0.608545</td>
      <td>0.993036</td>
    </tr>
    <tr>
      <th>CDC42EP3</th>
      <td>52.644710</td>
      <td>0.256493</td>
      <td>0.589972</td>
      <td>0.434755</td>
      <td>0.663740</td>
      <td>0.993036</td>
    </tr>
    <tr>
      <th>CHCHD2</th>
      <td>167.927239</td>
      <td>-0.497018</td>
      <td>0.449454</td>
      <td>-1.105825</td>
      <td>0.268802</td>
      <td>0.993036</td>
    </tr>
    <tr>
      <th>RPS6KA5</th>
      <td>107.058310</td>
      <td>0.379362</td>
      <td>0.458318</td>
      <td>0.827725</td>
      <td>0.407826</td>
      <td>0.993036</td>
    </tr>
  </tbody>
</table>
<p>1481 rows × 6 columns</p>
</div>



The results can be visualized using a volcano plot.


    
![png](rna_psbk_files/rna_psbk_37_0.png)
    


After performing DEA, we can use the resulting gene-level statistics for enrichment analysis.
While any statistic can be used, we recommend using `t-values` rather than `log2FCs`, as `t-values` account for the
significance of the change.
We will transform the obtained `t-values`, stored in the column `stat`, into a wide matrix format so that it can be
used with `decoupler`.




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
      <th>TRAC</th>
      <th>UCP2</th>
      <th>CMBL</th>
      <th>PJA2</th>
      <th>RNF125</th>
      <th>DNAJB6</th>
      <th>RIPOR2</th>
      <th>SERBP1</th>
      <th>RPL18AP3</th>
      <th>TM9SF3</th>
      <th>...</th>
      <th>CXCR4</th>
      <th>RPS19</th>
      <th>LIMD2</th>
      <th>CSNK1A1</th>
      <th>EIF3H</th>
      <th>LYRM7</th>
      <th>PNN</th>
      <th>CDC42EP3</th>
      <th>CHCHD2</th>
      <th>RPS6KA5</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>disease.vs.normal</th>
      <td>-1.714261</td>
      <td>-1.297142</td>
      <td>1.243269</td>
      <td>0.266893</td>
      <td>2.094537</td>
      <td>0.830098</td>
      <td>-1.210973</td>
      <td>-1.968079</td>
      <td>0.961525</td>
      <td>0.684806</td>
      <td>...</td>
      <td>1.351265</td>
      <td>-0.388153</td>
      <td>-0.25986</td>
      <td>-0.109398</td>
      <td>-0.566955</td>
      <td>0.120008</td>
      <td>0.512152</td>
      <td>0.434755</td>
      <td>-1.105825</td>
      <td>0.827725</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 1481 columns</p>
</div>



## Enrichment analysis

Enrichment analysis tests whether a specific set of omics features is "overrepresented" or "coordinated"
in the measured data compared to a background distribution. These sets are predefined based on existing
biological knowledge and may vary depending on the omics technology used.

Enrichment analysis requires the use of an enrichment method, and several options are available.
In the original manuscript of `decoupler` {cite:p}`decoupler`, we benchmarked multiple methods
and found that the univariate linear model (`ulm`) outperformed the others. Therefore, we will use
`ulm` in this vignette. 

The scores from {func}`decoupler.mt.ulm` should be interpreted such that larger magnitudes indicate
greater significance, while the sign reflects whether the features in the set are overrepresented
(positive) or underrepresented (negative) compared to the background.

### Transcription factor scoring from gene regulatory networks

Transcription factors (TFs) are genes that, once translated into proteins, bind to DNA and regulate
the expression of other genes by either promoting or inhibiting their transcription. Gene
Regulatory Networks (GRNs) capture these TF-gene interactions and can be constructed from prior
knowledge or inferred from omics data. The fundamental unit of a GRN is a TF and its associated target
genes, collectively known as a *regulon*. Each regulon functions as a gene set in enrichment analysis.

Although TFs are measured in transcriptomic data, their transcript levels often do not reflect their actual
activity in a given cell. Instead, scoring TFs through enrichment analysis based on the expression of
their target genes provides a more accurate representation of their regulatory activity
{cite:p}`grn_review`.

#### CollecTRI network
[CollecTRI](https://github.com/saezlab/CollecTRI) is a comprehensive resource containing a curated
collection of TFs and their transcriptional targets compiled from 12 different resources {cite:p}`collectri`. This collection
provides an increased coverage of transcription factors and a superior performance in identifying perturbed
TFs compared to other literature based GRNs such as 
[DoRothEA](https://saezlab.github.io/dorothea/) {cite:p}`dorothea`.
Similar to DoRothEA, interactions are weighted by their mode of regulation
(activation or inhibition).

In this tutorial we will use the human version but other organisms are available.
We can use `decoupler` to retrieve it from the [OmniPath](https://omnipathdb.org/) server {cite:p}`omnipath`.

<div class="alert alert-info">

**Note**
In this tutorial we use the network CollecTRI, but we could use any other GRN coming from an inference method such as [CellOracle](https://morris-lab.github.io/CellOracle.documentation/), [pySCENIC](https://pyscenic.readthedocs.io/en/latest/) or [SCENIC+](https://scenicplus.readthedocs.io/en/latest/). 
</div> 


    ---------------------------------------------------------------------------

    SSLError                                  Traceback (most recent call last)

    SSLError: [SSL: UNEXPECTED_EOF_WHILE_READING] EOF occurred in violation of protocol (_ssl.c:1010)

    
    The above exception was the direct cause of the following exception:


    MaxRetryError                             Traceback (most recent call last)

    File /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/requests/adapters.py:667, in HTTPAdapter.send(self, request, stream, timeout, verify, cert, proxies)
        666 try:
    --> 667     resp = conn.urlopen(
        668         method=request.method,
        669         url=url,
        670         body=request.body,
        671         headers=request.headers,
        672         redirect=False,
        673         assert_same_host=False,
        674         preload_content=False,
        675         decode_content=False,
        676         retries=self.max_retries,
        677         timeout=timeout,
        678         chunked=chunked,
        679     )
        681 except (ProtocolError, OSError) as err:


    File /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/urllib3/connectionpool.py:841, in HTTPConnectionPool.urlopen(self, method, url, body, headers, retries, redirect, assert_same_host, timeout, pool_timeout, release_conn, chunked, body_pos, preload_content, decode_content, **response_kw)
        839     new_e = ProtocolError("Connection aborted.", new_e)
    --> 841 retries = retries.increment(
        842     method, url, error=new_e, _pool=self, _stacktrace=sys.exc_info()[2]
        843 )
        844 retries.sleep()


    File /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/urllib3/util/retry.py:519, in Retry.increment(self, method, url, response, error, _pool, _stacktrace)
        518     reason = error or ResponseError(cause)
    --> 519     raise MaxRetryError(_pool, url, reason) from reason  # type: ignore[arg-type]
        521 log.debug("Incremented Retry for (url='%s'): %r", url, new_retry)


    MaxRetryError: HTTPSConnectionPool(host='zenodo.org', port=443): Max retries exceeded with url: /records/8192729/files/CollecTRI_regulons.csv?download=1 (Caused by SSLError(SSLEOFError(8, '[SSL: UNEXPECTED_EOF_WHILE_READING] EOF occurred in violation of protocol (_ssl.c:1010)')))

    
    During handling of the above exception, another exception occurred:


    SSLError                                  Traceback (most recent call last)

    Cell In[26], line 1
    ----> 1 collectri = dc.op.collectri(organism="human")
          2 collectri


    File /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/decoupler/op/_collectri.py:40, in collectri(organism, remove_complexes, license, verbose)
         18 """
         19 CollecTRI gene regulatory network :cite:p:`collectri`.
         20 
       (...)     37 and if available, the PMIDs supporting each interaction.
         38 """
         39 url = 'https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1'
    ---> 40 ct = _download(url, verbose=verbose)
         41 # Update resources
         42 resources = []


    File /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/decoupler/_download.py:21, in _download(url, verbose, **kwargs)
         19 _log(m, level='info', verbose=verbose)
         20 chunks = []
    ---> 21 with requests.get(url, stream=True) as r:
         22     r.raise_for_status()
         23     with tqdm(unit='B', unit_scale=True, desc="Progress", disable=not verbose) as pbar:


    File /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/requests/api.py:73, in get(url, params, **kwargs)
         62 def get(url, params=None, **kwargs):
         63     r"""Sends a GET request.
         64 
         65     :param url: URL for the new :class:`Request` object.
       (...)     70     :rtype: requests.Response
         71     """
    ---> 73     return request("get", url, params=params, **kwargs)


    File /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/requests/api.py:59, in request(method, url, **kwargs)
         55 # By using the 'with' statement we are sure the session is closed, thus we
         56 # avoid leaving sockets open which can trigger a ResourceWarning in some
         57 # cases, and look like a memory leak in others.
         58 with sessions.Session() as session:
    ---> 59     return session.request(method=method, url=url, **kwargs)


    File /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/requests/sessions.py:589, in Session.request(self, method, url, params, data, headers, cookies, files, auth, timeout, allow_redirects, proxies, hooks, stream, verify, cert, json)
        584 send_kwargs = {
        585     "timeout": timeout,
        586     "allow_redirects": allow_redirects,
        587 }
        588 send_kwargs.update(settings)
    --> 589 resp = self.send(prep, **send_kwargs)
        591 return resp


    File /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/requests/sessions.py:703, in Session.send(self, request, **kwargs)
        700 start = preferred_clock()
        702 # Send the request
    --> 703 r = adapter.send(request, **kwargs)
        705 # Total elapsed time of the request (approximately)
        706 elapsed = preferred_clock() - start


    File /opt/anaconda3/envs/scanpy/lib/python3.12/site-packages/requests/adapters.py:698, in HTTPAdapter.send(self, request, stream, timeout, verify, cert, proxies)
        694         raise ProxyError(e, request=request)
        696     if isinstance(e.reason, _SSLError):
        697         # This branch is for urllib3 v1.22 and later.
    --> 698         raise SSLError(e, request=request)
        700     raise ConnectionError(e, request=request)
        702 except ClosedPoolError as e:


    SSLError: HTTPSConnectionPool(host='zenodo.org', port=443): Max retries exceeded with url: /records/8192729/files/CollecTRI_regulons.csv?download=1 (Caused by SSLError(SSLEOFError(8, '[SSL: UNEXPECTED_EOF_WHILE_READING] EOF occurred in violation of protocol (_ssl.c:1010)')))


#### Scoring
TF scores can be easily computed by running the `ulm` method.


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    Cell In[27], line 2
          1 # Run
    ----> 2 tf_acts, tf_padj = dc.mt.ulm(data=data, net=collectri)
          4 # Filter by sign padj
          5 msk = (tf_padj.T < 0.05).iloc[:, 0]


    NameError: name 'collectri' is not defined


The obtained scores for the most active and inactive TFs can be visualized as follows.


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    Cell In[28], line 1
    ----> 1 dc.pl.barplot(data=tf_acts, name="disease.vs.normal", figsize=(4, 3))


    NameError: name 'tf_acts' is not defined


ATF3, KLF8 and ATF2 appear to be the most activated TFs in this treatment,
while GATA1, MYC and NFE2 appear to be inactivated.

A network of selected TFs (top and bottom ranked by activity) can also be visualized, with nodes colored by TF activity and
target gene expression.


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    Cell In[29], line 2
          1 dc.pl.network(
    ----> 2     net=collectri,
          3     data=data,
          4     score=tf_acts,
          5     sources=["ATF3", "KLF8", "GATA1", "MYC"],
          6     targets=5,
          7     figsize=(5, 5),
          8     vcenter=True,
          9     by_abs=True,
         10     size_node=15,
         11 )


    NameError: name 'collectri' is not defined


ATF3 appears to be active in disease cells, as its positively regulated targets are upregulated.

If needed, we can also look at a volcano plot of the target genes.


    
![png](rna_psbk_files/rna_psbk_49_0.png)
    


### Pathway Scoring

The same approach used for TF scoring can also be applied to pathways. Numerous
databases provide curated pathway gene sets, with one of the most well-known being MSigDB, which
includes several collections {cite:p}`msigdb`. 
These and many other resources can be accessed using the function {func}`decoupler.op.resource`.
To view the list of available databases, use {func}`decoupler.op.show_resources`.

#### PROGENy Pathway Genes
[PROGENy](https://saezlab.github.io/progeny/) is a comprehensive resource that provides a curated
collection of pathways and their target genes, along with weights for each interaction
{cite:p}`progeny`.

Below is a brief description of each pathway:

- **Androgen**: involved in the growth and development of the male reproductive organs
- **EGFR**: regulates growth, survival, migration, apoptosis, proliferation, and differentiation in mammalian cells
- **Estrogen**: promotes the growth and development of the female reproductive organs
- **Hypoxia**: promotes angiogenesis and metabolic reprogramming when O2 levels are low
- **JAK-STAT**: involved in immunity, cell division, cell death, and tumor formation
- **MAPK**: integrates external signals and promotes cell growth and proliferation
- **NFkB**: regulates immune response, cytokine production and cell survival
- **p53**: regulates cell cycle, apoptosis, DNA repair and tumor suppression
- **PI3K**: promotes growth and proliferation
- **TGFb**: involved in development, homeostasis, and repair of most tissues
- **TNFa**: mediates haematopoiesis, immune surveillance, tumour regression and protection from infection
- **Trail**: induces apoptosis
- **VEGF**: mediates angiogenesis, vascular permeability, and cell migration
- **WNT**: regulates organ morphogenesis during development and tissue repair

This is how to access to it.




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
      <th>source</th>
      <th>target</th>
      <th>weight</th>
      <th>padj</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Androgen</td>
      <td>TMPRSS2</td>
      <td>11.490631</td>
      <td>2.384806e-47</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Androgen</td>
      <td>NKX3-1</td>
      <td>10.622551</td>
      <td>2.205102e-44</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Androgen</td>
      <td>MBOAT2</td>
      <td>10.472733</td>
      <td>4.632376e-44</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Androgen</td>
      <td>KLK2</td>
      <td>10.176186</td>
      <td>1.944410e-40</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Androgen</td>
      <td>SARG</td>
      <td>11.386852</td>
      <td>2.790210e-40</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>62416</th>
      <td>p53</td>
      <td>ENPP2</td>
      <td>2.771405</td>
      <td>4.993215e-02</td>
    </tr>
    <tr>
      <th>62417</th>
      <td>p53</td>
      <td>ARRDC4</td>
      <td>3.494328</td>
      <td>4.996747e-02</td>
    </tr>
    <tr>
      <th>62418</th>
      <td>p53</td>
      <td>MYO1B</td>
      <td>-1.148057</td>
      <td>4.997905e-02</td>
    </tr>
    <tr>
      <th>62419</th>
      <td>p53</td>
      <td>CTSC</td>
      <td>-1.784693</td>
      <td>4.998864e-02</td>
    </tr>
    <tr>
      <th>62420</th>
      <td>p53</td>
      <td>NAA50</td>
      <td>-1.435013</td>
      <td>4.998884e-02</td>
    </tr>
  </tbody>
</table>
<p>62421 rows × 4 columns</p>
</div>



#### Scoring
Pathway scores can be readily computed by running the `ulm` method.




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
      <th>Estrogen</th>
      <th>Hypoxia</th>
      <th>NFkB</th>
      <th>PI3K</th>
      <th>TGFb</th>
      <th>TNFa</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>disease.vs.normal</th>
      <td>-3.008037</td>
      <td>3.779347</td>
      <td>4.047449</td>
      <td>-2.391001</td>
      <td>3.207506</td>
      <td>2.859247</td>
    </tr>
  </tbody>
</table>
</div>



The obtained scores for the most active and inactive pathways can be visualized as follows


    
![png](rna_psbk_files/rna_psbk_55_0.png)
    


They can also be visualized as follows


    
![png](rna_psbk_files/rna_psbk_57_0.png)
    


COVID patients exhibit increased activity of the NFkB, Hypoxia, TGFb, and TNFa pathways in T cells.

Conversely, the disease is associated with reduced activity in pathways such as Estrogen and PI3K.

Targets of the NFkB pathway can be visualized in a scatter plot.


    
![png](rna_psbk_files/rna_psbk_59_0.png)
    


The observed activation of the Hypoxia pathway is driven by the fact that many of its positively weighted
target genes have positive t-values.

We can also visualize the targets of NFkB as a leading edge plot.

Because PROGENy includes both positive and negative target genes,
the leading edge plot can be separated into positive and negative components.

    (+) leading edge: ['NR4A2' 'PLK3' 'P2RX7' 'JUN' 'SQSTM1']
    (-) leading edge: ['DCAF12' 'SNCA' 'ADIPOR1' 'ZMAT2' 'BSG']



    
![png](rna_psbk_files/rna_psbk_61_1.png)
    



    
![png](rna_psbk_files/rna_psbk_61_2.png)
    


#### Hallmark gene sets
[Hallmark](https://www.gsea-msigdb.org/gsea/msigdb/human/collection_details.jsp#H)
gene sets are curated collections of genes that represent specific, well-defined biological states or processes.
They are part of MSigDB and were developed to reduce redundancy and improve interpretability compared to older,
more overlapping gene set collections {cite:p}`msigdb`.

A total of 50 gene sets are provided, designed to be non-redundant, concise, and biologically coherent.

This is how to access them.




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
      <th>source</th>
      <th>target</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>IL2_STAT5_SIGNALING</td>
      <td>MAFF</td>
    </tr>
    <tr>
      <th>1</th>
      <td>COAGULATION</td>
      <td>MAFF</td>
    </tr>
    <tr>
      <th>2</th>
      <td>HYPOXIA</td>
      <td>MAFF</td>
    </tr>
    <tr>
      <th>3</th>
      <td>TNFA_SIGNALING_VIA_NFKB</td>
      <td>MAFF</td>
    </tr>
    <tr>
      <th>4</th>
      <td>COMPLEMENT</td>
      <td>MAFF</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>7313</th>
      <td>PANCREAS_BETA_CELLS</td>
      <td>STXBP1</td>
    </tr>
    <tr>
      <th>7314</th>
      <td>PANCREAS_BETA_CELLS</td>
      <td>ELP4</td>
    </tr>
    <tr>
      <th>7315</th>
      <td>PANCREAS_BETA_CELLS</td>
      <td>GCG</td>
    </tr>
    <tr>
      <th>7316</th>
      <td>PANCREAS_BETA_CELLS</td>
      <td>PCSK2</td>
    </tr>
    <tr>
      <th>7317</th>
      <td>PANCREAS_BETA_CELLS</td>
      <td>PAX6</td>
    </tr>
  </tbody>
</table>
<p>7318 rows × 2 columns</p>
</div>



#### Scoring
Pathway scores can be easily computed by running the `ulm` method.




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
      <th>HEME_METABOLISM</th>
      <th>INFLAMMATORY_RESPONSE</th>
      <th>MYC_TARGETS_V1</th>
      <th>OXIDATIVE_PHOSPHORYLATION</th>
      <th>REACTIVE_OXYGEN_SPECIES_PATHWAY</th>
      <th>TNFA_SIGNALING_VIA_NFKB</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>disease.vs.normal</th>
      <td>-5.764502</td>
      <td>2.88992</td>
      <td>-6.109761</td>
      <td>-5.686228</td>
      <td>-2.995499</td>
      <td>6.69648</td>
    </tr>
  </tbody>
</table>
</div>



The obtained scores for the most active and inactive pathways can be visualized as follows.


    
![png](rna_psbk_files/rna_psbk_67_0.png)
    


Or alternatively like this.


    
![png](rna_psbk_files/rna_psbk_69_0.png)
    


COVID patients exhibit increased activity for `TNFA_SIGNALING_VIA_NFKB`.

Conversely, the disease is associated with reduced activity in 
`MYC_TARGETS_V1`, `HEME_METABOLISM` and `OXIDATIVE_PHOSPHORYLATION`.

Targets of the `TNFA_SIGNALING_VIA_NFKB` can be visualized as a leading
edge plot.

    leading edge: ['NR4A2' 'JUN' 'SQSTM1' 'JUNB' 'IER5']



    
![png](rna_psbk_files/rna_psbk_71_1.png)
    


    /var/folders/dq/m9jbf2_x2_v6cnd1c0nt1rc40000gn/T/ipykernel_20072/3580451634.py:1: RuntimeWarning: Failed to import dependencies for application/vnd.jupyter.widget-view+json representation. (ModuleNotFoundError: No module named 'ipywidgets')
      sc.logging.print_versions()





<table class=table>
            <thead style="position: sticky; top: 0; background-color: var(--jp-layout-color0, var(--vscode-editor-background, white));">
        <tr><th>Package</th><th>Version</th></tr>
    </thead>
    <tbody>
        <tr><td><strong>pandas</strong></td><td>2.2.3</td></tr>
        <tr><td><strong>numpy</strong></td><td>2.2.5</td></tr>
        <tr><td><strong>scanpy</strong></td><td>1.11.1</td></tr>
        <tr><td><strong>decoupler</strong></td><td>2.0.0</td></tr>
        <tr><td><strong>anndata</strong></td><td>0.11.4</td></tr>
        <tr><td><strong>matplotlib</strong></td><td>3.10.3</td></tr>
        <tr><td><strong>pydeseq2</strong></td><td>0.5.3</td></tr>
    </tbody>
    <thead style="position: sticky; top: 0; background-color: var(--jp-layout-color0, var(--vscode-editor-background, white));">
        <tr><th>Component</th><th>Info</th></tr>
    </thead>
    <tbody>
        <tr><td>Python</td><td>3.12.9 | packaged by Anaconda, Inc. | (main, Feb  6 2025, 12:55:12) [Clang 14.0.6 ]</td></tr>
        <tr><td>OS</td><td>macOS-26.1-arm64-arm-64bit</td></tr>
        <tr><td>CPU</td><td>12 logical CPU cores, arm</td></tr>
        <tr><td>GPU</td><td>No GPU found</td></tr>
        <tr><td>Updated</td><td>2025-11-22 00:21</td></tr>
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
    <tr><td>llvmlite</td><td>0.44.0</td></tr>
    <tr><td>pytz</td><td>2025.2</td></tr>
    <tr><td>requests</td><td>2.32.3</td></tr>
    <tr><td>jedi</td><td>0.19.2</td></tr>
    <tr><td>patsy</td><td>1.0.1</td></tr>
    <tr><td>asttokens</td><td>3.0.0</td></tr>
    <tr><td>joblib</td><td>1.5.0</td></tr>
    <tr><td>pyzmq</td><td>26.2.0</td></tr>
    <tr><td>formulaic-contrasts</td><td>1.0.0</td></tr>
    <tr><td>pure-eval</td><td>0.2.2</td></tr>
    <tr><td>traitlets</td><td>5.14.3</td></tr>
    <tr><td>defusedxml</td><td>0.7.1</td></tr>
    <tr><td>typing_extensions</td><td>4.12.2</td></tr>
    <tr><td>packaging</td><td>24.2</td></tr>
    <tr><td>numba</td><td>0.61.2</td></tr>
    <tr><td>matplotlib-inline</td><td>0.1.6</td></tr>
    <tr><td>formulaic</td><td>1.2.1</td></tr>
    <tr><td>charset-normalizer</td><td>3.3.2</td></tr>
    <tr><td>session-info2</td><td>0.1.2</td></tr>
    <tr><td>threadpoolctl</td><td>3.6.0</td></tr>
    <tr><td>ipython</td><td>9.1.0</td></tr>
    <tr><td>cffi</td><td>1.17.1</td></tr>
    <tr><td>ptyprocess</td><td>0.7.0</td></tr>
    <tr><td>natsort</td><td>8.4.0</td></tr>
    <tr><td>certifi</td><td>2025.4.26 (2025.04.26)</td></tr>
    <tr><td>kiwisolver</td><td>1.4.8</td></tr>
    <tr><td>pycparser</td><td>2.21</td></tr>
    <tr><td>wrapt</td><td>2.0.1</td></tr>
    <tr><td>python-dateutil</td><td>2.9.0.post0</td></tr>
    <tr><td>executing</td><td>0.8.3</td></tr>
    <tr><td>marsilea</td><td>0.5.2</td></tr>
    <tr><td>PySocks</td><td>1.7.1</td></tr>
    <tr><td>jupyter_client</td><td>8.6.3</td></tr>
    <tr><td>docrep</td><td>0.3.2</td></tr>
    <tr><td>tqdm</td><td>4.67.1</td></tr>
    <tr><td>narwhals</td><td>2.12.0</td></tr>
    <tr><td>Pygments</td><td>2.19.1</td></tr>
    <tr><td>pyparsing</td><td>3.2.3</td></tr>
    <tr><td>tornado</td><td>6.4.2</td></tr>
    <tr><td>statsmodels</td><td>0.14.4</td></tr>
    <tr><td>Brotli</td><td>1.0.9</td></tr>
    <tr><td>PyYAML</td><td>6.0.2</td></tr>
    <tr><td>jupyter_core</td><td>5.7.2</td></tr>
    <tr><td>pyarrow</td><td>20.0.0</td></tr>
    <tr><td>igraph</td><td>0.11.8</td></tr>
    <tr><td>platformdirs</td><td>4.3.7</td></tr>
    <tr><td>setuptools</td><td>78.1.1</td></tr>
    <tr><td>stack-data</td><td>0.2.0</td></tr>
    <tr><td>comm</td><td>0.2.1</td></tr>
    <tr><td>interface-meta</td><td>1.3.0</td></tr>
    <tr><td>legacy-api-wrap</td><td>1.4.1</td></tr>
    <tr><td>pillow</td><td>11.2.1</td></tr>
    <tr><td>scikit-learn</td><td>1.5.2</td></tr>
    <tr><td>legendkit</td><td>0.3.4</td></tr>
    <tr><td>Cython</td><td>3.1.0</td></tr>
    <tr><td>prompt-toolkit</td><td>3.0.43</td></tr>
    <tr><td>urllib3</td><td>2.3.0</td></tr>
    <tr><td>leidenalg</td><td>0.10.2</td></tr>
    <tr><td>pexpect</td><td>4.8.0</td></tr>
    <tr><td>cycler</td><td>0.12.1</td></tr>
    <tr><td>scipy</td><td>1.15.3</td></tr>
    <tr><td>six</td><td>1.17.0</td></tr>
    <tr><td>adjustText</td><td>1.3.0</td></tr>
    <tr><td>appnope</td><td>0.1.3</td></tr>
    <tr><td>h5py</td><td>3.13.0</td></tr>
    <tr><td>ipykernel</td><td>6.29.5</td></tr>
    <tr><td>psutil</td><td>5.9.0</td></tr>
    <tr><td>texttable</td><td>1.7.0</td></tr>
    <tr><td>xgboost</td><td>3.0.1</td></tr>
    <tr><td>decorator</td><td>5.1.1</td></tr>
    <tr><td>wcwidth</td><td>0.2.5</td></tr>
    <tr><td>idna</td><td>3.7</td></tr>
    <tr><td>debugpy</td><td>1.8.11</td></tr>
    <tr><td>seaborn</td><td>0.13.2</td></tr>
</tbody>
    </table>
</div>
    </details>
        <details>
            <summary>Copyable Markdown</summary>
            <pre>| Package    | Version |
| ---------- | ------- |
| pandas     | 2.2.3   |
| numpy      | 2.2.5   |
| scanpy     | 1.11.1  |
| decoupler  | 2.0.0   |
| anndata    | 0.11.4  |
| matplotlib | 3.10.3  |
| pydeseq2   | 0.5.3   |

| Dependency          | Version                |
| ------------------- | ---------------------- |
| parso               | 0.8.4                  |
| llvmlite            | 0.44.0                 |
| pytz                | 2025.2                 |
| requests            | 2.32.3                 |
| jedi                | 0.19.2                 |
| patsy               | 1.0.1                  |
| asttokens           | 3.0.0                  |
| joblib              | 1.5.0                  |
| pyzmq               | 26.2.0                 |
| formulaic-contrasts | 1.0.0                  |
| pure-eval           | 0.2.2                  |
| traitlets           | 5.14.3                 |
| defusedxml          | 0.7.1                  |
| typing_extensions   | 4.12.2                 |
| packaging           | 24.2                   |
| numba               | 0.61.2                 |
| matplotlib-inline   | 0.1.6                  |
| formulaic           | 1.2.1                  |
| charset-normalizer  | 3.3.2                  |
| session-info2       | 0.1.2                  |
| threadpoolctl       | 3.6.0                  |
| ipython             | 9.1.0                  |
| cffi                | 1.17.1                 |
| ptyprocess          | 0.7.0                  |
| natsort             | 8.4.0                  |
| certifi             | 2025.4.26 (2025.04.26) |
| kiwisolver          | 1.4.8                  |
| pycparser           | 2.21                   |
| wrapt               | 2.0.1                  |
| python-dateutil     | 2.9.0.post0            |
| executing           | 0.8.3                  |
| marsilea            | 0.5.2                  |
| PySocks             | 1.7.1                  |
| jupyter_client      | 8.6.3                  |
| docrep              | 0.3.2                  |
| tqdm                | 4.67.1                 |
| narwhals            | 2.12.0                 |
| Pygments            | 2.19.1                 |
| pyparsing           | 3.2.3                  |
| tornado             | 6.4.2                  |
| statsmodels         | 0.14.4                 |
| Brotli              | 1.0.9                  |
| PyYAML              | 6.0.2                  |
| jupyter_core        | 5.7.2                  |
| pyarrow             | 20.0.0                 |
| igraph              | 0.11.8                 |
| platformdirs        | 4.3.7                  |
| setuptools          | 78.1.1                 |
| stack-data          | 0.2.0                  |
| comm                | 0.2.1                  |
| interface-meta      | 1.3.0                  |
| legacy-api-wrap     | 1.4.1                  |
| pillow              | 11.2.1                 |
| scikit-learn        | 1.5.2                  |
| legendkit           | 0.3.4                  |
| Cython              | 3.1.0                  |
| prompt-toolkit      | 3.0.43                 |
| urllib3             | 2.3.0                  |
| leidenalg           | 0.10.2                 |
| pexpect             | 4.8.0                  |
| cycler              | 0.12.1                 |
| scipy               | 1.15.3                 |
| six                 | 1.17.0                 |
| adjustText          | 1.3.0                  |
| appnope             | 0.1.3                  |
| h5py                | 3.13.0                 |
| ipykernel           | 6.29.5                 |
| psutil              | 5.9.0                  |
| texttable           | 1.7.0                  |
| xgboost             | 3.0.1                  |
| decorator           | 5.1.1                  |
| wcwidth             | 0.2.5                  |
| idna                | 3.7                    |
| debugpy             | 1.8.11                 |
| seaborn             | 0.13.2                 |

| Component | Info                                                                                |
| --------- | ----------------------------------------------------------------------------------- |
| Python    | 3.12.9 | packaged by Anaconda, Inc. | (main, Feb  6 2025, 12:55:12) [Clang 14.0.6 ] |
| OS        | macOS-26.1-arm64-arm-64bit                                                          |
| CPU       | 12 logical CPU cores, arm                                                           |
| GPU       | No GPU found                                                                        |
| Updated   | 2025-11-22 00:21                                                                    |</pre>
        </details>



    [NbConvertApp] Converting notebook rna_psbk.ipynb to markdown
    [NbConvertApp] Support files will be in rna_psbk_files/
    [NbConvertApp] Writing 33966 bytes to archive/rna_psbk.md


    [NbConvertApp] Converting notebook rna_psbk.ipynb to html
    [NbConvertApp] WARNING | Alternative text is missing on 20 image(s).
    [NbConvertApp] Writing 2460215 bytes to archive/rna_psbk.html

