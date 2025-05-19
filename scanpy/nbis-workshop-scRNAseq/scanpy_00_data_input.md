## data

- [nbis data](https://nbisweden.github.io/workshop-scRNAseq/other/data.html)


## scanpy read_10x_mtx

- read_10x_mtx 的 prefix 参数用于添加 GEOID这种前置的ID，简化编辑步骤
- gz文件无需解压即可读取（建议不要解压缩）

```python
### 非常简单，几本和软件说明书中的写法一样
adata = sc.read_10x_mtx(
    path="./data",       # 文件所在目录
    prefix="GSE149689_",
    var_names="gene_symbols",  # 用 gene symbol 作为基因名（features.tsv 的第2列）
    cache=True           # 是否缓存 .h5ad 格式（可加快后续读取）
)
```
### scanpy 教程数据

[ref](https://scanpy.readthedocs.io/en/stable/tutorials/experimental/pearson_residuals.html)

数据可以wget直接下载。
- [x] 经过了检查，没有和已经发表的数据对应起来。放弃研究可重复性（2025-5-19）

### GSE149689 实例

[GSE149689是实例数据来源](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689)，注意格式是mtx形式 

ATATCCTCAGGTTTAC-19	Ctrl	ctrl_19，检查barcode 19 表示了样品的编号
但还是需要仔细研究来确认

### GSE178318 CRC data (barcode, gene, mtx), pmid34489408, liver metasits

### GSE132465 nature genetics data, SMC

### GSE144735 nature genetics data, 6 Belgian colorectal cancer patients

### GSE136831 IPF data

### GSE159354 SscILD h5文件格式

### 
```
sceasy::convertFormat(data2, from="seurat", to="anndata",
                       outFile='GSE35893_data.h5ad')

轻松转化格式，也可以参考nbis数据原始处理与筛选的脚本，是否可以直接用R生成
```

### 更多实例，需要回去翻一下我工作站里的工作目录
