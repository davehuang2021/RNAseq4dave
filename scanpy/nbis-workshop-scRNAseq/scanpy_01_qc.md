- [ref:NBIS scanpy_01_qc](https://nbisweden.github.io/workshop-scRNAseq/labs/scanpy/scanpy_01_qc.html)

## 小结与体会

- 前置步骤还包括数据的整合和补充meta信息
- QC的步骤的关键在于去掉不好的细胞和样品，如高线粒体比例细胞，低核糖体比例的，血红蛋白污染，doublets样品等等
- 细胞检出的基因数量控制
- 过滤之后保存数据，这里特别注意要保存raw，因为这个步骤得到的数据是必定要后续使用的

## CN note： Get data
注意这里样品的文件名已经提示了样品的类型和编号

## CN note： Collate （整合数据）

- 注意合并数据的方法是很简单的，只要adata.concatenate即可
- 这里删除数据之前，记得先打印一下看看格式，熟悉熟悉一下数据结构


## CN note：Calculate QC

- 每个细胞的线粒体和核糖体基因百分比
    - 于穿孔细胞导致胞质 RNA 丢失。其原因在于，线粒体比单个转录分子更大。因此濒死细胞或低质量，线粒体转录本占比高
- 血红蛋白基因的比例可以指示红细胞的污染情况

- 先定义哪些是线粒体基因，ribosomal genes和hemoglobin genes.


## CN note： Plot QC

质量控制的绘图，能看到一些分布，但不提供确定的QC cut-off

## CN note: Filtering

- 至少检测到 200 个基因且基因需要在至少 3 个细胞中表达的细胞
- Extremely high number of detected genes could indicate doublets, 检测到的基因数量极高可能表明存在双联体

### doublets 预测是单独的notebook

## CN note: Mito/Ribo filtering

- pct_counts_mt 低较好，细胞无转录组渗漏
- pct_counts_ribo不能太低，但为什么不清楚


As the ribosomal proteins are highly expressed they will make up a larger proportion of the transcriptional landscape when fewer of the lowly expressed genes are detected。
按理说，高核糖体基因比例，显示了对低表达基因的检出情况较差。但为什么要维持在一个高的比例呢，有低表达基因的量可能是和doublets正相关？

- 这里的过滤方式直接用了python数据结构来操作


## CN note: Filter genes

### 保留的特征基因

- 这里做了移除基因的操作，原因在于他们的差异源于技术因素，且对功能差异的贡献较小。并且由于细胞穿孔，红细胞污染等因素，也会引入bias，因此移除这些基因。

### 性别偏差

- 性别差异当然会引入影响。但这里只给出了分析方法



## CN note: Predict doublets

在典型的 10 倍实验中，双联体的比例与上样细胞的数量呈线性关系。因此这是个根据趋势来预测的工作。

Most doublet detectors simulates doublets by merging cell counts and predicts doublets as cells that have similar embeddings as the simulated doublets. Most such packages need an assumption about the number/proportion of expected doublets in the dataset. The data you are using is subsampled, but the original datasets contained about **5 000 cells per sample**, hence we can assume that **they loaded about 9 000 cells** and should have a doublet rate at about 4%.

这里的假设是需要分析师来确认伤样的模式，即单次上样的细胞量。这里缺失了一张表的截图，展示 loaded call count 和 doublet 比例的对应关系。

这里要用到批次 batch 信息，这里对应的是 sample。

- 我们应该预期两个细胞比单个细胞检测到更多的基因，让我们检查一下我们预测的双联体是否总体上也检测到更多的基因。

这个信息可以画图展示


### 移除所有的预测 doublets 细胞并保存数据

