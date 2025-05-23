编写规则参考[[habits2design_easy_to_use_pipeline]]的规则。

## 方法

- [ ] 先写大纲，然后填写代码 (在这个文件里)
- [ ] 写一个PPT把所有的关键节点串联一下
- 大纲根据笔记本上梳理的流程图来编写
- 先使用课程数据写一个，再整理一个实际应用的数据

## example dataet

### QC

主要针对细胞，去掉几种类型的有问题的细胞

### integration

这部分提前是基于雨不做这个能看到批次效应。所以其实也可以先做聚类。
如何检查批次效应是否存在：
- 直接做聚类，看是否有批次特异性的类
- 看批次在聚类中的分布
- 分析批次特异性的高可变基因

```
### 注意这里基本聚类的工作就做完了，是建立在新选择的高可变基因基础上的。回避了批次特异性基因。

sce.pp.harmony_integrate(adata, 'sample')
# Then we calculate a new umap and tsne.
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.tl.tsne(adata, use_rep='X_pca_harmony')
sc.tl.leiden(adata, resolution=0.5)

```

### clutering

- 在 integration 的基础上做聚类
- 存在多种聚类方法，多尝试（问题：如何判断使用那种方法？）
  - 各种聚类参数可以添加标签来标记，可以多尝试，比较不通的效果
- 对细胞亚群，单独进行聚类要重新定义差异基因和恢复原始数据，也就是从.raw里重新提取。
  - integration也要重新做，本质上就是对子集做重新分群
  - 问题，结果如何整合在一起？

### DEG

### pseudoBulk

## real dataset 1

