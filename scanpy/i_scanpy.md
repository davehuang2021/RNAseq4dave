
### 重要的注意事项

- [ ] 从数据输入到完成聚类，中间经历若干步骤，务必确认加载数据类型（scale or log CPM or reads counts），以保证方法正确在使用
- [ ] **从.raw中提取**：QC, clustering等常用操作，因为涉及到基于选择后的基因，因此再利用数据时候，可能要经常从raw中提取，这个要注意！！
  - [ ] 未保存raw数据的h5ad，不会有这个数据结构
- [ ] scale后的数据，做表达量可视化会有误导性，需要注意。
  - [ ] 如何确定下载后的h5ad自带的数据是什么类型的？
- [ ] 下载的数据，得确认其数据格式（方法？）
  - [ ] 有必要的话，需要从矩阵开始重新分析
- [ ] 其他

### 代表性的项目

### 一个基本流程

### 重要的应用实例

#ref

- [NBIS scanpy](https://nbisweden.github.io/workshop-scRNAseq/home_contents.html): 
- [NBIS scRNAseq 2023](https://nbisweden.github.io/workshop-scrnaseq-2023/):这个和最新版本比，数据上差别不大
- [NBIS youtube](https://youtube.com/playlist?list=PLBsJUKzoJTHQA4Qg1yc1RRY2Km4t4vEeN&si=p37W8NwQREqf617q):
- [Pearson residuals for normalization of single-cell RNA-seq UMI data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02451-7):
- [scanpy Pseudo-bulk](https://decoupler-py.readthedocs.io/en/latest/notebooks/pseudobulk.html):

#notebook #GTDToolbox 以下文件添加这2个标签，用来管理任务和标记文件类型

### scanpy tasks

- [x] git文件夹在obsidian下好用的
- [ ] 下一步要测试在jupyter 系统中git 文件夹是否好用
