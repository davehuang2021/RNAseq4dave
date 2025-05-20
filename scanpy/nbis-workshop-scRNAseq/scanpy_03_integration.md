## CN note

数据集成的本质就是批次效应去除


## CN note: Data preparation

- The new variable gene selection should not be performed on the
scaled data matrix.
- 变异基因选择不能在缩放数据上做。这里记住之前提出的，缩放动作会覆盖表达量的这个事实。

## CN note: Detect variable genes

- 不处理批次效应，直接检测变异基因，会混入批次特意基因
- sc.pp.highly_variable_genes 的 batch_key 参数，在之前的分析中是未启用的
- Select all genes that are variable in at least 2 datasets and use for remaining analysis.
    - 选择至少 2 个数据集中可变的所有基因并用于剩余分析。这个解决方案就是排除掉单个批次的因素
- 注意流程到了哪个阶段：选择完 variable genes 后，把scale 和 PCA 做掉了
    - [ ] 这里看不懂了，**需要回顾**
- 这里还注意**保存旧umap&tSNE的方法**


## CN note: BBKNN 

知道有这个方法就行了。重点方法是 Harmony


## CN note: Harmony

- 最popular 的 batch effect remove 方法，但也就是它需要独立安装需要注意一下

效果如下：![scanpy harmony](scanpy_03_integration.png)

```python
for g in genes:
    plt.rcParams['axes.spines.left']   = True
    plt.rcParams['axes.spines.right']  = False
    plt.rcParams['axes.spines.top']    = False
    plt.rcParams['axes.spines.bottom'] = True
    
    plt.figure(figsize=(15, 5))
    
    colors = ['#1A63BE','#FBA002',  ]#  '#858CC4',
    # order = ['normal', 'EAC', 'ESCC']

    # 根据中位数排序
    value = g+", log2(TPM+1)"
    medians = exp.groupby('study')[value].median()
    sorted_categories = medians.sort_values(ascending=False).index
    print(sorted_categories[0:2])
    ax =sns.boxplot(y=value,x='study',data=exp[(~exp[value].apply(np.isinf))], 
                     palette=colors, 
                     order=sorted_categories,
                     showfliers=False, width=0.6
                    )
    interval = 2
    # plt.yticks(np.arange(0, np.max(exp[g+", log2TPM"]), interval))
    plt.yticks(np.arange(-4, 12, interval))

    plt.ylim(bottom=0,top=12)
    plt.xticks(rotation=90,fontsize=8)
    plt.title(g+' EXP in GTEx')
    ### 添加显著性检验的标记： Mann-Whitney： 
    pairs=[(sorted_categories[0], sorted_categories[1])]
    annotator = Annotator(ax, pairs, data=exp, x='study', y=g+", log2(TPM+1)", order=sorted_categories)
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside') ### 也叫Wilcoxon Rank Sum Test， 威尔科克森秩和检验
    annotator.apply_and_annotate()
    plt.savefig('.png', dpi=300, bbox_inches='tight')
    plt.show()
```
