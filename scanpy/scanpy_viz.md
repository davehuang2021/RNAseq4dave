
### scRNAseq gene exp corr

**虽然 z-score scale 不太会影响 R ，但可视化时候还是要看一下是否用的是 scale后的数据！！！**

```python
### 绘制 scRNAseq 基因表达相关性的图
sc.pl.scatter(
        Adata[
           (Adata.obs.cell_type=='tumor')
            &(Adata.obs.patient=='patient1'),], 
        x=gene1, 
        y=gene2,
        color='patient',
        palette=['#1A63BE'],
        size=50,show=False,title='',
    )
```

### scanpy 提取表达量到 obs

```python
### 本质上就是从adata.X中提取数据
Adata.obs[gene1] = list(Adata[: , [gene1]].X.T[0])
```

### scanpy 根据表达量来过滤细胞

```python
t_adata[((t_adata[: , [gene1_symbl]].X>0) & (t_adata[: , [gene2_symbl]].X>0)),:].obs
```

### scanpy densigy 绘图

- 这里注意，通常 scale 之后，会有在<0区间产生峰值；之前通常这个峰值在0的位置，也就是未测到表达量基因的信号。
- 这里的绘图是用seaborn实现的
```python
import seaborn as sns
sns.kdeplot(t_adata[: , gene_symbol].X.toarray())
### B2M 是普遍高表达的基因
sns.kdeplot(adata[: , 'B2M'].X.toarray())
### RPS27 是0位置和高表达都有双峰的基因
sns.kdeplot(adata[: , 'RPS27'].X.toarray())

```


### scanpy 整合已经做好的数据

#### tSNE 数据整合
有时候会单独提供可下载的table，而不是提供一个整合好的h5ad文件

```python
### 使用index来对齐，先用pandas整合进 obs
tsne = pd.read_csv('cluster.coords',sep='\t',skiprows=2,header=None).set_index(0)
tsne = tsne.rename(columns={1: 'Global_tSNE_1',2: 'Global_tSNE_2'})
adata.obs = pd.merge(adata.obs, tsne, left_index=True, right_index=True)

###生成替换用的数据
X_tsne_list = []
for i,j in adata.obs.iterrows():
    X_tsne_list.append([j.Global_tSNE_1,j.Global_tSNE_2])
#     X_tsne_list.append([j.Sub_tSNE_1,j.Sub_tSNE_2])
replace_X_tsne = np.array(X_tsne_list)
adata.obsm['X_tsne']=replace_X_tsne

### 检测一下画图结果对不对，测试时候得有个已经披露的数据可用
sc.pl.tsne(adata, color=['group'],size=60,ncols=3,wspace=0.5,palette='Set3')
```
