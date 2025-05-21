
```bash
# 保存文件夹结构

cd ./scanpy/nbis-workshop-scRNAseq/

tree data >data-dir_files.md

# 计算md5用于检查文件的一致性
for i in `find ./data `;do;md5sum $i;done >data-dir_files-md5.md

# data-dir_files-cmd.sh.md用于保存以上命令和其他乱七八糟的命令
```