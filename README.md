# ezbiocloud
A python script for uploading multiple 16S rRNA sequences to ezbiocloud.
# 安装
> ### 克隆仓库
```shell
#### 克隆仓库到本地 ####
git clone https://github.com/1996xjm/ezbiocloud.git

#### 安装python依赖 ####
# python版本最好3.7以上

# 安装 biopython
pip install biopython

```

# 使用


```shell
Usage: python ezbiocloud.py -f/--fastaFile <fasta file> -o/--outputPath <output annotation file path>

Example: python ezbiocloud.py -f example.fasta -o ./

Options:
  -h, --help            show this help message and exit
  -f FASTAFILE, --fastaFile=FASTAFILE
                        fastaFile
  -o OUTPUTPATH, --outputPath=OUTPUTPATH
                        output annotation file path
```
