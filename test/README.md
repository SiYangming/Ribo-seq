# test/ — chr20 端到端冒烟模块

本目录是整套流程的**环境与链路自检**数据集（GSE182201 的 chr20 子集 + 配套参考），用于验证当前机器依赖是否齐全、脚本能否跑通。

| 路径 | 内容 |
|------|------|
| `info.csv` | 样本表（fastq 路径相对仓库根） |
| `GSE182201/` | chr20 子集 FASTQ + `make_test_data.sh` |
| `reference/` | chr20 GENCODE / rRNA / tRNA（不含 STAR/RSEM 索引，运行时构建） |

## 一键自检

```bash
# 静态：语法 + 测试文件在位 + 可选依赖探测
bash test/run_smoke.sh

# 真跑 Totals 冒烟（需 conda 环境；首次会建索引，较慢）
./run.sh --pipeline Totals --list-steps
./run.sh --pipeline Totals \
  --input-csv test/info.csv \
  --fasta-dir test/reference \
  --genome-version v49 \
  --ref-suffix _chr20 \
  --output-dir results/test_smoke \
  --threads 4
```

环境变量等价写法：

```bash
export RIBO_SEQ_FASTA_DIR="$PWD/test/reference"
export RIBO_SEQ_REF_SUFFIX=_chr20
export GENOME_VERSION=v49
```

## 网络补全 FASTQ（若本地缺文件）

```bash
bash scripts/download_testdata.sh
# 默认写入 test/GSE182201/，按 test/info.csv（或根目录 info.csv）列名拉取 nf-core/test-datasets
```

## 说明

- 生产/论文数据请另建工作目录，不要覆盖 `test/`。
- 索引目录 `STAR_index/`、`rsem_bowtie2_index/` 由 `Shell_scripts/check_and_build_indices.sh` 生成，已 gitignore。
