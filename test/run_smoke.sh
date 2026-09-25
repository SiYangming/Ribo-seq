#!/usr/bin/env bash
# chr20 test 模块冒烟：文件在位 + 脚本语法 + 可选二进制探测（不强制真跑）
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST="$ROOT/test"
FAIL=0

ok() { echo "  OK  $*"; }
bad() { echo "  FAIL $*"; FAIL=1; }

echo "==> [1/4] test 模块文件"
test -f "$TEST/info.csv" || bad "missing test/info.csv"
test -d "$TEST/GSE182201" || bad "missing test/GSE182201"
test -d "$TEST/reference" || bad "missing test/reference"
fastq_n=$(find "$TEST/GSE182201" -name '*.fastq.gz' | wc -l | tr -d ' ')
[ "$fastq_n" -ge 1 ] || bad "no FASTQ in test/GSE182201"
ok "info.csv + reference + ${fastq_n} FASTQ"

echo "==> [2/4] 参考关键文件（chr20）"
for f in \
  "GENCODE/v49/gencode.v49.pc_transcripts_chr20_filtered.fa" \
  "GENCODE/v49/gencode.v49.annotation_chr20.gtf" \
  "GENCODE/v49/transcript_info/gencode.v49.pc_transcripts_chr20_region_lengths.csv" \
  "rRNA/sortmerna_rrna.fasta" \
  "tRNA/hg38-mature-tRNAs-dna.fasta"
do
  if [ -f "$TEST/reference/$f" ]; then
    ok "$f"
  else
    bad "missing reference/$f"
  fi
done

echo "==> [3/4] bash -n 关键脚本 + run.sh --help"
bash -n "$ROOT/run.sh" || bad "run.sh syntax"
bash -n "$ROOT/Shell_scripts/common_variables.sh" || bad "common_variables.sh syntax"
n=0
while IFS= read -r -d '' sh; do
  bash -n "$sh" || bad "syntax $sh"
  n=$((n + 1))
done < <(find "$ROOT/Shell_scripts" -name '*.sh' -print0)
ok "checked $n Shell_scripts/*.sh"
help_out="$TEST/.help_out.txt"
"$ROOT/run.sh" --help >"$help_out" 2>&1 || true
if grep -Fq 'ref-suffix' "$help_out"; then
  ok "run.sh --help"
else
  bad "run.sh --help missing --ref-suffix"
fi
rm -f "$help_out"

echo "==> [4/4] 可选依赖探测（缺失仅警告）"
for bin in fastqc cutadapt umi_tools bowtie2 bbmap STAR rsem-calculate-expression samtools Rscript python3; do
  if command -v "$bin" >/dev/null 2>&1; then
    ok "found $bin"
  else
    echo "  WARN missing $bin (install via conda envs in *.yml)"
  fi
done

if [ "$FAIL" -ne 0 ]; then
  echo "SMOKE FAILED"
  exit 1
fi
echo "SMOKE PASSED (static). For a real run see test/README.md"
