# Online enrichment helpers (optional; Downstream_5 still defaults to local fgsea)

| Script | Backend | Install |
|--------|---------|---------|
| [Enrichr.R](Enrichr.R) | Maayanlab Enrichr (`enrichR` R package) | `install.packages("enrichR")` |
| [gProfiler.R](gProfiler.R) | [gprofiler2](https://biit.cs.ut.ee/gprofiler/page/r) | `install.packages("gprofiler2")` |
| [../../Python_scripts/gProfiler.py](../../Python_scripts/gProfiler.py) | [g:Profiler API](https://biit.cs.ut.ee/gprofiler/page/apis) / `gprofiler-official` | `pip install gprofiler-official requests` |
| [fgsea.R](fgsea.R) | Local GMT (offline) | Bioconductor `fgsea` |

## Examples

```bash
# Enrichr (gene list)
Rscript R_scripts/gsea/Enrichr.R --genes my_genes.txt --treatment demo

# Enrichr from DESeq2 TE groups
Rscript R_scripts/gsea/Enrichr.R \
  --deseq2 results/Analysis/DESeq2_output/X_merged_DESeq2.csv \
  --treatment X --padj 0.1

# g:Profiler R
Rscript R_scripts/gsea/gProfiler.R \
  --genes my_genes.txt --organism hsapiens --treatment demo

# g:Profiler Python (REST or official client)
python3 Python_scripts/gProfiler.py \
  --genes my_genes.txt --organism hsapiens --treatment demo
```

Outputs default to `$RIBO_SEQ_PARENT_DIR/Analysis/{Enrichr,gProfiler}/` (or `./results/...`).
