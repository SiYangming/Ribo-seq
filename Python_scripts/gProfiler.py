#!/usr/bin/env python3
"""gProfiler.py — functional enrichment via g:Profiler REST API / official client.

API docs: https://biit.cs.ut.ee/gprofiler/page/apis
Python client: pip install gprofiler-official  (preferred)
Fallback: plain requests POST to /api/gost/profile/

Examples:
  python Python_scripts/gProfiler.py --genes genes.txt --organism hsapiens
  python Python_scripts/gProfiler.py --deseq2 merged_DESeq2.csv --te-group "TE up" \\
      --organism hsapiens --outdir results/Analysis/gProfiler
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Any

GPROFILER_GOST_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
DEFAULT_SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="g:Profiler g:GOSt enrichment (Python)")
    p.add_argument("--genes", help="Text file: one gene ID/symbol per line")
    p.add_argument("--deseq2", help="DESeq2 merged CSV with gene_sym (+ optional TE_group)")
    p.add_argument("--te-group", help="Filter DESeq2 rows by TE_group, e.g. 'TE up'")
    p.add_argument("--organism", default=os.environ.get("RIBO_SEQ_ORGANISM", "hsapiens"))
    p.add_argument(
        "--sources",
        default=",".join(DEFAULT_SOURCES),
        help="Comma-separated source IDs (default: GO:BP,GO:MF,GO:CC,KEGG,REAC,WP)",
    )
    p.add_argument("--user-threshold", type=float, default=0.05)
    p.add_argument(
        "--correction",
        default="g_SCS",
        choices=("g_SCS", "fdr", "bonferroni"),
    )
    p.add_argument(
        "--treatment",
        default=os.environ.get("RIBO_SEQ_TREATMENT", "run"),
        help="Label prefix for output files",
    )
    p.add_argument(
        "--outdir",
        default=None,
        help="Output directory (default: $RIBO_SEQ_PARENT_DIR/Analysis/gProfiler or ./results/...)",
    )
    p.add_argument(
        "--api-url",
        default=GPROFILER_GOST_URL,
        help="g:GOSt API endpoint",
    )
    p.add_argument(
        "--use-client",
        action="store_true",
        help="Force gprofiler-official client if installed",
    )
    return p.parse_args()


def read_gene_file(path: Path) -> list[str]:
    genes: list[str] = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        genes.append(s)
    # unique, preserve order
    seen: set[str] = set()
    out: list[str] = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def read_deseq2(path: Path, te_group: str | None) -> dict[str, list[str]]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or "gene_sym" not in rows[0]:
        raise SystemExit(f"[ERROR] {path} must contain gene_sym column")

    def collect(group: str | None) -> list[str]:
        seen: set[str] = set()
        out: list[str] = []
        for r in rows:
            if group is not None and r.get("TE_group") != group:
                continue
            g = (r.get("gene_sym") or "").strip()
            if g and g not in seen:
                seen.add(g)
                out.append(g)
        return out

    if te_group:
        label = "".join(c if c.isalnum() else "_" for c in te_group)
        return {label: collect(te_group)}

    out: dict[str, list[str]] = {}
    for grp in ("TE up", "TE down"):
        genes = collect(grp)
        if genes:
            label = "".join(c if c.isalnum() else "_" for c in grp)
            out[label] = genes
    if not out:
        out["all_genes"] = collect(None)
    return out


def load_queries(args: argparse.Namespace) -> dict[str, list[str]]:
    if args.genes:
        return {"default": read_gene_file(Path(args.genes))}
    if args.deseq2:
        return read_deseq2(Path(args.deseq2), args.te_group)
    raise SystemExit("[ERROR] Provide --genes or --deseq2")


def gost_via_client(
    genes: list[str],
    organism: str,
    sources: list[str],
    user_threshold: float,
    correction: str,
) -> list[dict[str, Any]]:
    from gprofiler import GProfiler  # type: ignore

    gp = GProfiler(return_all=False)
    return gp.profile(
        organism=organism,
        query=genes,
        sources=sources,
        user_threshold=user_threshold,
        significance_threshold_method=correction,
    )


def gost_via_requests(
    genes: list[str],
    organism: str,
    sources: list[str],
    user_threshold: float,
    correction: str,
    api_url: str,
) -> list[dict[str, Any]]:
    payload = {
        "organism": organism,
        "query": genes,
        "sources": sources,
        "user_threshold": user_threshold,
        "significance_threshold_method": correction,
        "no_evidences": True,
    }
    headers = {
        "User-Agent": "SiYangming-Ribo-seq/gProfiler.py",
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    body = json.dumps(payload).encode("utf-8")

    try:
        import requests

        r = requests.post(api_url, json=payload, headers=headers, timeout=120)
        r.raise_for_status()
        data = r.json()
    except ImportError:
        import urllib.error
        import urllib.request

        req = urllib.request.Request(api_url, data=body, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                data = json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            raise SystemExit(f"[ERROR] g:Profiler HTTP {e.code}: {e.read()[:500]!r}") from e
    except Exception as e:
        raise SystemExit(f"[ERROR] g:Profiler request failed: {e}") from e

    return data.get("result") or []


def flatten_row(row: dict[str, Any]) -> dict[str, Any]:
    flat: dict[str, Any] = {}
    for k, v in row.items():
        if isinstance(v, (list, dict)):
            flat[k] = json.dumps(v, ensure_ascii=False)
        else:
            flat[k] = v
    return flat


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    flat_rows = [flatten_row(r) for r in rows]
    # union of keys
    fieldnames: list[str] = []
    seen: set[str] = set()
    for r in flat_rows:
        for k in r:
            if k not in seen:
                seen.add(k)
                fieldnames.append(k)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(flat_rows)


def main() -> int:
    args = parse_args()
    parent = os.environ.get("RIBO_SEQ_PARENT_DIR") or "results"
    outdir = Path(args.outdir or os.path.join(parent, "Analysis", "gProfiler"))
    outdir.mkdir(parents=True, exist_ok=True)

    sources = [s.strip() for s in args.sources.split(",") if s.strip()]
    queries = load_queries(args)

    use_client = args.use_client
    if not use_client:
        try:
            import gprofiler  # noqa: F401

            use_client = True
        except ImportError:
            use_client = False

    for label, genes in queries.items():
        print(f"[INFO] g:Profiler [{label}]: {len(genes)} genes; organism={args.organism}")
        if len(genes) < 3:
            print(f"[WARN] skip {label}: <3 genes")
            continue
        if use_client:
            try:
                result = gost_via_client(
                    genes, args.organism, sources, args.user_threshold, args.correction
                )
            except Exception as e:
                print(f"[WARN] gprofiler-official failed ({e}); falling back to REST")
                result = gost_via_requests(
                    genes,
                    args.organism,
                    sources,
                    args.user_threshold,
                    args.correction,
                    args.api_url,
                )
        else:
            result = gost_via_requests(
                genes,
                args.organism,
                sources,
                args.user_threshold,
                args.correction,
                args.api_url,
            )

        # gprofiler-official returns list[dict]; REST same
        if not isinstance(result, list):
            result = list(result) if result else []

        prefix = f"{args.treatment}_{label}"
        out_csv = outdir / f"{prefix}_gost.csv"
        write_csv(out_csv, result)
        print(f"[INFO] wrote {out_csv} ({len(result)} terms)")

    print(f"[INFO] gProfiler.py done → {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
