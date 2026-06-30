#!/usr/bin/env python3
"""Build reciprocal-best-hit gene map from two Prokka annotations via tBLASTn."""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path


def parse_fasta_headers(faa_path: Path) -> list[str]:
    ids: list[str] = []
    with faa_path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def parse_tblastn_best_hits(
    blast_path: Path,
    min_pident: float,
    min_qcovs: float,
) -> dict[str, tuple[str, float, float, float]]:
    """Return query_id -> (subject_id, bitscore, pident, qcovs) for top hit per query."""
    hits_by_query: dict[str, list[tuple[str, float, float, float, float]]] = defaultdict(list)

    with blast_path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            qseqid, sseqid, pident, length, qstart, qend, sstart, send, evalue, bitscore = fields[:10]
            qlen = int(fields[10]) if len(fields) > 10 else None
            slen = int(fields[11]) if len(fields) > 11 else None
            qcovs = float(fields[12]) if len(fields) > 12 else None

            if qcovs is None and qlen:
                qcov = abs(int(qend) - int(qstart)) + 1
                qcovs = 100.0 * qcov / qlen

            pident_f = float(pident)
            bitscore_f = float(bitscore)
            qcovs_f = float(qcovs) if qcovs is not None else 0.0

            if pident_f < min_pident or qcovs_f < min_qcovs:
                continue

            hits_by_query[qseqid].append((sseqid, bitscore_f, pident_f, qcovs_f, float(evalue)))

    best: dict[str, tuple[str, float, float, float]] = {}
    for qseqid, hit_list in hits_by_query.items():
        hit_list.sort(key=lambda item: (-item[1], item[4], -item[2], -item[3]))
        subject, bitscore, pident, qcovs, _ = hit_list[0]
        best[qseqid] = (subject, bitscore, pident, qcovs)
    return best


def run_tblastn(
    query_faa: Path,
    subject_ffn: Path,
    out_path: Path,
    num_threads: int,
) -> None:
    db_prefix = out_path.with_suffix("")
    subprocess.run(
        ["makeblastdb", "-in", str(subject_ffn), "-dbtype", "nucl", "-out", str(db_prefix)],
        check=True,
    )
    subprocess.run(
        [
            "tblastn",
            "-query",
            str(query_faa),
            "-db",
            str(db_prefix),
            "-out",
            str(out_path),
            "-outfmt",
            "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qcovs",
            "-evalue",
            "1e-5",
            "-max_target_seqs",
            "5",
            "-num_threads",
            str(num_threads),
        ],
        check=True,
    )


def load_prokka_tsv(tsv_path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with tsv_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["locus_tag"]] = row
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-prefix", required=True, help="Label for query genome (e.g. ATCC8482)")
    parser.add_argument("--subject-prefix", required=True, help="Label for subject genome (e.g. CL09)")
    parser.add_argument("--query-faa", type=Path, required=True)
    parser.add_argument("--subject-faa", type=Path, required=True)
    parser.add_argument("--query-ffn", type=Path, required=True)
    parser.add_argument("--subject-ffn", type=Path, required=True)
    parser.add_argument("--query-tsv", type=Path, required=True)
    parser.add_argument("--subject-tsv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--map-out", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-pident", type=float, default=70.0)
    parser.add_argument("--min-qcovs", type=float, default=70.0)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    forward_blast = args.outdir / f"{args.query_prefix}_to_{args.subject_prefix}.tblastn.tsv"
    reverse_blast = args.outdir / f"{args.subject_prefix}_to_{args.query_prefix}.tblastn.tsv"

    run_tblastn(args.query_faa, args.subject_ffn, forward_blast, args.threads)
    run_tblastn(args.subject_faa, args.query_ffn, reverse_blast, args.threads)

    forward_best = parse_tblastn_best_hits(forward_blast, args.min_pident, args.min_qcovs)
    reverse_best = parse_tblastn_best_hits(reverse_blast, args.min_pident, args.min_qcovs)

    query_meta = load_prokka_tsv(args.query_tsv)
    subject_meta = load_prokka_tsv(args.subject_tsv)

    reciprocal_rows: list[dict[str, str | float]] = []
    for query_gene, (subject_gene, f_bitscore, f_pident, f_qcovs) in forward_best.items():
        reverse_hit = reverse_best.get(subject_gene)
        if reverse_hit is None:
            continue
        reverse_query, r_bitscore, r_pident, r_qcovs = reverse_hit
        if reverse_query != query_gene:
            continue

        qrow = query_meta.get(query_gene, {})
        srow = subject_meta.get(subject_gene, {})
        reciprocal_rows.append(
            {
                f"{args.query_prefix}_locus_tag": query_gene,
                f"{args.subject_prefix}_locus_tag": subject_gene,
                f"{args.query_prefix}_gene": qrow.get("gene", ""),
                f"{args.subject_prefix}_gene": srow.get("gene", ""),
                f"{args.query_prefix}_product": qrow.get("product", ""),
                f"{args.subject_prefix}_product": srow.get("product", ""),
                "forward_bitscore": f_bitscore,
                "reverse_bitscore": r_bitscore,
                "forward_pident": f_pident,
                "reverse_pident": r_pident,
                "forward_qcovs": f_qcovs,
                "reverse_qcovs": r_qcovs,
            }
        )

    reciprocal_rows.sort(key=lambda row: row[f"{args.query_prefix}_locus_tag"])

    fieldnames = list(reciprocal_rows[0].keys()) if reciprocal_rows else [
        f"{args.query_prefix}_locus_tag",
        f"{args.subject_prefix}_locus_tag",
        f"{args.query_prefix}_gene",
        f"{args.subject_prefix}_gene",
        f"{args.query_prefix}_product",
        f"{args.subject_prefix}_product",
        "forward_bitscore",
        "reverse_bitscore",
        "forward_pident",
        "reverse_pident",
        "forward_qcovs",
        "reverse_qcovs",
    ]

    with args.map_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(reciprocal_rows)

    print(f"Forward best hits: {len(forward_best)}")
    print(f"Reverse best hits: {len(reverse_best)}")
    print(f"Reciprocal best hits: {len(reciprocal_rows)}")
    print(f"Gene map written to {args.map_out}")


if __name__ == "__main__":
    main()
