"""Run Rearrangement Clustering without comparison algorithms."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import time
from pathlib import Path

import numpy as np

from core import fingerprint, load_or_build_graph, read_dataset, rearrange_sparse


ROOT = Path(__file__).resolve().parent
DATASETS = ("16s_01", "16s_02", "18s_01", "18s_02", "coi_01", "coi_02")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=DATASETS, default="16s_01")
    parser.add_argument("--data-dir", type=Path, help="Custom directory containing <name>.fasta and <name>.csv")
    parser.add_argument("--threshold", type=float, default=0.95, help="Similarity threshold in (0, 1); default: 0.95")
    parser.add_argument("--threads", type=int, default=8, help="Graph construction workers; default: 8")
    parser.add_argument("--cache-dir", type=Path, default=ROOT / "cache")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--vsearch", default="vsearch", help="VSEARCH executable used when building a graph")
    args = parser.parse_args()
    if not 0 < args.threshold < 1:
        parser.error("--threshold must be in (0, 1)")
    if args.threads < 1:
        parser.error("--threads must be positive")

    directory = (args.data_dir or ROOT / "data" / args.dataset).resolve()
    name = directory.name
    started = time.perf_counter()
    records = read_dataset(directory)
    print(f"Dataset: {name}; sequences: {len(records):,}; threshold: {args.threshold:g}", flush=True)
    graph = load_or_build_graph(
        records, args.cache_dir / f"{name}_sparse_neighbors.npz",
        executable=shutil.which(args.vsearch) or args.vsearch, threads=args.threads,
    )
    graph_seconds = time.perf_counter() - started
    result = rearrange_sparse(graph, args.threshold)
    output = args.output_dir or ROOT / "results" / name / f"theta_{args.threshold:g}"
    output.mkdir(parents=True, exist_ok=True)

    with (output / "assignments.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("seq_id", "cluster_id"))
        writer.writerows((record.sequence_id, int(label)) for record, label in zip(records, result.labels))
    with (output / "ordering.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("position", "seq_id", "cluster_id"))
        writer.writerows((position, records[int(node)].sequence_id, int(result.labels[node]))
                         for position, node in enumerate(result.order))
    with (output / "edges.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("left_seq_id", "right_seq_id", "distance", "similarity", "is_cut"))
        writer.writerows((records[int(left)].sequence_id, records[int(right)].sequence_id,
                          float(gap), 1 - float(gap), int(cut))
                         for left, right, gap, cut in zip(
                             result.order[:-1], result.order[1:], result.gaps, result.breaks))
    summary = {
        "dataset": name,
        "input_fingerprint": fingerprint(records),
        "sequence_count": len(records),
        "cluster_count": int(result.breaks.sum() + 1),
        "threshold": args.threshold,
        "distance_cutoff": 1 - args.threshold,
        "retained_path_edges": int((~result.breaks).sum()),
        "cut_path_edges": int(result.breaks.sum()),
        "objective": float(np.minimum(result.gaps, 1 - args.threshold).sum()),
        "graph_load_or_build_seconds": graph_seconds,
        "clustering_seconds": result.runtime_seconds,
        "total_seconds": time.perf_counter() - started,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"Clusters: {summary['cluster_count']:,}; clustering: {result.runtime_seconds:.2f}s")
    print(f"Output: {output.resolve()}")


if __name__ == "__main__":
    main()
