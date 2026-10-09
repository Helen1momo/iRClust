"""Small standalone checks; no external executables are needed."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from core import SparseGraph, _two_opt_sparse, fingerprint, read_dataset


ROOT = Path(__file__).resolve().parent


class ClusteringTests(unittest.TestCase):
    def test_cached_cli_outputs_and_threshold_boundary(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data, cache, output = root / "toy", root / "cache", root / "output"
            data.mkdir()
            cache.mkdir()
            (data / "toy.fasta").write_text(
                ">a\nAAAA\n>b\nAAAT\n>c\nTTTT\n>d\nTTTA\n")
            (data / "toy.csv").write_text(
                "seq_id,class,family,genus\na,A,A,A\nb,A,A,A\nc,B,B,B\nd,B,B,B\n")
            records = read_dataset(data)
            np.savez_compressed(cache / "toy_sparse_neighbors.npz",
                                edge_u=[0, 0, 0, 1, 1, 2], edge_v=[1, 2, 3, 2, 3, 3],
                                edge_distance=np.asarray([.125, .5, .5, .5, .5, .125], dtype=np.float32),
                                features=np.eye(4, dtype=np.float32))
            metadata = {"fingerprint": fingerprint(records), "sketch_k": 8,
                        "sketch_size": 256, "candidate_k": 128, "exact_k": 128}
            (cache / "toy_sparse_neighbors.json").write_text(json.dumps(metadata))
            command = [sys.executable, str(ROOT / "run_clustering.py"),
                       "--data-dir", str(data), "--cache-dir", str(cache),
                       "--output-dir", str(output), "--vsearch", "missing-vsearch-test"]
            for threshold, expected_k in (("0.75", 2), ("0.875", 4)):
                completed = subprocess.run(command + ["--threshold", threshold],
                                           cwd=root, capture_output=True, text=True)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                summary = json.loads((output / "summary.json").read_text())
                self.assertEqual(summary["cluster_count"], expected_k)
                with (output / "assignments.csv").open() as handle:
                    assignments = list(csv.DictReader(handle))
                self.assertEqual([row["seq_id"] for row in assignments], ["a", "b", "c", "d"])
                with (output / "ordering.csv").open() as handle:
                    self.assertEqual({row["seq_id"] for row in csv.DictReader(handle)}, {"a", "b", "c", "d"})
                with (output / "edges.csv").open() as handle:
                    edges = list(csv.DictReader(handle))
                self.assertEqual(len(edges), 3)
                self.assertEqual(sum(int(row["is_cut"]) for row in edges), expected_k - 1)
            failed = subprocess.run(command + ["--threshold", "1"], capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("--threshold must be in (0, 1)", failed.stderr)

    def test_two_opt_reduces_objective_and_preserves_path(self):
        graph = SparseGraph(
            np.asarray([0, 0, 0, 1, 1, 2]), np.asarray([1, 2, 3, 2, 3, 3]),
            np.asarray([.2, .05, .3, .05, .05, .2]), np.eye(4), ["AAAA"] * 4)
        before = np.arange(4)
        after = _two_opt_sparse(before, graph, .25)
        cost = lambda order: sum(min(graph.exact(int(a), int(b)), .25)
                                 for a, b in zip(order[:-1], order[1:]))
        self.assertEqual(sorted(after.tolist()), before.tolist())
        self.assertLess(cost(after), cost(before))


if __name__ == "__main__":
    unittest.main()
