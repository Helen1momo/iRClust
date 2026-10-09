"""Sparse distance graph and threshold-regularized rearrangement."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio import Align
from numba import njit


LEVELS = ("class", "family", "genus")
_BASE_BITS = {"A": 0, "C": 1, "G": 2, "T": 3}
_MASK64 = (1 << 64) - 1


@dataclass(frozen=True)
class Record:
    sequence_id: str
    sequence: str
    taxonomy: dict[str, str]


@dataclass(frozen=True)
class Result:
    order: np.ndarray
    labels: np.ndarray
    gaps: np.ndarray
    breaks: np.ndarray
    runtime_seconds: float


def read_dataset(directory: Path) -> list[Record]:
    name = directory.name
    with (directory / f"{name}.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {"seq_id", *LEVELS}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{name}: taxonomy CSV must contain {sorted(required)}")
    taxonomy = {
        row["seq_id"]: {level: row[level].strip() for level in LEVELS}
        for row in rows
    }
    if len(taxonomy) != len(rows):
        raise ValueError(f"{name}: duplicate taxonomy IDs")

    records: list[Record] = []
    sequence_id: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        nonlocal sequence_id, chunks
        if sequence_id is None:
            return
        sequence = "".join(chunks).replace(" ", "").upper().replace("U", "T")
        labels = taxonomy.get(sequence_id)
        if not sequence or labels is None or any(not value for value in labels.values()):
            raise ValueError(f"{name}: invalid FASTA/taxonomy row for {sequence_id!r}")
        records.append(Record(sequence_id, sequence, labels))
        sequence_id, chunks = None, []

    with (directory / f"{name}.fasta").open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                sequence_id = line[1:].split()[0]
            elif sequence_id is None:
                raise ValueError(f"{name}: sequence before FASTA header at line {line_number}")
            else:
                chunks.append(line)
    flush()
    ids = [record.sequence_id for record in records]
    if len(ids) != len(set(ids)) or set(ids) != set(taxonomy):
        raise ValueError(f"{name}: FASTA and taxonomy IDs must match one-to-one")
    return records


def fingerprint(records: list[Record]) -> str:
    digest = hashlib.sha256()
    for record in records:
        digest.update(record.sequence_id.encode())
        digest.update(b"\0")
        digest.update(record.sequence.encode("ascii"))
        digest.update(b"\0")
    return digest.hexdigest()


def _mix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & _MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & _MASK64
    return value ^ (value >> 31)


def _sketch(sequences: list[str], k: int, size: int) -> np.ndarray:
    features = np.zeros((len(sequences), size), dtype=np.float32)
    code_mask = (1 << (2 * k)) - 1
    for row, sequence in enumerate(sequences):
        code = valid = 0
        for base in sequence:
            bits = _BASE_BITS.get(base)
            if bits is None:
                code = valid = 0
                continue
            code = ((code << 2) | bits) & code_mask
            valid += 1
            if valid >= k:
                features[row, _mix64(code) % size] += 1
        norm = np.linalg.norm(features[row])
        if norm:
            features[row] /= norm
        else:
            features[row, _mix64(len(sequence)) % size] = 1
    return features


def _candidate_pairs(features: np.ndarray, neighbors: int, threads: int) -> np.ndarray:
    import hnswlib

    index = hnswlib.Index(space="cosine", dim=features.shape[1])
    index.init_index(max_elements=len(features), ef_construction=200, M=16, random_seed=42)
    index.add_items(features, np.arange(len(features)), num_threads=threads)
    index.set_ef(max(256, 2 * (neighbors + 1)))
    labels, _ = index.knn_query(features, k=min(neighbors + 1, len(features)), num_threads=threads)
    rows = np.repeat(np.arange(len(features), dtype=np.int64), labels.shape[1])
    columns = labels.reshape(-1).astype(np.int64)
    keep = rows != columns
    u, v = np.minimum(rows[keep], columns[keep]), np.maximum(rows[keep], columns[keep])
    encoded = np.unique(u * len(features) + v)
    return np.column_stack((encoded // len(features), encoded % len(features)))


def _vsearch_distances(
    sequences: list[str], pairs: np.ndarray, executable: str, threads: int
) -> tuple[np.ndarray, int]:
    by_query: dict[int, list[tuple[int, int]]] = {}
    for pair_index, (query, target) in enumerate(pairs):
        by_query.setdefault(int(query), []).append((pair_index, int(target)))
    queries = sorted(by_query)
    batches = [queries[start : start + 4] for start in range(0, len(queries), 4)]
    values = np.full(len(pairs), np.nan, dtype=np.float32)
    evaluations = 0
    for batch in batches:
        targets = {target for query in batch for _, target in by_query[query]}
        evaluations += len(batch) * len(targets)

    with tempfile.TemporaryDirectory(prefix="rearrangement_vsearch_") as temp_name:
        root = Path(temp_name)

        def run_batch(item: tuple[int, list[int]]) -> tuple[np.ndarray, np.ndarray]:
            batch_number, batch_queries = item
            requested: dict[tuple[int, int], int] = {}
            targets: set[int] = set()
            for query in batch_queries:
                for pair_index, target in by_query[query]:
                    requested[(query, target)] = pair_index
                    targets.add(target)
            query_path = root / f"q{batch_number}.fasta"
            target_path = root / f"t{batch_number}.fasta"
            query_path.write_text(
                "".join(f">{node}\n{sequences[node]}\n" for node in batch_queries),
                encoding="utf-8",
            )
            target_path.write_text(
                "".join(f">{node}\n{sequences[node]}\n" for node in sorted(targets)),
                encoding="utf-8",
            )
            command = [
                executable, "--usearch_global", str(query_path), "--db", str(target_path),
                "--id", "0", "--iddef", "2", "--maxaccepts", "0", "--maxrejects", "0",
                "--minwordmatches", "0", "--minseqlength", "1", "--strand", "plus",
                "--threads", "1", "--userout", "-", "--userfields",
                "query+target+ids+mism+gaps", "--quiet",
            ]
            completed = subprocess.run(command, capture_output=True, text=True)
            if completed.returncode:
                raise RuntimeError(completed.stderr.strip())
            indices, distances = [], []
            for line in completed.stdout.splitlines():
                query, target, matches, mismatches, gaps = map(int, line.split("\t"))
                pair_index = requested.get((query, target))
                if pair_index is None:
                    continue
                denominator = matches + mismatches + gaps
                indices.append(pair_index)
                distances.append(1 - matches / denominator if denominator else 1)
            if len(indices) != len(requested):
                raise RuntimeError("VSEARCH did not return every requested candidate pair")
            return np.asarray(indices), np.asarray(distances, dtype=np.float32)

        with ThreadPoolExecutor(max_workers=threads) as executor:
            for indices, distances in executor.map(run_batch, enumerate(batches)):
                values[indices] = distances
    if np.isnan(values).any():
        raise RuntimeError("VSEARCH candidate refinement contains missing values")
    return values, evaluations


def _retain_nearest(
    n: int, pairs: np.ndarray, distances: np.ndarray, neighbors: int
) -> tuple[np.ndarray, np.ndarray]:
    incident: list[list[tuple[float, int, int]]] = [[] for _ in range(n)]
    encoded = pairs[:, 0] * n + pairs[:, 1]
    for index, (u, v) in enumerate(pairs):
        item = (float(distances[index]), int(v), int(encoded[index]))
        incident[int(u)].append(item)
        incident[int(v)].append((item[0], int(u), item[2]))
    retained: set[int] = set()
    for items in incident:
        items.sort()
        retained.update(item[2] for item in items[:neighbors])
    retained = np.asarray(sorted(retained), dtype=np.int64)
    source_order = np.argsort(encoded)
    selected = source_order[np.searchsorted(encoded[source_order], retained)]
    return np.column_stack((retained // n, retained % n)), distances[selected]


def _aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -4
    aligner.internal_open_gap_score = -20
    aligner.internal_extend_gap_score = -2
    aligner.end_open_gap_score = -2
    aligner.end_extend_gap_score = -1
    aligner.wildcard = "N"
    return aligner


def _exact_distance(aligner: Align.PairwiseAligner, a: str, b: str) -> float:
    alignment = aligner.align(a, b)[0]
    indices = alignment.indices
    both = np.flatnonzero((indices[0] >= 0) & (indices[1] >= 0))
    if not len(both):
        return 1.0
    segment = indices[:, both[0] : both[-1] + 1]
    aligned = (segment[0] >= 0) & (segment[1] >= 0)
    a_bytes = np.frombuffer(a.encode("ascii"), dtype=np.uint8)[segment[0, aligned]]
    b_bytes = np.frombuffer(b.encode("ascii"), dtype=np.uint8)[segment[1, aligned]]
    matches = int(((a_bytes == b_bytes) | (a_bytes == ord("N")) | (b_bytes == ord("N"))).sum())
    denominator = len(a_bytes) + int((~aligned).sum())
    return 1 - matches / denominator if denominator else 1.0


class SparseGraph:
    def __init__(
        self,
        edge_u: np.ndarray,
        edge_v: np.ndarray,
        distances: np.ndarray,
        features: np.ndarray,
        sequences: list[str],
    ) -> None:
        self.n = len(sequences)
        self.edge_u = np.asarray(edge_u, dtype=np.int64)
        self.edge_v = np.asarray(edge_v, dtype=np.int64)
        self.distances = np.asarray(distances, dtype=np.float32)
        self.features = np.asarray(features, dtype=np.float32)
        self.sequences = sequences
        self._aligner = _aligner()
        self._lookup = {
            int(u) * self.n + int(v): float(distance)
            for u, v, distance in zip(self.edge_u, self.edge_v, self.distances)
        }
        self._exact_cache: dict[int, float] = {}
        lists: list[list[tuple[float, int]]] = [[] for _ in range(self.n)]
        for u, v, distance in zip(self.edge_u, self.edge_v, self.distances):
            lists[int(u)].append((float(distance), int(v)))
            lists[int(v)].append((float(distance), int(u)))
        self.neighbor_indices = []
        for items in lists:
            items.sort()
            self.neighbor_indices.append(np.asarray([item[1] for item in items], dtype=np.int64))

    def known(self, i: int, j: int) -> float | None:
        if i == j:
            return 0.0
        u, v = (i, j) if i < j else (j, i)
        return self._lookup.get(u * self.n + v)

    def exact(self, i: int, j: int) -> float:
        value = self.known(i, j)
        if value is not None:
            return value
        u, v = (i, j) if i < j else (j, i)
        key = u * self.n + v
        if key in self._exact_cache:
            return self._exact_cache[key]
        value = _exact_distance(self._aligner, self.sequences[u], self.sequences[v])
        self._exact_cache[key] = value
        return value


def load_or_build_graph(
    records: list[Record], cache_path: Path, *, executable: str, threads: int
) -> SparseGraph:
    parameters = {"sketch_k": 8, "sketch_size": 256, "candidate_k": 128, "exact_k": 128}
    metadata_path = cache_path.with_suffix(".json")
    expected_fingerprint = fingerprint(records)
    if cache_path.is_file() and metadata_path.is_file():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("fingerprint") == expected_fingerprint and all(
            metadata.get(key) == value for key, value in parameters.items()
        ):
            cached = np.load(cache_path)
            print(f"  cache: {cache_path.name} ({len(cached['edge_u']):,} edges)")
            return SparseGraph(
                cached["edge_u"], cached["edge_v"], cached["edge_distance"],
                cached["features"], [record.sequence for record in records],
            )

    if shutil.which(executable) is None:
        raise FileNotFoundError(
            f"A matching graph cache is unavailable. Install VSEARCH and "
            f"requirements-build.txt to build it; executable: {executable!r}"
        )
    started = time.perf_counter()
    sequences = [record.sequence for record in records]
    features = _sketch(sequences, parameters["sketch_k"], parameters["sketch_size"])
    pairs = _candidate_pairs(features, parameters["candidate_k"], threads)
    distances, evaluations = _vsearch_distances(sequences, pairs, executable, threads)
    retained, distances = _retain_nearest(len(records), pairs, distances, parameters["exact_k"])
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path, edge_u=retained[:, 0], edge_v=retained[:, 1],
        edge_distance=distances.astype(np.float32), features=features,
    )
    metadata_path.write_text(
        json.dumps(
            {
                "fingerprint": expected_fingerprint,
                "num_sequences": len(records),
                **parameters,
                "exact_distance": "1 - VSEARCH global identity (iddef=2)",
                "exact_alignment_evaluations": evaluations,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"  cache built: {len(retained):,} edges in {time.perf_counter() - started:.1f}s")
    return SparseGraph(retained[:, 0], retained[:, 1], distances, features, sequences)


@njit(cache=True)
def _path_forest(n: int, edge_u: np.ndarray, edge_v: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    parent = np.arange(n, dtype=np.int64)
    rank = np.zeros(n, dtype=np.int8)
    degree = np.zeros(n, dtype=np.int8)
    adjacency = np.full((n, 2), -1, dtype=np.int64)

    def find(value: int) -> int:
        root = value
        while parent[root] != root:
            root = parent[root]
        while parent[value] != value:
            previous = parent[value]
            parent[value] = root
            value = previous
        return root

    for u, v in zip(edge_u, edge_v):
        if degree[u] == 2 or degree[v] == 2:
            continue
        root_u, root_v = find(u), find(v)
        if root_u == root_v:
            continue
        adjacency[u, degree[u]], adjacency[v, degree[v]] = v, u
        degree[u] += 1
        degree[v] += 1
        if rank[root_u] < rank[root_v]:
            parent[root_u] = root_v
        elif rank[root_u] > rank[root_v]:
            parent[root_v] = root_u
        else:
            parent[root_v] = root_u
            rank[root_u] += 1
    return adjacency, degree


def _fragments(adjacency: np.ndarray, degree: np.ndarray) -> list[np.ndarray]:
    visited = np.zeros(len(degree), dtype=bool)
    output = []
    for start in np.flatnonzero(degree < 2):
        if visited[start]:
            continue
        path, previous, current = [], -1, int(start)
        while current >= 0 and not visited[current]:
            visited[current] = True
            path.append(current)
            next_node = next((int(node) for node in adjacency[current] if node >= 0 and node != previous), -1)
            previous, current = current, next_node
        output.append(np.asarray(path, dtype=np.int64))
    if not visited.all():
        raise RuntimeError("Path forest contains a cycle")
    return output


def _join_sparse(fragments: list[np.ndarray], graph: SparseGraph) -> np.ndarray:
    remaining = sorted(fragments, key=lambda item: (-len(item), int(item[0])))
    path = remaining.pop(0).tolist()
    while remaining:
        endpoints: dict[int, tuple[int, bool]] = {}
        for index, fragment in enumerate(remaining):
            endpoints[int(fragment[0])] = (index, False)
            endpoints[int(fragment[-1])] = (index, True)
        current = path[-1]
        selected = next((endpoints[int(node)] for node in graph.neighbor_indices[current] if int(node) in endpoints), None)
        if selected is None:
            candidates = np.fromiter(endpoints, dtype=np.int64)
            similarities = graph.features[candidates] @ graph.features[current]
            endpoint = int(candidates[np.lexsort((candidates, -similarities))[0]])
            graph.exact(current, endpoint)
            selected = endpoints[endpoint]
        index, reverse = selected
        fragment = remaining.pop(index)
        path.extend(map(int, fragment[::-1] if reverse else fragment))
    return np.asarray(path, dtype=np.int64)


def _two_opt_sparse(order: np.ndarray, graph: SparseGraph, delta: float, passes: int = 25) -> np.ndarray:
    order = order.copy()
    position = np.empty(len(order), dtype=np.int64)
    position[order] = np.arange(len(order))
    for _ in range(passes):
        moves = 0
        for i in range(len(order) - 2):
            a, b = int(order[i]), int(order[i + 1])
            old_left = graph.exact(a, b)
            for c in graph.neighbor_indices[a][:32]:
                c = int(c)
                j = int(position[c])
                if j <= i + 1 or j >= len(order) - 1:
                    continue
                d = int(order[j + 1])
                old_right, new_left, new_right = graph.exact(c, d), graph.known(a, c), graph.known(b, d)
                if None in (new_left, new_right):
                    continue
                old_cost = min(old_left, delta) + min(float(old_right), delta)
                new_cost = min(float(new_left), delta) + min(float(new_right), delta)
                if new_cost + 1e-12 < old_cost:
                    order[i + 1 : j + 1] = order[i + 1 : j + 1][::-1]
                    position[order[i + 1 : j + 1]] = np.arange(i + 1, j + 1)
                    moves += 1
                    break
        if not moves:
            break
    return order


def rearrange_sparse(graph: SparseGraph, threshold: float) -> Result:
    if not 0 < threshold < 1:
        raise ValueError("threshold must be in (0, 1)")
    started = time.perf_counter()
    delta = 1 - threshold
    keep = graph.distances < delta
    u, v, distance = graph.edge_u[keep], graph.edge_v[keep], graph.distances[keep]
    encoded = u * graph.n + v
    edge_order = np.lexsort((encoded, distance))
    adjacency, degree = _path_forest(graph.n, u[edge_order], v[edge_order])
    order = _join_sparse(_fragments(adjacency, degree), graph)
    order = _two_opt_sparse(order, graph, delta)
    gaps = np.asarray([graph.exact(int(a), int(b)) for a, b in zip(order[:-1], order[1:])])
    breaks = gaps >= delta
    ordered_labels = np.r_[0, np.cumsum(breaks)].astype(np.int64)
    labels = np.empty(graph.n, dtype=np.int64)
    labels[order] = ordered_labels
    if len(np.unique(order)) != graph.n:
        raise RuntimeError("Rearrangement did not produce a complete path")
    return Result(order, labels, gaps, breaks, time.perf_counter() - started)


def rearrange_dense(distance: np.ndarray, threshold: float) -> Result:
    """Small full-matrix solver used only by the 30-sequence simulation."""

    started = time.perf_counter()
    n, delta = len(distance), 1 - threshold
    u, v = np.triu_indices(n, 1)
    keep = distance[u, v] < delta
    u, v = u[keep], v[keep]
    edge_order = np.lexsort((u * n + v, distance[u, v]))
    adjacency, degree = _path_forest(n, u[edge_order], v[edge_order])
    remaining = sorted(_fragments(adjacency, degree), key=lambda item: (-len(item), int(item[0])))
    path = remaining.pop(0).tolist()
    while remaining:
        left = np.asarray([item[0] for item in remaining])
        right = np.asarray([item[-1] for item in remaining])
        left_i, right_i = int(np.argmin(distance[path[-1], left])), int(np.argmin(distance[path[-1], right]))
        reverse = distance[path[-1], right[right_i]] < distance[path[-1], left[left_i]]
        index = right_i if reverse else left_i
        fragment = remaining.pop(index)
        path.extend(map(int, fragment[::-1] if reverse else fragment))
    order = np.asarray(path, dtype=np.int64)
    improved = True
    while improved:
        improved = False
        for i in range(n - 2):
            for j in range(i + 2, n - 1):
                a, b, c, d = order[i], order[i + 1], order[j], order[j + 1]
                old = min(distance[a, b], delta) + min(distance[c, d], delta)
                new = min(distance[a, c], delta) + min(distance[b, d], delta)
                if new + 1e-12 < old:
                    order[i + 1 : j + 1] = order[i + 1 : j + 1][::-1]
                    improved = True
                    break
            if improved:
                break
    gaps = distance[order[:-1], order[1:]]
    breaks = gaps >= delta
    labels = np.empty(n, dtype=np.int64)
    labels[order] = np.r_[0, np.cumsum(breaks)]
    return Result(order, labels, gaps, breaks, time.perf_counter() - started)
