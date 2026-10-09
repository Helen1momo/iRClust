# Rearrangement Clustering

Standalone nucleotide sequence clustering with six bundled datasets.

## Install

Use Python 3.10–3.12. In the repository directory, run:

```bash
python -m pip install -r requirements.txt
```

## Run

```bash
python run_clustering.py
```

This runs `16s_01` at similarity threshold `0.95`. The bundled datasets include precomputed graphs, so VSEARCH is not required.

Choose a dataset and threshold:

```bash
python run_clustering.py --dataset 18s_02 --threshold 0.97
```

Available datasets: `16s_01`, `16s_02`, `18s_01`, `18s_02`, `coi_01`, `coi_02`.

Specify an output directory or view all options:

```bash
python run_clustering.py --dataset coi_01 --threshold 0.90 --output-dir results/my_run
python run_clustering.py --help
```

## Output

Results are saved to `results/<dataset>/theta_<threshold>/` by default:

- `assignments.csv`: sequence IDs and cluster IDs in input order.
- `ordering.csv`: sequences in rearranged path order.
- `edges.csv`: adjacent sequence distances and cut indicators.
- `summary.json`: cluster count, threshold, and runtime.

Positions and cluster IDs start at zero. Repeating a run overwrites its output files; use `--output-dir` to keep separate runs.

## Custom data

Install VSEARCH and make it available on `PATH`, then install the additional dependencies:

```bash
python -m pip install -r requirements-build.txt
```

Prepare a directory `my_sequences/` containing `my_sequences.fasta` and `my_sequences.csv`. The CSV must contain `seq_id,class,family,genus`, with IDs matching the FASTA headers. Use `NA` for unavailable taxonomy labels.

```bash
python run_clustering.py --data-dir my_sequences --threshold 0.95
```

The first run builds a graph cache; subsequent runs reuse it. Use `--vsearch /path/to/vsearch` if VSEARCH is not on `PATH`.
