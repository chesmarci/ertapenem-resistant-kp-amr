#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

CONDA_EXE="${CONDA_EXE:-$(command -v conda || true)}"
PYTHON="${PYTHON:-python3}"
AMRFINDER_ENV="${AMRFINDER_ENV:-kp-amr-processing}"
DATABASE_VERSION="${AMRFINDER_DATABASE_VERSION:-2026-08-07.1}"
THREADS="${THREADS:-4}"
PHENOTYPE_DATA="${PHENOTYPE_DATA:-data/raw/phenotype.parquet}"
GENOME_DIR="data/processed/ertapenem_R_imipenem_S/ncbi_dataset/data"
DATABASE_ROOT="data/local_amrfinder_db"
DATABASE_DIR="$DATABASE_ROOT/$DATABASE_VERSION"
AMRFINDER_DIR="results/ertapenem_R_imipenem_S/amrfinder"
MERGED_TABLE="results/ertapenem_R_imipenem_S/combined_amr.tsv"
SUMMARY="results/ertapenem_R_imipenem_S/summary_report.txt"
HEATMAP="figures/carbapenemase_heatmap_ertapenem_R_imipenem_S.svg"

if [[ -z "$CONDA_EXE" ]]; then
    echo "Conda was not found. Install Miniconda/Anaconda or set CONDA_EXE." >&2
    exit 1
fi

if ! "$CONDA_EXE" env list --json | "$PYTHON" -c \
    'import json,sys; raise SystemExit(0 if any(p.rstrip("/").endswith("/envs/" + sys.argv[1]) for p in json.load(sys.stdin)["envs"]) else 1)' \
    "$AMRFINDER_ENV"; then
    echo "Conda environment '$AMRFINDER_ENV' was not found." >&2
    echo "Create it with: conda create -n $AMRFINDER_ENV -c conda-forge -c bioconda python=3.12 ncbi-amrfinderplus=4.2.7" >&2
    exit 1
fi

if ! "$CONDA_EXE" run -n "$AMRFINDER_ENV" amrfinder --version >/dev/null; then
    echo "AMRFinderPlus is unavailable in '$AMRFINDER_ENV'." >&2
    exit 1
fi

if ! "$CONDA_EXE" run -n "$AMRFINDER_ENV" python -c \
    'import duckdb, matplotlib, sklearn' >/dev/null 2>&1; then
    echo "Python dependencies are missing from '$AMRFINDER_ENV'." >&2
    echo "Install them with: conda run -n $AMRFINDER_ENV python -m pip install -r requirements.txt" >&2
    exit 1
fi

if [[ ! -d "$GENOME_DIR" ]] || ! find "$GENOME_DIR" -type f -name '*.fna' -print -quit | grep -q .; then
    if [[ ! -f "$PHENOTYPE_DATA" ]]; then
        echo "No assembly FASTAs found under '$GENOME_DIR' and phenotype data is absent at '$PHENOTYPE_DATA'." >&2
        echo "Obtain the public phenotype parquet or set PHENOTYPE_DATA to its location." >&2
        exit 1
    fi
    "$CONDA_EXE" run -n "$AMRFINDER_ENV" python scripts/fetch_ertapenem_R_imipenem_S_genomes.py \
        --phenotype-data "$PHENOTYPE_DATA"
fi

"$CONDA_EXE" run -n "$AMRFINDER_ENV" python -c '
import csv
import sys
from pathlib import Path

genome_dir = Path(sys.argv[1])
manifest = Path(sys.argv[2])
with manifest.open(newline="") as handle:
    expected = {row["assembly_ID"].strip() for row in csv.DictReader(handle, delimiter="\t") if row["assembly_ID"].strip()}
actual = {path.parent.name for path in genome_dir.glob("GC[A,F]_*/*.fna")}
missing = sorted(expected - actual)
if missing:
    raise SystemExit("Missing FASTA assemblies for accessions: " + ", ".join(missing))
if not actual:
    raise SystemExit("No genome FASTAs found in " + str(genome_dir))
print(f"Validated {len(actual)} downloaded assemblies against the phenotype manifest.")
' "$GENOME_DIR" "data/processed/ertapenem_R_imipenem_S/phenotype_pairs.tsv"

if [[ ! -d "$DATABASE_DIR" ]]; then
    echo "Preparing AMRFinderPlus database $DATABASE_VERSION..."
    "$CONDA_EXE" run -n "$AMRFINDER_ENV" amrfinder_update --database "$DATABASE_ROOT"
fi
if [[ ! -d "$DATABASE_DIR" ]]; then
    echo "Expected database version $DATABASE_VERSION was not installed at '$DATABASE_DIR'." >&2
    echo "Set AMRFINDER_DATABASE_VERSION to an installed version or install the requested database." >&2
    exit 1
fi

echo "Running AMRFinderPlus on local assemblies (THREADS=$THREADS)."
"$CONDA_EXE" run -n "$AMRFINDER_ENV" python scripts/run_amrfinder_cohort.py \
    --input-dir "$GENOME_DIR" \
    --output-dir "$AMRFINDER_DIR" \
    --database "$DATABASE_DIR" \
    --threads "$THREADS"

"$CONDA_EXE" run -n "$AMRFINDER_ENV" python scripts/generate_amr_results.py \
    --input-dir "$AMRFINDER_DIR" \
    --output "$MERGED_TABLE" \
    --summary-output "$SUMMARY"

"$CONDA_EXE" run -n "$AMRFINDER_ENV" python scripts/carbapenemase_heatmap.py \
    --input "$MERGED_TABLE" \
    --output "$HEATMAP"

echo "Pipeline complete."
echo "Merged AMR results: $MERGED_TABLE"
echo "Summary: $SUMMARY"
echo "Carbapenemase heatmap: $HEATMAP"
