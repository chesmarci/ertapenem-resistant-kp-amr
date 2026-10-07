#!/usr/bin/env python3
"""Run AMRFinderPlus on every assembly FASTA in a cohort directory."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run AMRFinderPlus --plus for every NCBI assembly FASTA in a cohort directory."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/processed/ertapenem_R_imipenem_S/ncbi_dataset/data"),
        help="Directory containing one NCBI Datasets accession folder per assembly.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/ertapenem_R_imipenem_S/amrfinder"),
        help="Directory for per-genome AMRFinderPlus TSV outputs.",
    )
    parser.add_argument(
        "--database",
        type=Path,
        default=Path("data/local_amrfinder_db/2026-08-07.1"),
        help="AMRFinderPlus database version directory.",
    )
    parser.add_argument("--threads", type=int, default=4, help="Threads used per genome.")
    parser.add_argument("--amrfinder", default="amrfinder", help="AMRFinderPlus executable path.")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run genomes that already have a valid output table.",
    )
    return parser.parse_args()


def has_valid_output(path: Path) -> bool:
    if not path.is_file():
        return False
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return bool(reader.fieldnames and {"Element symbol", "Type"}.issubset(reader.fieldnames))


def main() -> None:
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be at least 1.")
    amrfinder = shutil.which(args.amrfinder)
    if amrfinder is None:
        raise FileNotFoundError(
            f"AMRFinderPlus executable not found: {args.amrfinder!r}. "
            "Add it to PATH or pass its full path with --amrfinder "
            "(for example, <conda-environment>/bin/amrfinder)."
        )
    if not args.database.is_dir():
        raise FileNotFoundError(f"AMRFinderPlus database not found: {args.database}")

    assemblies = sorted(
        (fasta.parent.name, fasta)
        for prefix in ("GCA_", "GCF_")
        for fasta in args.input_dir.glob(f"{prefix}*/*.fna")
    )
    if not assemblies:
        raise FileNotFoundError(f"No GCA assembly FASTAs found under {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    skipped = 0
    for index, (accession, fasta) in enumerate(assemblies, start=1):
        output = args.output_dir / f"{accession}_genomic_amr.tsv"
        if not args.force and has_valid_output(output):
            skipped += 1
            print(f"[{index}/{len(assemblies)}] {accession}: output exists; skipping", flush=True)
            continue

        temporary_output = output.with_suffix(output.suffix + ".tmp")
        command = [
            amrfinder,
            "--database",
            str(args.database),
            "--nucleotide",
            str(fasta),
            "--organism",
            "Klebsiella_pneumoniae",
            "--plus",
            "--threads",
            str(args.threads),
            "--output",
            str(temporary_output),
        ]
        print(f"[{index}/{len(assemblies)}] {accession}", flush=True)
        try:
            subprocess.run(command, check=True)
            if not has_valid_output(temporary_output):
                raise ValueError(f"AMRFinderPlus created an invalid output table for {accession}.")
            temporary_output.replace(output)
        finally:
            temporary_output.unlink(missing_ok=True)

    version = subprocess.run(
        [amrfinder, "--database", str(args.database), "--database_version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    metadata = args.output_dir / "run_metadata.txt"
    metadata.write_text(
        "\n".join(
            [
                f"Completed (UTC): {datetime.now(timezone.utc).isoformat()}",
                f"Assemblies processed: {len(assemblies)}",
                f"Existing outputs skipped: {skipped}",
                f"AMRFinderPlus database: {version}",
                "Options: --plus --organism Klebsiella_pneumoniae",
                f"Threads per assembly: {args.threads}",
                f"Forced re-run: {args.force}",
                f"Input directory: {args.input_dir}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(f"Completed {len(assemblies)} genomes.")
    print(f"AMRFinderPlus tables: {args.output_dir}")
    print(f"Run metadata: {metadata}")


if __name__ == "__main__":
    main()
