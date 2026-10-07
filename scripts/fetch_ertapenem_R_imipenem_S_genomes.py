#!/usr/bin/env python3
"""Find matched K. pneumoniae AST records and download their public NCBI assemblies."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import urllib.request
import zipfile
from pathlib import Path

import duckdb


NCBI_DATASETS_DOWNLOAD = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accessions}/download"
ASSEMBLY_ACCESSION = re.compile(r"^GC[AF]_\d+\.\d+$")

PAIR_QUERY = """
WITH kp AS (
    SELECT *
    FROM read_parquet(?)
    WHERE lower(species) LIKE '%klebsiella pneumoniae%'
),
paired_ast AS (
    SELECT DISTINCT
        e.BioSample_ID,
        e.SRA_accession,
        e.assembly_ID,
        e.database AS ertapenem_source,
        e.country,
        e.collection_year,
        e.measurement AS ertapenem_mic,
        e.measurement_sign AS ertapenem_mic_sign,
        e.measurement_units AS ertapenem_mic_units,
        e.ast_standard AS ertapenem_ast_standard,
        i.database AS imipenem_source,
        i.measurement AS imipenem_mic,
        i.measurement_sign AS imipenem_mic_sign,
        i.measurement_units AS imipenem_mic_units,
        i.ast_standard AS imipenem_ast_standard
    FROM kp AS e
    JOIN kp AS i ON e.BioSample_ID = i.BioSample_ID
    WHERE lower(trim(e.antibiotic_name)) = 'ertapenem'
      AND lower(trim(i.antibiotic_name)) = 'imipenem'
      AND lower(trim(e.resistance_phenotype)) = 'resistant'
      AND lower(trim(i.resistance_phenotype)) = 'susceptible'
      AND e.BioSample_ID IS NOT NULL
      AND trim(e.BioSample_ID) <> ''
)
SELECT *
FROM paired_ast
ORDER BY assembly_ID, BioSample_ID
"""


def fetch_pairs(phenotype_path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not phenotype_path.is_file():
        raise FileNotFoundError(
            f"Phenotype parquet not found: {phenotype_path}. "
            "Obtain the public source dataset first; see data/README.md."
        )

    connection = duckdb.connect()
    try:
        result = connection.execute(PAIR_QUERY, [str(phenotype_path)])
        fieldnames = [column[0] for column in result.description]
        rows = [
            {key: "" if value is None else str(value) for key, value in zip(fieldnames, record)}
            for record in result.fetchall()
        ]
    finally:
        connection.close()

    if not rows:
        raise ValueError("No paired ertapenem-resistant/imipenem-susceptible K. pneumoniae records found.")
    return fieldnames, rows


def write_manifest(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def archive_accessions(archive_path: Path) -> set[str]:
    with zipfile.ZipFile(archive_path) as archive:
        return {
            Path(member).parts[2]
            for member in archive.namelist()
            if len(Path(member).parts) > 3
            and Path(member).parts[:2] == ("ncbi_dataset", "data")
            and Path(member).suffix == ".fna"
        }


def download_assemblies(accessions: list[str], output_dir: Path, archive_path: Path) -> None:
    expected = set(accessions)
    if archive_path.exists():
        found = archive_accessions(archive_path)
        missing = expected - found
        if missing:
            raise ValueError(
                f"Cached NCBI package {archive_path} is missing assemblies: {', '.join(sorted(missing))}. "
                "Move the stale package and rerun to download a fresh one."
            )
    else:
        url = NCBI_DATASETS_DOWNLOAD.format(accessions=",".join(accessions))
        url += "?include_annotation_type=GENOME_FASTA"
        request = urllib.request.Request(url, headers={"User-Agent": "ertapenem-kp-amr-profiling/1.0"})
        temporary_archive = archive_path.with_suffix(archive_path.suffix + ".part")
        try:
            with urllib.request.urlopen(request, timeout=300) as response, temporary_archive.open("wb") as output:
                shutil.copyfileobj(response, output)
            found = archive_accessions(temporary_archive)
            missing = expected - found
            if missing:
                raise ValueError(f"NCBI package is missing expected assemblies: {', '.join(sorted(missing))}")
            temporary_archive.replace(archive_path)
        finally:
            temporary_archive.unlink(missing_ok=True)

    base = output_dir.resolve()
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.namelist():
            destination = (output_dir / member).resolve()
            if destination != base and base not in destination.parents:
                raise ValueError(f"Unsafe path in NCBI package: {member}")
        archive.extractall(output_dir)

    extracted = {
        fasta.parent.name
        for fasta in output_dir.glob("ncbi_dataset/data/*/*.fna")
    }
    missing = expected - extracted
    if missing:
        raise ValueError(f"Assemblies were not extracted: {', '.join(sorted(missing))}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter paired AST records and download matching public K. pneumoniae genome assemblies."
    )
    parser.add_argument(
        "--phenotype-data",
        type=Path,
        default=Path("data/raw/phenotype.parquet"),
        help="Local copy of the public phenotype parquet source.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/processed/ertapenem_R_imipenem_S"),
        help="Ignored local output directory for the TSV manifest and NCBI FASTA package.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fieldnames, rows = fetch_pairs(args.phenotype_data)
    accessions = sorted(
        {
            row["assembly_ID"].strip()
            for row in rows
            if ASSEMBLY_ACCESSION.fullmatch(row["assembly_ID"].strip())
        }
    )
    if not accessions:
        raise ValueError("Matched AST records contain no valid NCBI assembly accessions.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / "phenotype_pairs.tsv"
    archive_path = args.output_dir / "ncbi_genome_fastas.zip"
    write_manifest(manifest_path, fieldnames, rows)
    download_assemblies(accessions, args.output_dir, archive_path)

    print(f"Matched AST record pairs: {len(rows)}")
    print(f"Distinct assemblies downloaded: {len(accessions)}")
    print(f"Phenotype manifest: {manifest_path}")
    print(f"Genome FASTAs: {args.output_dir / 'ncbi_dataset' / 'data'}")
    print(f"NCBI package archive: {archive_path}")


if __name__ == "__main__":
    main()
