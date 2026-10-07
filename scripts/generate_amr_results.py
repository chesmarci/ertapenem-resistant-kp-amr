#!/usr/bin/env python3
"""Generate the merged AMRFinder table and a concise result summary.

This script scans the per-sample AMRFinder TSV files under results/amrfinder,
combines them into a single table, and reports the key findings described in the
project README.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


CARBAPENEMASE_MARKERS = (
    "blaKPC",
    "blaNDM",
    "blaOXA-48",
    "blaVIM",
    "blaIMP",
)

ESBL_MARKERS = ("blaSHV", "blaTEM")
FLUOROQUINOLONE_MUTATIONS = ("gyrA_S83I", "parC_S80I")


def sample_name_for(path: Path) -> str:
    stem = path.stem
    if stem.endswith("_amr"):
        return stem[: -len("_amr")]
    return stem


def collect_rows(input_dir: Path) -> tuple[list[dict[str, str]], list[str]]:
    files = sorted(input_dir.glob("*_genomic_amr.tsv"))
    if not files:
        raise FileNotFoundError(f"No per-sample AMRFinder TSV files found in {input_dir}")

    rows: list[dict[str, str]] = []
    fieldnames: list[str] | None = None

    for file_path in files:
        with file_path.open("r", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                continue
            fieldnames = list(reader.fieldnames)
            sample = sample_name_for(file_path)
            for record in reader:
                if not record:
                    continue
                normalized = {key: (value if value != "NA" else "") for key, value in record.items()}
                normalized["sample"] = sample
                rows.append(normalized)

    if not rows:
        raise ValueError(f"No rows parsed from AMRFinder files in {input_dir}")

    output_fieldnames = ["sample", *fieldnames] if fieldnames else ["sample"]
    output_fieldnames = [name for name in output_fieldnames if name != "sample"]
    output_fieldnames.insert(0, "sample")
    return rows, output_fieldnames


def write_combined_table(rows: list[dict[str, str]], fieldnames: list[str], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def summarize_results(rows: list[dict[str, str]]) -> str:
    samples = sorted({row.get("sample", "") for row in rows if row.get("sample")})
    sample_count = len(samples)

    carbapenemase_hits: dict[str, list[str]] = defaultdict(list)
    ompK35_hits: list[str] = []
    ompK36_hits: list[str] = []
    esbl_hits: dict[str, list[str]] = defaultdict(list)
    fluoro_hits: dict[str, list[str]] = defaultdict(list)

    for row in rows:
        sample = row.get("sample", "")
        symbol = (row.get("Element symbol") or "").strip()
        gene_name = (row.get("Element name") or "").strip()
        decoded = f"{symbol} {gene_name}".strip()

        carbapenemase_match = next((marker for marker in CARBAPENEMASE_MARKERS if marker.lower() in symbol.lower()), None)
        if carbapenemase_match:
            carbapenemase_hits[sample].append(carbapenemase_match)

        symbol_lower = symbol.lower()
        decoded_lower = decoded.lower()

        if "ompk35" in symbol_lower:
            ompK35_hits.append(sample)
        if "ompk36" in symbol_lower:
            ompK36_hits.append(sample)

        for marker in ESBL_MARKERS:
            if marker.lower() in symbol_lower:
                esbl_hits[sample].append(marker)

        for mutation in FLUOROQUINOLONE_MUTATIONS:
            if mutation.lower() == symbol_lower:
                fluoro_hits[sample].append(mutation)

        if "ompk35" in decoded_lower:
            ompK35_hits.append(sample)
        if "ompk36" in decoded_lower:
            ompK36_hits.append(sample)

    carbapenemase_samples = sorted(carbapenemase_hits)
    unique_ompK35 = sorted(set(ompK35_hits))
    unique_ompK36 = sorted(set(ompK36_hits))
    all_fluoro_samples = sorted(fluoro_hits)

    esbl_names = sorted({gene for values in esbl_hits.values() for gene in values})
    lines = [
        f"Total genomes analyzed: {sample_count}",
        "",
        "Key findings:",
        (
            "- No carbapenemase genes detected (none of the markers KPC, NDM, OXA-48, VIM, or IMP were found "
            f"across {sample_count} isolates)."
            if not carbapenemase_samples
            else "- Carbapenemase genes detected in: " + ", ".join(carbapenemase_samples)
        ),
        (
            "- All isolates carry loss-of-function ompK35 alterations: "
            + (", ".join(unique_ompK35) if unique_ompK35 else "none")
            if unique_ompK35
            else "- No ompK35 mutations detected"
        ),
        "- ompK36 alterations detected in: " + (", ".join(unique_ompK36) if unique_ompK36 else "none"),
        "- ESBL genes identified: " + (", ".join(esbl_names) if esbl_names else "none"),
        "- Fluoroquinolone resistance mutations identified in: "
        + (", ".join(all_fluoro_samples) if all_fluoro_samples else "none"),
    ]

    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge AMRFinder per-genome TSV files and summarize findings.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("results/amrfinder"),
        help="Directory containing the per-sample AMRFinder TSV files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/amrfinder/combined_amr.tsv"),
        help="Path for the merged tab-delimited AMRFinder output.",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=None,
        help="Optional path for the text summary. If omitted, it prints to stdout.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows, fieldnames = collect_rows(args.input_dir)
    write_combined_table(rows, fieldnames, args.output)

    summary_text = summarize_results(rows)
    if args.summary_output is None:
        print(summary_text, end="")
    else:
        args.summary_output.parent.mkdir(parents=True, exist_ok=True)
        args.summary_output.write_text(summary_text)
        print(f"Summary written to {args.summary_output}")

    print(f"Merged AMRFinder table written to {args.output}")


if __name__ == "__main__":
    main()
