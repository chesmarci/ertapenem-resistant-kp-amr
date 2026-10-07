#!/usr/bin/env python3
"""Compare two genome AMR profiles and render an SVG presence/absence heatmap."""

from __future__ import annotations

import argparse
import csv
import html
from collections import defaultdict
from pathlib import Path


def resolve_sample(samples: set[str], accession: str) -> str:
    matches = sorted(sample for sample in samples if sample == accession or sample.startswith(f"{accession}_"))
    if not matches:
        raise ValueError(f"No sample matching genome accession {accession!r} was found.")
    if len(matches) > 1:
        raise ValueError(f"Genome accession {accession!r} matches multiple samples: {', '.join(matches)}")
    return matches[0]


def load_amr_features(path: Path) -> tuple[set[str], dict[str, set[str]], dict[str, str]]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "Element symbol", "Type", "Class"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")

        samples: set[str] = set()
        features: dict[str, set[str]] = defaultdict(set)
        classes: dict[str, str] = {}
        for row in reader:
            sample = (row.get("sample") or "").strip()
            symbol = (row.get("Element symbol") or "").strip()
            if not sample or not symbol or (row.get("Type") or "").strip().upper() != "AMR":
                continue
            samples.add(sample)
            features[sample].add(symbol)
            classes[symbol] = (row.get("Class") or "Unclassified").strip() or "Unclassified"

    if not samples:
        raise ValueError(f"No AMR features were found in {path}.")
    return samples, features, classes


def render_svg(
    samples: tuple[str, str],
    features: dict[str, set[str]],
    classes: dict[str, str],
    output_path: Path,
) -> None:
    first, second = samples
    feature_rows = sorted(
        features[first] | features[second],
        key=lambda symbol: (classes.get(symbol, "Unclassified").casefold(), symbol.casefold()),
    )

    margin_left = 270
    margin_top = 112
    cell_width = 190
    row_height = 27
    width = margin_left + 2 * cell_width + 36
    height = margin_top + row_height * len(feature_rows) + 72
    first_label = "_".join(first.split("_")[:2])
    second_label = "_".join(second.split("_")[:2])
    shared = len(features[first] & features[second])
    first_only = len(features[first] - features[second])
    second_only = len(features[second] - features[first])

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#172b4d}.title{font-size:20px;font-weight:700}",
        ".subtitle{font-size:12px;fill:#526579}.header{font-size:13px;font-weight:700}",
        ".label{font-size:12px}.class{font-size:10px;fill:#526579}.legend{font-size:11px}",
        ".cell-label{font-size:10px}.present{fill:#ffffff}.absent{fill:#778899}",
        "</style>",
        f'<text x="24" y="32" class="title">AMR profile comparison: {html.escape(first_label)} vs {html.escape(second_label)}</text>',
        '<text x="24" y="55" class="subtitle">AMR-typed AMRFinderPlus elements; cells indicate detected element symbols</text>',
        f'<text x="24" y="78" class="subtitle">Shared: {shared}  |  {html.escape(first_label)} only: {first_only}  |  {html.escape(second_label)} only: {second_only}</text>',
        f'<text x="{margin_left + cell_width / 2}" y="101" text-anchor="middle" class="header">{html.escape(first_label)}</text>',
        f'<text x="{margin_left + cell_width + cell_width / 2}" y="101" text-anchor="middle" class="header">{html.escape(second_label)}</text>',
    ]

    for index, symbol in enumerate(feature_rows):
        y = margin_top + index * row_height
        feature_class = classes.get(symbol, "Unclassified")
        lines.append(
            f'<text x="24" y="{y + 16}" class="label">{html.escape(symbol)}</text>'
            f'<text x="{margin_left - 12}" y="{y + 16}" text-anchor="end" class="class">{html.escape(feature_class)}</text>'
        )
        for column, sample in enumerate(samples):
            x = margin_left + column * cell_width + 16
            present = symbol in features[sample]
            fill = "#168a72" if present else "#edf1f5"
            stroke = "#11705d" if present else "#d2dae3"
            label = "present" if present else "not detected"
            lines.append(
                f'<rect x="{x}" y="{y}" width="158" height="20" rx="4" fill="{fill}" stroke="{stroke}"/>'
                f'<text x="{x + 79}" y="{y + 14}" text-anchor="middle" '
                f'class="cell-label {"present" if present else "absent"}">{label}</text>'
            )

    legend_y = margin_top + len(feature_rows) * row_height + 26
    lines.extend(
        [
            f'<rect x="24" y="{legend_y - 12}" width="14" height="14" rx="3" fill="#168a72"/>',
            f'<text x="44" y="{legend_y}" class="legend">Detected</text>',
            f'<rect x="120" y="{legend_y - 12}" width="14" height="14" rx="3" fill="#edf1f5" stroke="#d2dae3"/>',
            f'<text x="140" y="{legend_y}" class="legend">Not detected</text>',
            "</svg>",
        ]
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare two genomes' AMRFinder AMR profiles as an SVG heatmap.")
    parser.add_argument("--input", type=Path, default=Path("results/amrfinder/combined_amr.tsv"))
    parser.add_argument("--genome-a", default="GCA_000633755.1")
    parser.add_argument("--genome-b", default="GCA_000694015.1")
    parser.add_argument("--output", type=Path, default=Path("figures/genome_amr_comparison.svg"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    samples, features, classes = load_amr_features(args.input)
    first = resolve_sample(samples, args.genome_a)
    second = resolve_sample(samples, args.genome_b)
    if first == second:
        raise ValueError("Choose two different genomes.")
    render_svg((first, second), features, classes, args.output)
    print(f"AMR comparison visualization written to {args.output}")


if __name__ == "__main__":
    main()
