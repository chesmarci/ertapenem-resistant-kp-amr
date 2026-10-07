#!/usr/bin/env python3
"""Plot detected carbapenemase gene families across all genomes in a merged AMR TSV."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from sklearn.preprocessing import MultiLabelBinarizer


CARBAPENEMASE_FAMILIES = (
    ("KPC", re.compile(r"^blaKPC(?:[-_]?\d+)?$", re.IGNORECASE)),
    ("NDM", re.compile(r"^blaNDM(?:[-_]?\d+)?$", re.IGNORECASE)),
    (
        "OXA-48-like",
        re.compile(r"^blaOXA[-_]?(?:48|181|204|232|244|245|247|436|484)(?:[-_][A-Za-z0-9]+)?$", re.IGNORECASE),
    ),
    ("VIM", re.compile(r"^blaVIM(?:[-_]?\d+)?$", re.IGNORECASE)),
    ("IMP", re.compile(r"^blaIMP(?:[-_]?\d+)?$", re.IGNORECASE)),
)


def load_carbapenemases(path: Path) -> tuple[list[str], dict[tuple[str, str], set[str]]]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "Element symbol", "Type"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")

        genomes: set[str] = set()
        calls: dict[tuple[str, str], set[str]] = {}
        for row in reader:
            genome = (row.get("sample") or "").strip()
            symbol = (row.get("Element symbol") or "").strip()
            if not genome:
                continue
            genomes.add(genome)
            if (row.get("Type") or "").strip().upper() != "AMR":
                continue
            for family, pattern in CARBAPENEMASE_FAMILIES:
                if pattern.fullmatch(symbol):
                    calls.setdefault((family, genome), set()).add(symbol)
                    break

    if not genomes:
        raise ValueError(f"No genomes found in {path}.")
    return sorted(genomes), calls


def render_heatmap(genomes: list[str], calls: dict[tuple[str, str], set[str]], output: Path) -> None:
    family_names = [family for family, _ in CARBAPENEMASE_FAMILIES]
    genome_families = [
        [family for family in family_names if calls.get((family, genome))]
        for genome in genomes
    ]
    encoder = MultiLabelBinarizer(classes=family_names)
    presence_matrix = encoder.fit_transform(genome_families)

    genome_labels = ["_".join(genome.split("_")[:2]) for genome in genomes]
    figure_width = max(8, 2.5 + 1.8 * len(genomes))
    figure, axis = plt.subplots(figsize=(figure_width, 5.2))
    figure.subplots_adjust(left=0.18, right=0.98, top=0.78, bottom=0.25)
    colors = ListedColormap(["#edf1f5", "#c2413b"])
    axis.imshow(presence_matrix.T, cmap=colors, vmin=0, vmax=1, aspect="auto")

    axis.set_xticks(range(len(genomes)), labels=genome_labels, rotation=25, ha="right")
    axis.set_yticks(range(len(family_names)), labels=family_names)
    axis.set_xlabel("Genome")
    axis.set_ylabel("Carbapenemase family")
    axis.set_xticks([index - 0.5 for index in range(1, len(genomes))], minor=True)
    axis.set_yticks([index - 0.5 for index in range(1, len(family_names))], minor=True)
    axis.grid(which="minor", color="white", linewidth=3)
    axis.tick_params(which="minor", bottom=False, left=False)

    for family_index, family in enumerate(family_names):
        for genome_index, genome in enumerate(genomes):
            alleles = sorted(calls.get((family, genome), set()))
            annotation = ", ".join(alleles) if alleles else "not detected"
            text_color = "white" if alleles else "#738294"
            axis.text(
                genome_index,
                family_index,
                annotation,
                ha="center",
                va="center",
                color=text_color,
                fontsize=9,
                weight="bold" if alleles else "normal",
                wrap=True,
            )

    axis.legend(
        handles=[
            Patch(facecolor="#c2413b", label="Detected"),
            Patch(facecolor="#edf1f5", edgecolor="#d2dae3", label="Not detected"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.27),
        ncol=2,
        frameon=False,
    )
    figure.suptitle(
        "Carbapenemase gene heatmap across genomes",
        x=0.04,
        y=0.97,
        ha="left",
        fontsize=16,
        weight="bold",
    )
    figure.text(
        0.04,
        0.90,
        f"AMRFinderPlus calls | {len(genomes)} genomes | "
        f"{int(presence_matrix.sum())} detected family/genome combinations",
        ha="left",
        fontsize=10,
        color="#526579",
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)
    svg_lines = output.read_text(encoding="utf-8").splitlines()
    output.write_text("\n".join(line.rstrip() for line in svg_lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot a carbapenemase-family heatmap for every genome in a TSV.")
    parser.add_argument("--input", type=Path, default=Path("results/amrfinder/combined_amr.tsv"))
    parser.add_argument("--output", type=Path, default=Path("figures/carbapenemase_heatmap.svg"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    genomes, calls = load_carbapenemases(args.input)
    render_heatmap(genomes, calls, args.output)
    print(f"Carbapenemase heatmap for {len(genomes)} genomes written to {args.output}")


if __name__ == "__main__":
    main()
