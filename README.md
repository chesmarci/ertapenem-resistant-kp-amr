# Ertapenem-resistant *Klebsiella pneumoniae*: WGS-based AMR profiling

## Overview

This repository contains scripts and visualizations for whole-genome sequencing
(WGS)-based antimicrobial resistance (AMR) profiling of ertapenem-resistant
*Klebsiella pneumoniae*.

Public phenotype and genome data are referenced but are not redistributed in
this repository. See [Data sources](data/README.md) for dataset provenance and
accession links.

## Data source

The local phenotype export is a public composite table whose provenance
includes NCBI antibiogram, NDARO, PATRIC, CDC, and CABBAGE PubMed records.
Genome assemblies are retrieved from NCBI Datasets. Dataset files, assembly
downloads, and intermediate sequence-analysis files are excluded from Git;
source details and accession links are listed in [data/README.md](data/README.md).

The current merged AMRFinder results cover three genomes. The source workspace
mentions a larger pilot subset, but only these three genomes have per-genome
AMRFinder outputs in the merged results used for the figures.

## Methods

1. Obtain the public CABBAGE phenotype data and document its release, download
   date, and filtering steps.
2. Retrieve the corresponding genome assemblies from NCBI Datasets.
3. Run AMRFinderPlus v4.2.7 with `--plus` for each genome.
4. Place per-genome AMRFinder TSV outputs in `results/amrfinder/` and merge
   them with `scripts/generate_amr_results.py`.
5. Generate the heatmaps with the scripts below.

The source snapshot available locally yielded 21 distinct assembly accessions
with exact paired ertapenem-resistant/imipenem-susceptible calls; this is a
subset defined by phenotype annotations and assembly availability, not a claim
that all records across the contributing databases have been included.

To expand the WGS cohort to public *K. pneumoniae* isolates with paired
ertapenem-resistant and imipenem-susceptible phenotypes, first obtain the local
public composite phenotype parquet, then run:

```bash
python3 scripts/fetch_ertapenem_R_imipenem_S_genomes.py
```

This creates a local matched-phenotype manifest and downloads the matching
NCBI assembly FASTAs. Details and limitations are in [data/README.md](data/README.md).

## Tools

- NCBI Datasets CLI
- AMRFinderPlus v4.2.7 (database 2026-08-07.1)
- Python 3
- AMRFinderPlus v4.2.7 (database 2026-08-07.1)
- scikit-learn (`MultiLabelBinarizer`) and Matplotlib

Install the Python dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

## Key findings

- Carbapenemase genes were detected in all three genomes in the current merged
  results: *bla*KPC-3 in GCA_000633755.1 and GCA_000694015.1, and *bla*KPC-2 in
  GCA_902704155.1.
- No NDM, OXA-48-like, VIM, or IMP genes were detected in those results.
- All three genomes have an *ompK35* frameshift call. GCA_902704155.1 also has
  an *ompK36* alteration (`ompK36_T136TDT`).
- ESBL-family calls include *bla*SHV and *bla*TEM.
- Fluoroquinolone mutations include `gyrA_S83I` and `parC_S80I`.

## Interpretation

The current AMRFinder profiles include KPC carbapenemase calls, alongside
porin alterations and ESBL-family calls. Interpret these findings in the context
of the source data and AMRFinderPlus database version; they do not support a
non-carbapenemase-producing resistance conclusion for these three genomes.

## Repository structure

- `data/README.md` — public dataset references and local input instructions (data files excluded)
- `results/README.md` — explanation of excluded generated results
- `figures/` — generated visualizations
- `scripts/` — analysis scripts, including `generate_amr_results.py` to merge AMRFinder outputs, `compare_genomes.py` to plot a two-genome AMR comparison, and `carbapenemase_heatmap.py` to compare carbapenemase gene families across all genomes

Generate the default AMR comparison heatmap (GCA_000633755.1 vs GCA_000694015.1) with:

```bash
python3 scripts/compare_genomes.py
```

The SVG is written to `figures/genome_amr_comparison.svg`. Change the genomes or output path with `--genome-a`, `--genome-b`, and `--output`.

Generate the carbapenemase-family heatmap for all genomes in the merged table with:

```bash
python3 scripts/carbapenemase_heatmap.py
```

The SVG is written to `figures/carbapenemase_heatmap.svg`. It shows KPC, NDM,
OXA-48-like, VIM, and IMP families; pale cells indicate that the family was not
detected in the available AMRFinder results. The script uses scikit-learn to
encode gene-family presence/absence and Matplotlib to render the heatmap.

## Reproducibility

The pipeline can be extended to more genomes by running AMRFinderPlus on
additional assemblies and adding their per-genome TSV outputs to the local
`results/amrfinder/` directory. Record the source dataset release and filtering
criteria so results can be traced to the exact public data snapshot.

## Data policy

Public source datasets, downloaded assemblies, and generated result tables are
not committed. `.gitignore` excludes the `data/` and `results/` contents except
for their explanatory README files. The plots in `figures/` are derived
visualizations.

## License

No code license has been selected for this repository.

## Author

[Cheska Marciales]
