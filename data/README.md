# Data sources (not included)

This repository intentionally does not redistribute the public source datasets,
downloaded genome assemblies, phenotype tables, or intermediate sequence
analysis files. Those files are excluded from Git by the repository-level
`.gitignore`.

## Phenotype cohort

The project provenance identifies the phenotype data as the public CABBAGE
database accessed through the EMBL-EBI antimicrobial-resistance portal. Retrieve
the current public data from its source and record the dataset release/version,
download date, and filtering steps when reproducing the analysis. The source
workspace did not retain a stable release identifier or a direct download URL,
so no specific snapshot is claimed here.

## Genome assemblies

Genome assemblies are public NCBI Datasets records. The three genomes present in
the current merged AMRFinder results are:

- [GCA_000633755.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000633755.1/)
- [GCA_000694015.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000694015.1/)
- [GCA_902704155.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_902704155.1/)

Use the [NCBI Datasets portal](https://www.ncbi.nlm.nih.gov/datasets/) to
download the assembly package and consult the associated metadata for the
current record version.

## Local inputs expected by the scripts

To regenerate the results, obtain the source data and run AMRFinderPlus locally.
Place per-genome AMRFinder TSV results in `results/amrfinder/` using filenames
matching `*_genomic_amr.tsv`. Then, from the repository root, run:

```bash
python3 scripts/generate_amr_results.py
python3 scripts/carbapenemase_heatmap.py
python3 scripts/compare_genomes.py
```

All source files and generated result tables in these directories remain local
and are not tracked by Git.
