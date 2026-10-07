# Data sources (not included)

This repository intentionally does not redistribute the public source datasets,
downloaded genome assemblies, phenotype tables, or intermediate sequence
analysis files. Those files are excluded from Git by the repository-level
`.gitignore`.

## Phenotype cohort

The local `phenotype.parquet` input is a public, composite antimicrobial
phenotype export. Its `database` provenance field contains sources including
NCBI antibiogram, NDARO, PATRIC, CDC, and CABBAGE PubMed records. The public
[NCBI Pathogen Detection portal](https://www.ncbi.nlm.nih.gov/pathogens/)
describes the pathogen-genomics and AMR surveillance resource. This repository
does not include the phenotype export; obtain the public records from their
source portal and retain the export's release/version and download date. The
original workspace did not preserve a stable URL or snapshot identifier for
its parquet file, so this repository does not claim a particular release.

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

## Expanding the paired phenotype cohort

If a local copy of the public CABBAGE phenotype parquet is available at
`data/raw/phenotype.parquet`, discover *K. pneumoniae* records with exact AST
labels `ertapenem = resistant` and `imipenem = susceptible` (imipenem-relebactam
is not included) and download assemblies with:

```bash
python3 -m pip install -r requirements.txt
python3 scripts/fetch_ertapenem_R_imipenem_S_genomes.py
```

The script writes the matched AST rows to
`data/processed/ertapenem_R_imipenem_S/phenotype_pairs.tsv` and retrieves
genome FASTAs from the NCBI Datasets API into the same ignored local directory.
It includes only records with a valid NCBI assembly accession in the download.
Records with no assembly accession remain visible in the manifest but cannot
be downloaded as assemblies. The source snapshot available in the original
workspace yielded 32 paired AST/provenance rows, 25 rows with an assembly
accession, and 21 distinct genome assemblies. Seven pairs without assembly
accessions are listed in the phenotype manifest; assemblies were not
reconstructed from their SRA runs. AST source labels, MIC values, and standards
are preserved as supplied in the manifest and are not harmonized by the script.
Outputs are not committed to Git.

### Process the downloaded assemblies

Install AMRFinderPlus 4.2.7 and run the `--plus` workflow on each retrieved
assembly:

```bash
conda create -n kp-amr-processing -c conda-forge -c bioconda python=3.12 ncbi-amrfinderplus=4.2.7
conda run -n kp-amr-processing python -m pip install -r requirements.txt
conda run -n kp-amr-processing amrfinder_update --database data/local_amrfinder_db
bash scripts/process_ertapenem_R_imipenem_S.sh
```

If the named Conda environment already exists, install Python and pip into it:

```bash
conda install -n kp-amr-processing -c conda-forge python=3.12 pip
```

Then install the requirements as shown above.

This saves a per-genome TSV and run metadata under
`results/ertapenem_R_imipenem_S/amrfinder/`. The pipeline resumes completed
AMRFinder runs if interrupted. Source data, sequence assemblies, local database
files, and generated result tables remain local and are not committed.
