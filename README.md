# UCSD BENG 203 — Assignment 2

This folder contains code, data, and outputs for the Assignment 2 analysis.

## Contents

- `A2_p1_helpers.py`
  - Helper functions for ecDNA segment plotting, gene/exon parsing, expression signal averaging, and transcript binning.
- `A2_p1_main.py`
  - Main script for Part 1 analysis.
  - Generates an ecDNA cartoon, computes gene expression scores, identifies the best LPIN1 transcript, performs Poisson bin testing, and reports expressed non-coding regions.
- `A2_p3_c.py`
  - Computes the minimum Bloom filter size in kilobytes for several false positive error probabilities.
  - Writes results to `BF_size_error.csv` with values rounded to 3 decimals.
- `A2_p3_d.py`
  - Builds a Bloom filter from a FASTQ subset and evaluates query sequences for human-kmer membership.

## Data files

- `ecDNA_segments.bed` — ecDNA segment coordinates.
- `ecDNA_genes.bed` — gene intervals on ecDNA.
- `ecDNA_expression.bg` — background expression signal data.
- `UCSC_genome_table.csv` — refGene transcript table for gene/exon lookups.
- `UCSC_genome_table_d.csv` — non-coding interval table for Part 1(d).
- `BF-database/SRR873426.R1_sub.fastq.txt` — FASTQ subset used to build the Bloom filter in Part 3(d).
- `BF-Queries.txt` — query sequences tested against the Bloom filter.

## Outputs

- `ecDNA_cartoon.png`
  - Visual representation of ecDNA segment layout.
- `LPIN1_transcript_profile.png`
  - LPIN1 exon expression profile with significant bins highlighted.
- `BF_size_error.csv`
  - Bloom filter size results for target false positive rates.

## Requirements

- Python 3
- `pandas`
- `numpy`
- `scipy`
- `matplotlib`
- `mmh3`

Install dependencies with:

```bash
pip install pandas numpy scipy matplotlib mmh3
```

## Usage

Run the main scripts from this directory:

```bash
python A2_p1_main.py
python A2_p3_c.py
python A2_p3_d.py
```

> Run each script from the assignment root so file paths resolve correctly.
