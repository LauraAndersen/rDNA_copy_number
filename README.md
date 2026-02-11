
# Calculate rDNA copy number

This repository contains scripts for estimating ribosomal DNA (rDNA) copy number from whole-genome sequencing data, focusing on the 5S and 45S rDNA.

The analysis workflow consists of the following three steps:

1. **Reference and background preparation**  
   Retrieval and construction of reference rDNA and background sequences.

2. **Sequencing depth estimation**  
   Calculation of read depth within rDNA regions.

3. **Copy number normalization**  
   Derivation of normalized rDNA copy number estimates.


## Workflow

### Step 1 — Build reference and background sequences *(run once)*

This step prepares the rDNA reference sequences (45S and 5S) and a set of
background regions (exonic + intronic) used for normalization/comparison.

**What it does**
- Downloads and modifies the 45S and 5S rDNA sequences
- Downloads reference genome data
- Calls `build_background_ref.py` to extract exonic and intronic background regions

**Run (Slurm)**
```bash
sbatch get_reference_seq.sh
```

---

### Step 2 — Estimate rDNA coverage and background read depth

This step extracts reads originating from rDNA regions and computes position-wise coverage for the 5S and 45S arrays, normalized by background read depth (BRD).

**What it does**
- Slices rDNA regions from the input BAM file
- Converts sliced reads to FASTQ
- Maps reads to rDNA reference sequences
- Estimates background read depth (BRD) from single-copy exons and introns
- Computes position-wise depth normalized by BRD

**Run**
```bash
bash calc_rDNA_depth.sh <bam_file> <project_name>
```

---

### Step 3 — Collect normalized rDNA copy number estimates

This step aggregates depth estimates and derives rDNA copy number for
the 5S and 45S subunits across samples.

**What it does**
- Computes normalized rDNA copy number estimates
- Outputs a summary table for downstream analysis

**Run**
```bash
python scripts/collect_rDNA_CN.py --folder <results_folder> --output <output.tsv>
```

