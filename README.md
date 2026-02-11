
# Calculate rDNA copy number

This repository contains scripts for estimating ribosomal DNA (rDNA) copy number from whole-genome sequencing data, focusing on the 5S and 45S rDNA.

The analysis workflow consists of the following three steps:

1. **Reference and background preparation**  
   Retrieval and construction of reference rDNA and background sequences.

2. **Sequencing depth estimation**  
   Calculation of read depth within rDNA regions.

3. **Copy number normalization**  
   Derivation of normalized rDNA copy number estimates.

---

## Workflow

### Step 1 — Build reference and background sequences *(only run once)*

This step prepares the rDNA reference sequences (45S and 5S) and a set of
background regions (exonic + intronic) used for normalization/comparison.

**What it does**
- Downloads and modifies the 45S and 5S rDNA sequences
- Downloads reference genome data
- Calls `build_background_ref.py` to extract exonic and intronic background read depth (BRD) regions

**Run:**
```bash
sbatch get_reference_seq.sh
```

**Outputs:**
```
ref/
├── repeatmasker_hg38.bed
├── rDNA_slice_regions.hg38.bed
├── rDNA_combined.fasta
├── rDNA_combined.fasta.amb
├── rDNA_combined.fasta.ann
├── rDNA_combined.fasta.bwt
├── rDNA_combined.fasta.pac
├── rDNA_combined.fasta.sa
├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
├── Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
├── hg38_blastdb/
├── GCF_000001405.40_GRCh38.p14_genomic.gff
├── bg_chr1_exons_introns.bed
├── bg_acro_exons_introns.bed
├── 5S_X12811.1.fasta
└── 45S_U13369.1_modified_16kb.fasta
```

---

### Step 2 — Estimate rDNA coverage and background read depth

This step extracts reads originating from rDNA regions and computes position-wise coverage for the 5S and 45S arrays, normalized by BRD.

**What it does**
- Slices rDNA regions from the input BAM file
- Converts sliced reads to FASTQ
- Maps reads to rDNA reference sequences
- Estimates BRD from single-copy exons and introns
- Computes position-wise BRD normalized depth for 5S and 45S rDNA

**Inputs:**
- `sample_name.bam` — WGS BAM file
- `project_name` — Project identifier (e.g. `TCGA_BLCA`). Used to name the output directory.

**Run:**
```bash
bash calc_rDNA_depth.sh sample_name.bam project_name
```

**Outputs:**
```
└── output/
    └── project_name/
        ├── sample_name.5S.BRD_norm_depth.tsv
        └── sample_name.45S.BRD_norm_depth.tsv
```

These files contain position-wise rDNA depth estimates normalized by BRD for the 5S and 45S arrays.

---

### Step 3 — Collect normalized rDNA copy number estimates

This step aggregates depth estimates and derives rDNA copy number for
the 5S and 45S subunits across samples.

**What it does**
- Computes normalized rDNA copy number estimates for 5s and 45S rDNA subunits
- Outputs a summary table for downstream analysis

**Inputs:**
- `project_name` — Project identifier (e.g. `TCGA_BLCA`).
- `output.tsv` — filename for output summary tsv file.

**Run:**
```bash
python scripts/collect_rDNA_CN.py --folder project_name --output output.tsv
```

**Outputs:**
```
└── output/
    ├── summary_files
        └── output.tsv
    └── project_name/
        ├── sample_name.5S.BRD_norm_depth.tsv
        └──  sample_name.45S.BRD_norm_depth.tsv
```
