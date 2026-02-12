#!/bin/bash

##############################################################
# Run as: bash calc_rDNA_depth.sh bam project
# test bam: /faststorage/project/CancerEvolution_shared/Projects/chongming/test_10_BLCA/TCGA-4Z-AA7S-10A.bam

# For multiple runs run as:
# for bam in /project/CancerEvolution/Datasets/TCGA_WGS/randiip/data/TCGA-OV/*.bam; do
#     sample=$(basename "$bam" .bam)

#     sbatch -o logs/${sample}.log -e logs/${sample}.log \
#            scripts/calc_rDNA_depth.sh "$bam" TCGA_OV
# done
##############################################################

# Load conda environment
set -euo pipefail

THREADS=8 # $SLURM_CPUS_PER_TASK

# Define input
BAM="$1"
PROJECT_ID="$2"

# Define output
REF_DIR="ref"
# LOGS="logs"
OUT_DIR="output/${PROJECT_ID}"

# Create directories 
# mkdir -p "${LOGS}"
mkdir -p "${OUT_DIR}"

# Sample name from BAM filename
BAM_BASENAME="$(basename "${BAM}")"
SAMPLE="${BAM_BASENAME%%.*}"

# Reference files
RDNA_REF="${REF_DIR}/rDNA_combined.fasta"
BG_CHR1_BED="ref/bg_chr1_exons_introns.bed"          # used for 5S CN
BG_OTHER_BED="ref/bg_acro_exons_introns.bed"         # used for 45S CN

# rDNA slice coordinates in the *genome-aligned* BAM (as in the paper)
RDNASCAF_45S="chrUn_GL000220v1"                             # 45S reads mostly here
FIVE_S_REGION="chr1:226743523-231781906"                    # 1q42 +/- 2 Mb around 5S locus

# Names of rDNA contigs inside rDNA_combined.fasta
RDNA_5S_CONTIG="X12811.1"
RDNA_45S_CONTIG="45S_rDNA_reference"

# Check if output already exists
FINAL_5S="${OUT_DIR}/${SAMPLE}.5S.BRD_norm_depth.tsv"
FINAL_45S="${OUT_DIR}/${SAMPLE}.45S.BRD_norm_depth.tsv"

if [[ -f "${FINAL_5S}" && -f "${FINAL_45S}" ]]; then
  echo "[INFO] Final output already exists for ${SAMPLE}. Skipping."
  exit 0
fi

# Check all files are there
for f in "${BAM}" "${RDNA_REF}" "${BG_CHR1_BED}" "${BG_OTHER_BED}"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: Required file not found: $f" >&2
    exit 1
  fi
done

# Check if the bam file is indexed
if [[ ! -f "${BAM%.bam}.bai" && ! -f "${BAM}.bai" ]]; then 
  samtools index "${BAM}"; 
fi

echo "========================================"
echo "[START] $(date)"
echo "[INFO] Processing sample: ${SAMPLE}"
echo "[INFO] Input BAM: ${BAM}"
echo "========================================"

###############################
# Step 1: Slice rDNA regions
###############################

SAMPLE_FQ="${OUT_DIR}/${SAMPLE}.rDNA_slice.fastq"

echo "[INFO] Slicing BAM"
samtools view --fetch-pairs  -b "${BAM}" "$RDNASCAF_45S" "$FIVE_S_REGION" \
  | samtools fastq -N -o "${SAMPLE_FQ}"


########################################
# Step 2: Map rDNA-sliced reads to rDNA reference
########################################

RDNA_BAM="${OUT_DIR}/${SAMPLE}.rDNA_onRef.sorted.bam"

echo "[INFO] Mapping FASTQ to rDNA reference"
bwa mem -t "$THREADS" -p ${RDNA_REF} ${SAMPLE_FQ} \
  | samtools sort -o "${RDNA_BAM}"
samtools index "${RDNA_BAM}"


########################################
# Step 3. Estimate BRD for single-copy regions
########################################

echo "[INFO] Estimating background read depth (BRD)"
BRD_CHR1=$(samtools depth -b "${BG_CHR1_BED}" "${BAM}" \
  | awk '{sum+=$3; n++} END{print sum/n}')

BRD_OTHER=$(samtools depth -b "${BG_OTHER_BED}" "${BAM}" \
  | awk '{sum+=$3; n++} END{print sum/n}')

echo "[INFO] BRD chr1   (5S) : ${BRD_CHR1}"
echo "[INFO] BRD chr13+ (45S): ${BRD_OTHER}"


########################################
# Step 4: Get position wise depth for the aligned files
########################################

echo "[INFO] Computing normalized position-wise coverage on rDNA reference"

# 5S: normalize by BRD from chr1 single-copy regions
samtools depth -a -r "${RDNA_5S_CONTIG}" "${RDNA_BAM}" \
  | awk -v brd="${BRD_CHR1}" 'BEGIN{OFS="\t"} {print $1,$2,$3/brd}' \
  > "${OUT_DIR}/${SAMPLE}.5S.BRD_norm_depth.tsv"

# 45S: normalize by BRD from chr13/14/15/21/22 single-copy regions
samtools depth -a -r "${RDNA_45S_CONTIG}" "${RDNA_BAM}" \
  | awk -v brd="${BRD_OTHER}" 'BEGIN{OFS="\t"} {print $1,$2,$3/brd}' \
  > "${OUT_DIR}/${SAMPLE}.45S.BRD_norm_depth.tsv"

rm -f "${SAMPLE_FQ}"
rm -f "${RDNA_BAM}" "${RDNA_BAM}.bai"

echo "[DONE]"



