#!/usr/bin/env bash

set -euo pipefail

mkdir -p ref

#########################################################################
# Get the 45S rDNA reference
#########################################################################

# Download the reference file
echo "Downloading U13369.1 (45S rRNA) from GenBank..."
efetch -db nucleotide -id U13369.1 -format fasta > ref/U13369.1.fasta

# Index file, extract and concatenate regions
echo "Concatenating to make rearranged reference..."
samtools faidx ref/U13369.1.fasta
samtools faidx ref/U13369.1.fasta U13369.1:41021-42999 > ref/part1.fa
samtools faidx ref/U13369.1.fasta U13369.1:1-14000   > ref/part2.fa

# Merge the sequences
{
  echo ">45S_rDNA_reference (U13369.1_41021-42999_1-14000)"
  cat "ref/part1.fa" "ref/part2.fa" | grep -v "^>" \
    | awk '
        { seq = seq $0 }
        END {
          for (i = 1; i <= length(seq); i += 70)
            print substr(seq, i, 70)
        }
      '
} > "ref/45S_U13369.1_modified_16kb.fasta"

rm -f \
  "ref/U13369.1.fasta" \
  "ref/U13369.1.fasta.fai" \
  "ref/part1.fa" \
  "ref/part2.fa"

echo "45S rDNA reference built successfully:"
ls -lh ref/45S_U13369.1_modified_16kb.fasta

#########################################################################
# Get the 5S rDNA reference
#########################################################################

echo "Downloading X12811.1 (5S rRNA) from GenBank..."
efetch -db nucleotide -id X12811.1 -format fasta > ref/5S_X12811.1.fasta

echo "5S rDNA reference downloaded successfully:"
ls -lh ref/5S_X12811.1.fasta

# Combine the two fasta files and index
cat ref/5S_X12811.1.fasta ref/45S_U13369.1_modified_16kb.fasta > ref/rDNA_combined.fasta
bwa index ref/rDNA_combined.fasta

#########################################################################
# Download relevant reference files
#########################################################################

# Download the Ensembl reference genome (hg38)
wget -O ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

samtools faidx ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Download the genome annotation (GTF)
wget -O ref/GCF_000001405.40_GRCh38.p14_genomic.gff.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
gunzip ref/GCF_000001405.40_GRCh38.p14_genomic.gff.gz

# Download repetitive sequences from RepeatMasker 
wget -O ref/rmsk_hg38.txt.gz \
  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
zcat ref/rmsk_hg38.txt.gz | awk -v OFS='\t' '{sub(/^chr/, "", $6); print $6, $7, $8, $12}' > ref/repeatmasker_hg38.bed

rm ref/rmsk_hg38.txt.gz

#########################################################################
# Get exonic and intronic regions
#########################################################################

python scripts/build_background_ref.py \
  --genome ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --gff ref/GCF_000001405.40_GRCh38.p14_genomic.gff \
  --rmsk ref/repeatmasker_hg38.bed \
  --outdir ref

cat ref/bg_chr1_exon.bed ref/bg_chr1_intron.bed > ref/bg_chr1_exons_introns.bed
cat ref/bg_acro_exon.bed ref/bg_acro_intron.bed > ref/bg_acro_exons_introns.bed

rm ref/bg_chr1_exon.bed ref/bg_chr1_intron.bed ref/bg_acro_exon.bed ref/bg_acro_intron.bed
