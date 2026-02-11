#!/usr/bin/env bash
#SBATCH --account=CancerEvolution
#SBATCH --time=4:00:00
#SBATCH --mem=32g
#SBATCH --cpus-per-task=1

set -euo pipefail

host_dir=/faststorage/project/CancerEvolution_shared/Projects/laura/rDNA
mkdir -p ${host_dir}/ref

#########################################################################
# Get the 45S rDNA reference
#########################################################################

# Download the reference file
echo "Downloading U13369.1 (45S rRNA) from GenBank..."
efetch -db nucleotide -id U13369.1 -format fasta > ${host_dir}/ref/U13369.1.fasta

# Index file, extract and concatenate regions
echo "Concatenating to make rearranged reference..."
samtools faidx ${host_dir}/ref/U13369.1.fasta
samtools faidx ${host_dir}/ref/U13369.1.fasta U13369.1:41021-42999 > ${host_dir}/ref/part1.fa
samtools faidx ${host_dir}/ref/U13369.1.fasta U13369.1:1-14000   > ${host_dir}/ref/part2.fa

# Merge the sequences
{
  echo ">45S_rDNA_reference (U13369.1_41021-42999_1-14000)"
  cat "${host_dir}/ref/part1.fa" "${host_dir}/ref/part2.fa" | grep -v "^>" \
    | awk '
        { seq = seq $0 }
        END {
          for (i = 1; i <= length(seq); i += 70)
            print substr(seq, i, 70)
        }
      '
} > "${host_dir}/ref/45S_U13369.1_modified_16kb.fasta"

rm -f \
  "${host_dir}/ref/U13369.1.fasta" \
  "${host_dir}/ref/U13369.1.fasta.fai" \
  "${host_dir}/ref/part1.fa" \
  "${host_dir}/ref/part2.fa"

echo "45S rDNA reference built successfully:"
ls -lh ${host_dir}/ref/45S_U13369.1_modified_16kb.fasta

#########################################################################
# Get the 5S rDNA reference
#########################################################################

echo "Downloading X12811.1 (5S rRNA) from GenBank..."
efetch -db nucleotide -id X12811.1 -format fasta > ${host_dir}/ref/5S_X12811.1.fasta

echo "5S rDNA reference downloaded successfully:"
ls -lh ${host_dir}/ref/5S_X12811.1.fasta

# Combine the two fasta files and index
cat ${host_dir}/ref/5S_X12811.1.fasta ${host_dir}/ref/45S_U13369.1_modified_16kb.fasta > ${host_dir}/ref/rDNA_combined.fasta
bwa index ${host_dir}/ref/rDNA_combined.fasta

#########################################################################
# Download relevant reference files
#########################################################################

# Download the Ensembl reference genome (hg38)
wget -O ${host_dir}/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
  https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip ${host_dir}/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

samtools faidx ${host_dir}/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# Download the genome annotation (GTF)
wget -O ${host_dir}/ref/GCF_000001405.40_GRCh38.p14_genomic.gff.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
gunzip ${host_dir}/ref/GCF_000001405.40_GRCh38.p14_genomic.gff.gz

# Download repetitive sequences from RepeatMasker 
wget -O ${host_dir}/ref/rmsk_hg38.txt.gz \
  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
zcat ${host_dir}/ref/rmsk_hg38.txt.gz | awk -v OFS='\t' '{sub(/^chr/, "", $6); print $6, $7, $8, $12}' > ${host_dir}/ref/repeatmasker_hg38.bed

rm ${host_dir}/ref/rmsk_hg38.txt.gz

#########################################################################
# Get exonic and intronic regions
#########################################################################

python scripts/build_background_ref.py \
  --genome ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --gff ref/GCF_000001405.40_GRCh38.p14_genomic.gff \
  --rmsk ref/repeatmasker_hg38.bed \
  --outdir /faststorage/project/CancerEvolution_shared/Projects/laura/rDNA/ref

cat ${host_dir}/ref/bg_chr1_exon.bed ${host_dir}/ref/bg_chr1_intron.bed > ${host_dir}/ref/bg_chr1_exons_introns.bed
cat ${host_dir}/ref/bg_acro_exon.bed ${host_dir}/ref/bg_acro_intron.bed > ${host_dir}/ref/bg_acro_exons_introns.bed
