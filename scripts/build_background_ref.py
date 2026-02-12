#!/usr/bin/env python3
import argparse, os, subprocess
import pyranges as pr
import pandas as pd
import pysam
import numpy as np

##############################################################
# Functions
##############################################################

def compute_introns_from_transcripts_minus_exons(gff_gr, exons_all):
    # Filter transcripts
    tx = gff_gr[gff_gr.Feature == "mRNA"].copy()

    results = []
    cols = ["Chromosome","Start","End","Strand"]

    # Loop through genes in the df and their transcripts 
    for gene, gene_sub in tx.df.groupby("gene"):
        # Get all exons of the gene
        ex_gene_sub = pr.PyRanges(exons_all.df[exons_all.df["gene"] == gene])
        if ex_gene_sub.empty:
            continue

        gene_introns = []

        # For each transcript, compute introns
        for trans, trans_sub in gene_sub.groupby('ID'):
            # Get exons of the transcript
            sub_ex = ex_gene_sub[ex_gene_sub.Parent == trans]

            if sub_ex.empty:
                continue
            # Compute introns
            intron = pr.PyRanges(trans_sub[cols]).subtract(sub_ex[cols])

            # Get exons from other transcripts, and remove introns that overlaps other exons
            other_trans = ex_gene_sub[ex_gene_sub.Parent != trans]
            if not other_trans.empty:
                intron = intron.overlap(other_trans[cols], invert=True)
            # Add the intron to the list
            if not intron.empty:
                gene_introns.append(intron)

        if gene_introns:
            gene_introns = pr.concat(gene_introns).merge()
            gene_introns = gene_introns.df
            gene_introns["gene"] = gene
            results.append(gene_introns)

    return pr.PyRanges(pd.concat(results, ignore_index=True)).sort()

def pick_largest_exon_per_overlap(gr_exon):
    """
    # Inputs: PyRanges with feature == 'exon'; has gene column
    """
    # Cluster overlapping exons per chromosome (PyRanges.cluster). We also keep gene_id in the clustering key.
    clusters = gr_exon.cluster(by="gene")  # adds 'Cluster' per chromosome
    df = clusters.df.copy()

    # Keep the largest exon within each (gene_id, Cluster)
    df["len"] = (df["End"] - df["Start"]).astype(int)
    idx = df.groupby(["gene", "Cluster"])["len"].idxmax()
    picked = pr.PyRanges(df.loc[idx, ["Chromosome", "Start", "End", "Strand", "gene"]].sort_values(["Chromosome", "Start", "End"]))

    return picked

def load_repeatmasker_bed(rmsk_path):
    df = pd.read_csv(rmsk_path, 
                     sep="\t", 
                     header=None,
                     usecols=[0, 1, 2],
                     dtype={0: "string", 1: np.int32, 2: np.int32},
                     names=["Chromosome", "Start", "End"])
    
    return pr.PyRanges(df).merge()

def load_gff_and_get_coding_transcripts(gff_path):
    refseq_to_chr = {
        "NC_000001.11": "1", "NC_000002.12": "2", "NC_000003.12": "3", "NC_000004.12": "4",
        "NC_000005.10": "5", "NC_000006.12": "6", "NC_000007.14": "7", "NC_000008.11": "8",
        "NC_000009.12": "9", "NC_000010.11": "10", "NC_000011.10": "11", "NC_000012.12": "12",
        "NC_000013.11": "13", "NC_000014.9": "14", "NC_000015.10": "15", "NC_000016.10": "16",
        "NC_000017.11": "17", "NC_000018.10": "18", "NC_000019.10": "19", "NC_000020.11": "20",
        "NC_000021.9": "21", "NC_000022.11": "22", "NC_000023.11": "X", "NC_000024.10": "Y"
    }

    gff_df = pr.read_gff3(gff_path, as_df=True)

    # Rename RefSeq accessions to chromosome numbers
    gff_df["Chromosome"] = gff_df["Chromosome"].replace(refseq_to_chr)

    # Get a list of protein coding genes
    genes_pc = gff_df[gff_df.gene_biotype == "protein_coding"]
    gene_ids = set(genes_pc["ID"].dropna())

    # Filter mRNAs from protein coding genes (these are the ones than contains a coding region (CDS), which transcripts does not)
    tx = gff_df[gff_df.Feature == "mRNA"].copy()
    tx = tx[tx["Parent"].isin(gene_ids)]
    tx_ids = set(tx["ID"].dropna())

    # Make a final df with mRNA and exons
    tx_exons = gff_df[gff_df.Feature.isin(["mRNA", "exon"])].copy()
    tx_exons = tx_exons[(tx_exons['ID'].isin(tx_ids)) | (tx_exons['Parent'].isin(tx_ids))]

    tx_exons = tx_exons[['Chromosome', 'Feature', 'Start', 'End', 'Strand', 'ID', 'Parent', 'gene']].copy()

    return pr.PyRanges(tx_exons)

def build_candidates(gff_path, rmsk_path):
    # Load the gtf with genome annotation of Ensembl
    gff_gr = load_gff_and_get_coding_transcripts(gff_path)

    # Extracted exons and introns on chromosome 1, 13, 14, 15, 21 and 22
    gff_gr = gff_gr[gff_gr.Chromosome.isin(["1", "13","14","15","21","22"])]

    # Extract exons of protein coding transcripts (mRNAs)
    exons_all = gff_gr[gff_gr.Feature == "exon"]

    # Pick largest exon in overlapping isoforms (by gene + cluster)
    exons = pick_largest_exon_per_overlap(exons_all)

    # Extract introns from transcript spans minus union of exons
    introns = compute_introns_from_transcripts_minus_exons(gff_gr, exons_all)

    # Build exon+intron candidates per set
    ex_chr1 = exons[exons.Chromosome == '1']
    in_chr1 = introns[introns.Chromosome == '1']

    acro_chroms = ["13","14","15","21","22"]
    ex_acro = exons[exons.Chromosome.isin(acro_chroms)]
    in_acro = introns[introns.Chromosome.isin(acro_chroms)]

    # Subtract repeats
    rmsk_gr = load_repeatmasker_bed(rmsk_path)

    cands_chr1_exon, cands_chr1_intron = ex_chr1.overlap(rmsk_gr, invert=True).merge(), in_chr1.overlap(rmsk_gr, invert=True).merge()
    cands_acro_exon, cands_acro_intron = ex_acro.overlap(rmsk_gr, invert=True).merge(), in_acro.overlap(rmsk_gr, invert=True).merge()

    return cands_chr1_exon, cands_chr1_intron, cands_acro_exon, cands_acro_intron

def len_filter(gr, minlen=300, maxlen=10000):
    df = gr.df.copy()
    df["len"] = (df["End"] - df["Start"]).astype(int)
    df = df[(df["len"] >= minlen) & (df["len"] <= maxlen)]
    return pr.PyRanges(df[["Chromosome","Start","End"]]).sort()

def trim(gr, trim_bp=50):
    df = gr.df.copy()
    df["Start"] = df["Start"].astype(int) + trim_bp
    df["End"]   = df["End"].astype(int)   - trim_bp
    df = df[df["End"] > df["Start"]]
    return pr.PyRanges(df[["Chromosome","Start","End"]]).sort()

def write_fasta_from_pyranges(gr, genome_fa, out_fa):
    """Write a FASTA where each record header is 'chrom:start-end' matching BED coords."""
    fa = pysam.FastaFile(genome_fa)
    with open(out_fa, "w") as out:
        for _, r in gr.df.iterrows():
            chrom = str(r.Chromosome)
            start = int(r.Start)
            end   = int(r.End)
            seq = fa.fetch(chrom, start, end)
            out.write(f">{chrom}:{start}-{end}\n")
            # wrap sequence 60 nt per line
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")
    fa.close()

def run_makeblastdb(genome_fa, db_prefix):
    os.makedirs(os.path.dirname(db_prefix), exist_ok=True)
    subprocess.run([
        "makeblastdb",
        "-in", genome_fa,
        "-dbtype", "nucl",
        "-parse_seqids",
        "-out", db_prefix
    ], check=True)

def run_blastn(query_fa, db_prefix, out_tsv, threads=8, evalue=1e-6):
    subprocess.run([
        "blastn",
        "-query", query_fa,
        "-db", db_prefix,
        "-out", out_tsv,
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore qstart qend sstart send",
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-task", "megablast"
    ], check=True)

def filter_single_copy_from_blast(gr, blast_tsv):
    # Load the BLAST output file and get the ones that only aligns with itself
    df = pd.read_csv(blast_tsv, sep="\t", header=None,
                     names=["qseqid","sseqid","pident","length","evalue","bitscore","qstart","qend","sstart","send"])
    counts = df.groupby("qseqid").size()
    keep_ids = set(counts[counts == 1].index)

    gdf = gr.df.copy()
    gdf["id"] = gdf["Chromosome"].astype(str) + ":" + gdf["Start"].astype(str) + "-" + gdf["End"].astype(str)
    gdf = gdf[gdf["id"].isin(keep_ids)]

    return pr.PyRanges(gdf[["Chromosome","Start","End"]]).sort()


##############################################################
# Code
##############################################################

def main():
    parser = argparse.ArgumentParser(description="Build hg38 background regions (single-copy exons/introns)")
    parser.add_argument("--genome", required=True, help="GRCh38 FASTA (indexed; *.fai present)")
    parser.add_argument("--gff",    required=True, help="Annotation GFF for GRCh38")
    parser.add_argument("--rmsk",   required=True, help="RepeatMasker BED for hg38")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1) Build exon/intron candidates and subtract repeats
    print("Building exon/intron candidates...", flush=True)
    cands_chr1_exon, cands_chr1_intron, cands_acro_exon, cands_acro_intron = build_candidates(args.gff, args.rmsk)

    # 2) Run makeblastdb (once per genome path)
    print("Running blast...", flush=True)
    db_prefix = os.path.join(args.outdir, "hg38_blastdb", "hg38")
    if not (os.path.exists(db_prefix + ".nhr") or os.path.exists(db_prefix + ".nin") or os.path.exists(db_prefix + ".nsq")):
        run_makeblastdb(args.genome, db_prefix)

    # 3) Write FASTA files as input for BLAST (headers: chrom:start-end) and run BLAST
    candidates = {
        "chr1_exon":  cands_chr1_exon,
        "chr1_intron": cands_chr1_intron,
        "acro_exon":   cands_acro_exon,
        "acro_intron": cands_acro_intron,
    }
    for name, gr in candidates.items():
        out_fa = f"{args.outdir}/{name}_candidates.fa"
        out_tsv = f"{args.outdir}/{name}_candidates.tsv"
        out_bed = f"{args.outdir}/bg_{name}.bed"

        # Run blast
        write_fasta_from_pyranges(gr, args.genome, out_fa)
        run_blastn(out_fa, db_prefix, out_tsv)

        bg_untrim = filter_single_copy_from_blast(gr, out_tsv)
        filter_df = len_filter(bg_untrim)
        filter_trim_df = trim(filter_df)

        final_df = filter_trim_df.df[["Chromosome","Start","End"]]
        final_df["Chromosome"] = 'chr' + final_df["Chromosome"].astype(str)
        final_df.to_csv(out_bed, sep="\t", header=False, index=False)

        # Clean up intermediate files
        os.remove(out_fa)
        os.remove(out_tsv)

        print(f"[ok] Wrote:\n  {out_bed}")


if __name__ == "__main__":
    main()


