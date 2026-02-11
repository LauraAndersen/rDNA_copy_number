import argparse, os, subprocess
from pathlib import Path
import pandas as pd
import numpy as np

"""
Run as:
python scripts/collect_rDNA_CN.py --folder output/TCGA_PRAD/ --output output/summary_files/TCGA_PRAD_rDNA_copy_number.tsv
"""

def annotate_45s_unit(pos: pd.Series) -> pd.Series:
    """
    Vectorized mapping of positions to 45S subunit labels.
    Returns a categorical Series with values in {18S, 5.8S, 28S} or NaN.
    """
    LABELS = ["18S", "5.8S", "28S"]

    # Use np.select for 3 ranges (fast, explicit)
    conditions = [
        pos.between(5636, 7506), # 18S
        pos.between(8602, 8758), # 5.8S
        pos.between(9914, 14948), # 28S
    ]
    return pd.Series(np.select(conditions, LABELS, default=None), index=pos.index)

def collect_rDNA_CN(project_folder):

    # Define patient IDs
    patient_ids = list({f.name.split(".")[0] for f in Path(project_folder).glob("*.tsv")})

    rows = []

    # Loop over every patient 
    for pt in patient_ids:

        file_5s = f"{project_folder}/{pt}.5S.BRD_norm_depth.tsv"
        file_45s = f"{project_folder}/{pt}.45S.BRD_norm_depth.tsv"

        # --- 5S ---
        df_5S = pd.read_csv(
            file_5s,
            sep="\t",
            header=None,
            names=["name", "pos", "cov"],
        )
        CN_5S = df_5S["cov"].mean()
        rows.append({"sample_name": pt, "rDNA_gene": "5S", "CN": CN_5S})

        # --- 45S ---
        df_45S = pd.read_csv(
            file_45s,
            sep="\t",
            header=None,
            names=["name", "pos", "cov"],
        )

        df_45S["rDNA_gene"] = annotate_45s_unit(df_45S["pos"])
        df_45S = df_45S.dropna(subset=["rDNA_gene"])
        
        # Mean CN per subunit
        subunit_means = (
            df_45S.groupby("rDNA_gene", observed=True)["cov"]
            .mean()
            .rename("CN")
            .reset_index()
        )

        # Add rows for 18S/5.8S/28S
        for _, r in subunit_means.iterrows():
            rows.append({"sample_name": pt, "rDNA_gene": r["rDNA_gene"], "CN": r["CN"]})

        # Total 45S = mean of the three subunit means (matches your intent)
        rows.append({"sample_name": pt, "rDNA_gene": "45S", "CN": subunit_means["CN"].mean()})

    return pd.DataFrame(rows)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", required=True, help="")
    parser.add_argument("--output", required=True, help="")
    args = parser.parse_args()

    CN_df = collect_rDNA_CN(args.folder)

    # Save to tsv file
    CN_df.to_csv(args.output,
                 sep="\t",
                 index=False
                 )

if __name__ == "__main__":
    main()

