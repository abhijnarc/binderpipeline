# Parser to analyze Rosetta Interface Analyzer concatenated outout
# for protein-protein complex
# In the input:
# Line 1 : complex

import pandas as pd
import numpy as np
import subprocess
import os
import argparse
from pathlib import Path

DCUT = 4.5
PREFIX = os.environ.get("PREFIX", "bench20")
MODEL_COUNT = int(os.environ.get("MODEL_COUNT", "5"))
if MODEL_COUNT not in range(1, 6):
    raise SystemExit("MODEL_COUNT must be between 1 and 5")
ROOT = Path(__file__).resolve().parents[1]

def prodigy_analyze(compname, nmod):

    # empty dictionary
    metrics = {}
    # make pdb file name
    pdb = compname + '_' + nmod + '_0002.pdb'
    print(f"Analyzing {pdb} ....")

    # Run PRODIGY once. Its dissociation constant is reported in molar units.
    args = [
        "/home/wetlab/miniforge3/envs/ppi-affinity/bin/prodigy", pdb,
        "--selection", "A", "B", "--contact_list",
        "--distance-cutoff", str(DCUT),
    ]
    pres = subprocess.run(args, capture_output=True, text=True, check=False)

    energy = kd = np.nan  # Failed parses replace the old 9999 sentinel with NaN.
    if pres.returncode == 0:
        affinity_line = next((line for line in pres.stdout.splitlines()
                              if "Predicted binding affinity" in line), None)
        kd_line = next((line for line in pres.stdout.splitlines()
                        if "Predicted dissociation constant" in line), None)
        try:
            energy = float(affinity_line.rsplit(":", 1)[1].strip())
            kd = float(kd_line.rsplit(":", 1)[1].strip().split()[0])
            if kd <= 0:
                energy = kd = np.nan
        except (AttributeError, ValueError, IndexError):
            energy = kd = np.nan

    metrics["complex"] = compname
    metrics["nconf"] = nmod
    metrics["logKd"] = np.log10(kd) if pd.notna(kd) else np.nan
    metrics['energy'] = energy

    return metrics

# main code
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--crystal", action="store_true", help="run once on native crystal structures")
    parser.add_argument("--selection", default=str(ROOT / "bench20/selection.tsv"))
    args = parser.parse_args()
    comp_file = os.environ.get("COMP_LIST", str(ROOT / "comp_list"))

    # empty lists for prodigy metric
    prod_list = []

    # read the complexes into a list
    with open(comp_file, 'r') as cf:
        comp_list = [line.rstrip().removesuffix('.cif') for line in cf.readlines()]
    
    if args.crystal:
        import csv
        with open(args.selection, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        for row in rows:
            compname = row["id"]
            pdb = str(ROOT / "pdbbind/P-P" / f"{compname}_complex.pdb")
            args_cmd = ["/home/wetlab/miniforge3/envs/ppi-affinity/bin/prodigy", pdb, "--selection", row["chain1"], row["chain2"], "--contact_list", "--distance-cutoff", str(DCUT)]
            result = subprocess.run(args_cmd, capture_output=True, text=True, check=False)
            energy = kd = np.nan
            if result.returncode == 0:
                a = next((x for x in result.stdout.splitlines() if "Predicted binding affinity" in x), None)
                k = next((x for x in result.stdout.splitlines() if "Predicted dissociation constant" in x), None)
                try:
                    energy = float(a.rsplit(":", 1)[1].strip()); kd = float(k.rsplit(":", 1)[1].strip().split()[0])
                    if kd <= 0: energy = kd = np.nan
                except (AttributeError, ValueError, IndexError): energy = kd = np.nan
            prod_list.append({"complex": compname, "energy": energy, "prodigy_logKd": np.log10(kd) if pd.notna(kd) else np.nan})
    else:
        for compname in comp_list:
            for ncon in range(MODEL_COUNT):
                prod_list.append(prodigy_analyze(compname, str(ncon)))

    # convert the list to dataframe
    prod_master = pd.DataFrame(prod_list)

    # now we agregate the metrics over the complex
    cols_to_mean = prod_master.columns[2:]
    prod_complex = ( prod_master.groupby(prod_master.columns[0])[cols_to_mean]
                        .mean().reset_index() )

    # save the df to csv
    if args.crystal:
        prefix_path = Path(PREFIX)
        output_dir = prefix_path if prefix_path.is_absolute() or len(prefix_path.parts) > 1 else ROOT / prefix_path
        output_file = str(output_dir / "prodigy_crystal.csv")
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_file = f"{PREFIX}_prodigy.csv"
    prod_complex.to_csv(output_file, encoding='utf-8', index=False)
    print(f"Done ... saved PRODIGY metrics to {output_file}...")
