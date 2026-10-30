# Parser to analyze Rosetta Interface Analyzer concatenated outout
# for protein-protein complex
# In the input:
# Line 1 : complex

import pandas as pd
import numpy as np
import subprocess
import os

DCUT = 4.5
PREFIX = os.environ.get("PREFIX", "pilot10")

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

    # complex list file name
    comp_file = os.environ.get("COMP_LIST", "../comp_list")

    # empty lists for prodigy metric
    prod_list = []

    # read the complexes into a list
    with open(comp_file, 'r') as cf:
        comp_list = [line.rstrip().removesuffix('.cif') for line in cf.readlines()]
    
    # run a loop over the complexes
    for compname in comp_list:
        # nested loop over 5 AF3 run configs
        for ncon in range(5):
            prod_metric = prodigy_analyze(compname, str(ncon))

            # append the dataframes to lists
            prod_list.append(prod_metric)

    # convert the list to dataframe
    prod_master = pd.DataFrame(prod_list)

    # now we agregate the metrics over the complex
    cols_to_mean = prod_master.columns[2:]
    prod_complex = ( prod_master.groupby(prod_master.columns[0])[cols_to_mean]
                        .mean().reset_index() )

    # save the df to csv
    output_file = f"{PREFIX}_prodigy.csv"
    prod_complex.to_csv(output_file, encoding='utf-8', index=False)
    print(f"Done ... saved PRODIGY metrics to {output_file}...")
