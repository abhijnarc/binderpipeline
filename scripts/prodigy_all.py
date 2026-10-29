# Parser to analyze Rosetta Interface Analyzer concatenated outout
# for protein-protein complex
# In the input:
# Line 1 : complex

import pandas as pd
import numpy as np
import subprocess

DCUT = 4.5

def prodigy_analyze(compname, nmod):

    # empty dictionary
    metrics = {}
    # make pdb file name
    pdb = compname + '_' + nmod + '_0002.pdb'
    print(f"Analyzing {pdb} ....")

    # run prodigy as a system command
    args = ('/home/user/miniconda3/bin/prodigy ' + pdb +
            ' --selection A B --contact_list --distance-cutoff ' + str(DCUT))
    pres = subprocess.run(args, capture_output = True, shell = True,
                             encoding = "utf-8").stdout
    perr = subprocess.run(args, capture_output = True, shell = True,
                             encoding = "utf-8").stderr

    # get the values
    if len(perr) > 0:
        energy = kd = 9999
    else:
        p_array = pres.splitlines()
        kd = float(p_array[-1].split(':')[1])*1e09
        energy = float(p_array[-2].split(':')[1])

    metrics["complex"] = compname
    metrics["nconf"] = nmod
    metrics["logKdnM"] = np.log(kd)
    metrics['energy'] = energy

    return metrics

# main code
if __name__ == "__main__":

    # complex list file name
    comp_file = "../comp_list"

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
    prod_complex.to_csv('bench20_prod.csv', encoding = 'utf-8', index = False)
    print("Done ... saved prodigy metric values to csv...")
