# Parser to analyze Rosetta Interface Analyzer concatenated outout
# for protein-protein complex
# In the input:
# Line 1 : complex

import pandas as pd

def parse_interface_analyzer(compname, nmod):
    # make ifmetric file name
    filename = compname + '_' + nmod + '_ifmetric.sc'
    print(f"Analyzing file {filename} ....")

    with open(filename) as f:
        lines = [x.strip() for x in f if x.startswith("SCORE:")]

    header = lines[0].split()[1:]
    values = [l.split()[1:] for l in lines[1:]]

    df = pd.DataFrame(values, columns=header)

    numeric_cols = df.columns[:-1]

    df[numeric_cols] = df[numeric_cols].astype(float)

    comp = df.iloc[0]

    # initialize dictoionar for all metrics
    metrics = {}

    metrics["complex"] = compname
    metrics["nconf"] = nmod
    metrics["dG"] = comp["dG_separated"]
    metrics["dSASA"] = comp["dSASA_int"]
    metrics["HB"] = comp["hbonds_int"]
    metrics["SC"] = comp["sc_value"]
    metrics["unsatHB"] = comp["delta_unsatHbonds"]

    return metrics

# main code
if __name__ == "__main__":

    # complex list file name
    comp_file = "../../comp_list"

    # empty lists for IF metric
    ifmetric_list = []

    # read the complexes into a list
    with open(comp_file, 'r') as cf:
        comp_list = [line.rstrip().removesuffix('.cif') for line in cf.readlines()]
    
    # run a loop over the complexes
    for compname in comp_list:
        # nested loop over 5 AF3 run configs
        for ncon in range(5):
            ifmetric = parse_interface_analyzer(compname, str(ncon))

            # append the dataframes to lists
            ifmetric_list.append(ifmetric)

    # convert the list to dataframe
    ifmetric_master = pd.DataFrame(ifmetric_list)

    # save the df to csv
    ifmetric_master.to_csv('bench20_ifmetric.csv', encoding = 'utf-8', index = False)
    print("Done ... saved IF metric values to csv...")
