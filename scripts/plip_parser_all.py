#!/usr/bin/env python3



import xml.etree.ElementTree as ET
import pandas as pd
from collections import defaultdict
import glob
import os

PREFIX = os.environ.get("PREFIX", "pilot10")
INTERACTION_TYPES = ("HydrogenBond", "Hydrophobic", "SaltBridge", "PiStack")
INTERACTION_COLUMNS = (
    "interaction", "protein_chain", "protein_residue", "protein_resnr",
    "ligand_chain", "ligand_residue", "distance", "donor",
)

##########################################################################
# Parse PLIP XML
##########################################################################

def text(node, tag, default=""):
    """Safely retrieve child text."""
    child = node.find(tag)
    if child is None or child.text is None:
        return default
    return child.text.strip()


##########################################################################
# Main parser
##########################################################################

def parse_plip_xml(xmlfile):

    tree = ET.parse(xmlfile)
    root = tree.getroot()

    interactions = []

    summary = defaultdict(int)

    residue_stats = defaultdict(lambda: defaultdict(int))

    ######################################################################
    # Loop over every binding site
    ######################################################################

    for site in root.findall(".//bindingsite"):

        #ligand_chain = text(site, "chain")

        ##################################################################
        # Hydrogen Bonds
        ##################################################################

        for hb in site.findall(".//hydrogen_bonds/hydrogen_bond"):

            row = {}

            row["interaction"] = "HydrogenBond"

            row["protein_chain"] = text(hb, "reschain")

            row["protein_residue"] = text(hb, "restype")

            row["protein_resnr"] = text(hb, "resnr")

            row["ligand_chain"] = text(hb, "reschain_lig")

            row["ligand_residue"] = text(hb, "restype_lig")

            row["distance"] = text(hb, "dist_h-a")

            row["donor"] = text(hb, "protisdon")

            interactions.append(row)

        ##################################################################
        # Hydrophobic
        ##################################################################

        for hp in site.findall(".//hydrophobic_interactions/hydrophobic_interaction"):

            row = {}

            row["interaction"] = "Hydrophobic"

            row["protein_chain"] = text(hp, "reschain")

            row["protein_residue"] = text(hp, "restype")

            row["protein_resnr"] = text(hp, "resnr")

            row["ligand_chain"] = text(hp, "reschain_lig")

            row["ligand_residue"] = text(hp, "restype_lig")

            row["distance"] = text(hp, "dist")

            interactions.append(row)

        ##################################################################
        # Salt Bridges
        ##################################################################

        for sb in site.findall(".//salt_bridges/salt_bridge"):

            row = {}

            row["interaction"] = "SaltBridge"

            row["protein_chain"] = text(sb, "reschain")

            row["protein_residue"] = text(sb, "restype")

            row["protein_resnr"] = text(sb, "resnr")

            row["ligand_chain"] = text(sb, "reschain_lig")

            row["ligand_residue"] = text(sb, "restype_lig")

            row["distance"] = text(sb, "dist")

            interactions.append(row)

        ##################################################################
        # Pi stacking
        ##################################################################

        for pi in site.findall(".//pi_stacks/pi_stack"):

            row = {}

            row["interaction"] = "PiStack"
            row["protein_chain"] = text(pi, "reschain")
            row["protein_residue"] = text(pi, "restype")
            row["protein_resnr"] = text(pi, "resnr")
            row["ligand_chain"] = text(pi, "reschain_lig")
            row["ligand_residue"] = text(pi, "restype_lig")
            row["distance"] = text(pi, "centdist")

            interactions.append(row)

    ######################################################################
    # Convert to dataframe
    ######################################################################

    df = pd.DataFrame(interactions, columns=INTERACTION_COLUMNS)

    ######################################################################
    # Summary
    ######################################################################

    if df.empty:
        # plip_loop.sh explicitly supplies partners [[A], [B]]. Emit a complete
        # zero-count grid so downstream score/entropy aggregation keeps the pose.
        summary_df = pd.DataFrame(
            [(chain, interaction, 0) for chain in ("A", "B")
             for interaction in INTERACTION_TYPES],
            columns=["protein_chain", "interaction", "count"],
        )
    else:
        summary_df = (
            df.groupby(["protein_chain", "interaction"])
              .size()
              .reset_index(name="count")
        )


    return df, summary_df

##########################################################################
# Run
##########################################################################

if __name__ == "__main__":

    # complex list file name
    comp_file = os.environ.get("COMP_LIST", "../../comp_list")

    # empty lists for interactions and summary
    interactions_list = []
    summary_list = []

    # read the complexes into a list
    with open(comp_file, 'r') as cf:
        comp_list = [line.rstrip().removesuffix('.cif') for line in cf.readlines()]
    
    # run a loop over the complexes
    for compname in comp_list:
        # nested loop over 5 AF3 run configs
        for ncon in range(5):
            # make xml file name
            xmlf = compname + '_' + str(ncon) + '.xml'
            print(f"Analyzing {xmlf} ....")
            interactions, summary = parse_plip_xml(xmlf)
            # insert the complex name and AF3 config in dataframes
            interactions.insert(0, 'complex', compname)
            summary.insert(0, 'complex', compname)
            interactions.insert(1, 'nconf', ncon)
            summary.insert(1, 'nconf', ncon)

            # append the dataframes to lists
            interactions_list.append(interactions)
            summary_list.append(summary)

    # convert the list to dataframes
    interactions_df = pd.concat(interactions_list, ignore_index=True) 
    summary_df = pd.concat(summary_list, ignore_index=True) 

    # save the df's to csv
    interactions_path = f"{PREFIX}_interactions.csv"
    summary_path = f"{PREFIX}_intsummary.csv"
    interactions_df.to_csv(interactions_path, encoding='utf-8', index=False)
    summary_df.to_csv(summary_path, encoding='utf-8', index=False)

    print(f"Done ... wrote {interactions_path} and {summary_path}...")

    # Thats all
