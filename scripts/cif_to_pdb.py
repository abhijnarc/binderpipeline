#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  1 18:26:59 2025

Code to convert a set of cif files to pdb

@author: user
"""
# libraries to import
from Bio.PDB import MMCIFParser, PDBIO
import pandas as pd
import glob
import os


# function for the conversion

# function for the conversion
def cif2pdb(cif_file, pdb_file):

    # Initialize the mmCIF parser
    parser = MMCIFParser()

    # Load the structure from the CIF file
    structure = parser.get_structure("complex", cif_file)

    # Initialize the PDB writer
    io = PDBIO()

    # Set the structure and save as a PDB file
    io.set_structure(structure)
    io.save(pdb_file)

    print(f"Converted {cif_file} to {pdb_file} ....")


# main code

if __name__ == "__main__":
    
    # initialize counter
    pdbcount = 0
    # iterate through all 4-letter name directories
    for dname in glob.glob("????/"):
        # clone names and its lower case
        clone = os.path.dirname(dname)
        lclone = clone.lower()
        for cif in glob.glob(os.path.join(dname, "*.cif")):
            bcif = os.path.basename(cif)
            model_num = bcif.removesuffix(".cif").split('_')[-1]
            pdb = "AF3_models/" + clone + "_" + model_num + ".pdb"
            cif2pdb(cif, pdb)
            # increment counter
            pdbcount += 1

    print(f"All done ... saved {pdbcount} pdb files..")


# Thats all
