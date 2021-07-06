import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
import fileinput
from shutil import copyfile
import os

import get_data
import solvate
import intercalate_n_merge

# Get info about the mutation
[nres, resname] = get_data.get_res_resname('../windows_res')

# Get info about the amyloid structure
og_structure = PandasPdb().read_pdb('./leap/5kk3.pdb')
[res_min, res_max, nchains] = get_data.get_amyloid_info(og_structure)
tiatoms = get_data.get_tiatoms(og_structure, nres, resname)
tiatoms_total = tiatoms*nchains

# Arrange how chains will be mutated and create mutants
n_mutants = tiatoms_total//500 + (tiatoms_total % 500 > 0)
mutants = get_data.create_mutants(og_structure, n_mutants, nchains, nres, resname)
og_structure = PandasPdb().read_pdb('./leap/5kk3.pdb')

# Calculate charges and get order of solvation
charges = [get_data.get_charge(mutant) for mutant in mutants]
charges = [get_data.get_charge(og_structure)] + charges
mutants = [og_structure] + mutants
order = get_data.get_order(charges)

# Solvate and strip water and ions
for index, value in enumerate(order):
    if index == 0:
        if index == value:
            structure1 = [mutants[index], "WT"]
            structure2 = [mutants[index+1], f"mutant_{index+1}"]
        else:
            structure2 = [mutants[index], "WT"]
            structure1 = [mutants[index+1], f"mutant_{index+1}"]
    if index != 0:
        if index == value:
            structure1 = [mutants[index], f"mutant_{index}"]
            structure2 = [mutants[index+1], f"mutant_{index+1}"]
        else:
            structure2 = [mutants[index], f"mutant_{index}"]
            structure1 = [mutants[index+1], f"mutant_{index+1}"]
    solvate.solvate(structure1, structure2)
    solvate.strip(structure1, structure2)
    solvate.solvate_ions_pair(structure1, structure2)
    
    # Intercalate and merged structures
    intercalate_n_merge.intercalate(structure1, structure2)
    intercalate_n_merge.merge(nres, res_min, res_max, nchains, structure1, structure2)
    [reswt_str, resmut_str] = intercalate_n_merge.in_files_setup(nres, resname, res_min, res_max, nchains, structure1, structure2)
    
    intercalate_n_merge.create_free_energy_dir(nres, resname, reswt_str, resmut_str, structure1, structure2)