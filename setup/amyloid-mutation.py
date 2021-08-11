from biopandas.pdb import PandasPdb

import amyloid
import get_data
import solvate
import intercalate_n_merge

# Get info about the mutation
[nres, resname, boundary_chains] = get_data.get_res_resname('../windows_res')
boundary_chains = boundary_chains.split(',')

# Get info about the amyloid structure
og_structure = PandasPdb().read_pdb('./leap/5kk3.pdb')
amyloid_WT = amyloid.Amyloid(og_structure, nres, resname, boundary_chains)

# Arrange how chains will be mutated and create mutants
mutants_pdb, nchains_p_mutant = get_data.create_mutants(amyloid_WT)
mutated_chains = get_data.num_to_chains(nchains_p_mutant)
mutants = [amyloid.Mutant(mutants_pdb[i], nres, resname, boundary_chains, mutated_chains[i]) for i in range(len(mutants_pdb))]
og_structure = PandasPdb().read_pdb('./leap/5kk3.pdb')
amyloid_WT = amyloid.Amyloid(og_structure, nres, resname, boundary_chains)

# Get order of solvation
charges = [x.charge for x in mutants]
charges = [amyloid_WT.charge] + charges
mutants = [amyloid_WT] + mutants
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

    # Intercalate and merge structures
    intercalate_n_merge.intercalate(structure1, structure2)
    intercalate_n_merge.merge(structure1, structure2)
    [reswt_str, resmut_str, tiwt_str, timut_str] = intercalate_n_merge.in_files_setup(structure1, structure2)
    
    intercalate_n_merge.create_free_energy_dir(reswt_str, resmut_str, tiwt_str, timut_str, structure1, structure2)