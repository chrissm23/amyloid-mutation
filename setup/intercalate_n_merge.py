import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
from shutil import copyfile
import os
import subprocess

import get_data
import solvate

def recover_transformation_order(structure1, structure2):
    """Recover the order in which the alchemical transformation will be carried out"""
    if structure1[1] == "WT":
        amyloid1 = structure1
        amyloid2 = structure2
    elif structure2[1] == "WT":
        amyloid1 = structure2
        amyloid2 = structure1
    else:
        ordered_structures = [structure1[1], structure2[1]]
        ordered_structures.sort()
        if structure1[1] == ordered_structures[0]:
            amyloid1 = structure1
            amyloid2 = structure2
        elif structure2[1] == ordered_structures[0]:
            amyloid1 = structure2
            amyloid2 = structure1
    
    return [amyloid1, amyloid2]

def intercalate(structure1, structure2):
    """Create PDB with chains of structures intercalated"""
    dir = f'{structure1[1]}_{structure2[1]}'
    newdir = f'./{dir}/chains_to_intercalate'
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    num_to_alph = {
        1: "A",
        2: "B",
        3: "C",
        4: "D",
        5: "E",
        6: "F",
        7: "G",
        8: "H",
        9: "I",
        10: "J",
        11: "K",
        12: "L",
        13: "M",
        14: "N",
        15: "O",
        16: "P",
        17: "Q",
        18: "R",
        19: "S",
        20: "T",
        21: "U",
        22: "V",
        23: "W",
        24: "X",
        25: "Y",
        26: "Z",
        27: "AA",
        28: "AB",
        29: "AC",
        30: "AD",
        31: "AE",
        32: "AF"
    }

    [amyloid1, amyloid2] = recover_transformation_order(structure1, structure2)
    amyloid1_pdb = PandasPdb().read_pdb(f'./{dir}/{amyloid1[1]}_strip.pdb')
    amyloid2_pdb = PandasPdb().read_pdb(f'./{dir}/{amyloid2[1]}_strip.pdb')
    nchains = structure1[0].df['ATOM']['chain_id'].nunique()
    nres_pchain = structure1[0].df['ATOM']['residue_number'].nunique()
    nres_min1 = amyloid1_pdb.df['ATOM']['residue_number'].min()
    nres_min2 = amyloid2_pdb.df['ATOM']['residue_number'].min()
    chains_to_intercalate = ''
    chains = ''
    for i in range(nchains):      
        temp_chain1 = amyloid1_pdb
        temp_chain1.df['ATOM'] = temp_chain1.df['ATOM'][(temp_chain1.df['ATOM']['residue_number'] >= i*nres_pchain + nres_min1) & 
            (temp_chain1.df['ATOM']['residue_number'] < (i+1)*nres_pchain + nres_min1)]
        temp_chain1.df['OTHERS'] = temp_chain1.df['OTHERS'].iloc[[i+1,-1]]
        temp_chain1.to_pdb(f'{newdir}/chain1_{num_to_alph[i+1]}.pdb')
        amyloid1_pdb = PandasPdb().read_pdb(f'./{dir}/{amyloid1[1]}_strip.pdb')
        temp_chain2 = amyloid2_pdb
        temp_chain2.df['ATOM'] = temp_chain2.df['ATOM'][(temp_chain2.df['ATOM']['residue_number'] >= i*nres_pchain + nres_min2) & 
            (temp_chain2.df['ATOM']['residue_number'] < (i+1)*nres_pchain + nres_min2)]
        temp_chain2.df['OTHERS'] = temp_chain2.df['OTHERS'].iloc[[i+1,-1]]
        temp_chain2.to_pdb(f'{newdir}/chain2_{num_to_alph[i+1]}.pdb')
        amyloid2_pdb = PandasPdb().read_pdb(f'./{dir}/{amyloid2[1]}_strip.pdb')

        chains_to_intercalate += f'a1_{num_to_alph[i+1]} = loadpdb $basedir/chains_to_intercalate/chain1_{num_to_alph[i+1]}.pdb\n'
        chains_to_intercalate += f'a2_{num_to_alph[i+1]} = loadpdb $basedir/chains_to_intercalate/chain2_{num_to_alph[i+1]}.pdb\n'
        chains += f'a1_{num_to_alph[i+1]} a2_{num_to_alph[i+1]} '

    replace_dict = {
        "%dir%": dir,
        "%chains_to_intercalate%": chains_to_intercalate,
        "%chains%": chains
    }
    path_intercalate = f'./{dir}/intercalate.sh'
    copyfile('./solvate/intercalate.tmpl', path_intercalate)
    solvate.replace_in_file(path_intercalate, replace_dict)
    print("Intercalating and creating dual topology with ions...")
    solvate.make_executable(path_intercalate)
    subprocess.call(path_intercalate)
    print("Intercalation finished.")

def merge(nres, res_min, res_max, nchains, structure1, structure2):
    """Get masks to be used in tiMerge"""
    n_resids = res_max - res_min + 1
    nres_new = nres - res_min + 1
    chains_WT = [f'{i*n_resids + 1}-{(i+1)*n_resids}' for i in range(2*nchains) if i%2 == 0]
    chains_WT_str = ','.join(chains_WT)
    chains_MUT = [f'{i*n_resids + 1}-{(i+1)*n_resids}' for i in range(2*nchains) if (i+1)%2 == 0]
    chains_MUT_str = ','.join(chains_MUT)
    resids_WT = [f'{i*n_resids + nres_new}' for i in range(2*nchains) if i%2 == 0]
    resids_WT_str = ','.join(resids_WT)
    resids_MUT = [f'{i*n_resids + nres_new}' for i in range(2*nchains) if (i+1)%2 == 0]
    resids_MUT_str = ','.join(resids_MUT)

    dir = f'{structure1[1]}_{structure2[1]}'
    path_merge = f'./{dir}/merge.sh'
    copyfile('./solvate/merge.tmpl', path_merge)

    replace_dict = {
        "%dir%": dir,
        "%chains_WT%": chains_WT_str,
        "%chains_MUT%": chains_MUT_str,
        "%resids_WT%": resids_WT_str,
        "%resids_MUT%": resids_MUT_str
    }

    solvate.replace_in_file(path_merge, replace_dict)
    print("Merging topologies...")
    solvate.make_executable(path_merge)
    subprocess.call(path_merge)
    print("Merging finished.")

def in_files_setup(nres, resname, res_min, res_max, nchains, structure1, structure2):
    """Copies and fills out input files for solvation correction and minimization"""
    [amyloid1, amyloid2] = recover_transformation_order(structure1, structure2)
    charge1 = get_data.get_charge(amyloid1[0])
    charge2 = get_data.get_charge(amyloid2[0])

    n_resids = res_max - res_min + 1
    nres_new = nres - res_min + 1

    reswt = [f'{i*(n_resids + 1) + nres_new}' for i in range(nchains)]
    reswt_str = ','.join(reswt)
    resmut = [f'{(i+1)*(n_resids + 1)}' for i in range(nchains)]
    resmut_str = ','.join(resmut)

    dir = f'{structure1[1]}_{structure2[1]}'
    dual_topology = PandasPdb().read_pdb(f'./{dir}/dual_topology_ions.pdb')

    if charge1*charge2 > 0:
        if abs(charge1) < abs(charge2): #Get first charge2 - charge1 ions and add them to resmut_str
            d_charge = abs(charge2) - abs(charge1)
            positive_ions = 'Na+' in dual_topology.df['ATOM']['residue_name'].unique()
            if positive_ions:
                ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Na+']
                ions_masks = list(ions.loc[ions.index[[0, d_charge]], 'residue_number'])
                ions_masks = [str(x) for x in ions_masks]
                ions_masks = ','.join(ions_masks)
                ions_masks = ', ' + ions_masks
                resmut_str += ions_masks
            negative_ions = 'Cl-' in dual_topology.df['ATOM']['residue_name'].unique()
            if negative_ions:
                ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Cl-']
                ions_masks = list(ions.loc[ions.index[[0, d_charge]], 'residue_number'])
                ions_masks = [str(x) for x in ions_masks]
                ions_masks = ','.join(ions_masks)
                ions_masks = ', ' + ions_masks
                resmut_str += ions_masks
        if abs(charge1) > abs(charge2): #Get first charge1 - charge2 ions and add them to reswt_str
            d_charge = abs(charge1) - abs(charge2)
            positive_ions = 'Na+' in dual_topology.df['ATOM']['residue_name'].unique()
            if positive_ions:
                ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Na+']
                ions_masks = list(ions.loc[ions.index[[0, d_charge]], 'residue_number'])
                ions_masks = [str(x) for x in ions_masks]
                ions_masks = ','.join(ions_masks)
                ions_masks = ', ' + ions_masks
                reswt_str += ions_masks
            negative_ions = 'Cl-' in dual_topology.df['ATOM']['residue_name'].unique()
            if negative_ions:
                ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Cl-']
                ions_masks = list(ions.loc[ions.index[[0, d_charge]], 'residue_number'])
                ions_masks = [str(x) for x in ions_masks]
                ions_masks = ','.join(ions_masks)
                ions_masks = ', ' + ions_masks
                reswt_str += ions_masks
    elif charge1*charge2 < 0: #Get ions balancing charge of amyloid1 and add them to reswt_str. Get ions balancing charge amyloid2 and add them to resmut_str
        if charge1 < 0:
            ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Na+']
            ions_masks = list(ions.loc['residue_number'])
            ions_masks = [str(x) for x in ions_masks]
            ions_masks = ','.join(ions_masks)
            ions_masks = ', ' + ions_masks
            reswt_str += ions_masks
            ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Cl-']
            ions_masks = list(ions.loc['residue_number'])
            ions_masks = [str(x) for x in ions_masks]
            ions_masks = ','.join(ions_masks)
            ions_masks = ', ' + ions_masks
            resmut_str += ions_masks
        if charge1 > 0:
            ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Cl-']
            ions_masks = list(ions.loc['residue_number'])
            ions_masks = [str(x) for x in ions_masks]
            ions_masks = ','.join(ions_masks)
            ions_masks = ', ' + ions_masks
            reswt_str += ions_masks
            ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Na+']
            ions_masks = list(ions.loc['residue_number'])
            ions_masks = [str(x) for x in ions_masks]
            ions_masks = ','.join(ions_masks)
            ions_masks = ', ' + ions_masks
            resmut_str += ions_masks
    elif charge1*charge2 == 0:
        if charge1 < 0:
            ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Na+']
            ions_masks = list(ions.loc['residue_number'])
            ions_masks = [str(x) for x in ions_masks]
            ions_masks = ','.join(ions_masks)
            ions_masks = ', ' + ions_masks
            reswt_str += ions_masks
        if charge1 > 0:
            ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Cl-']
            ions_masks = list(ions.loc['residue_number'])
            ions_masks = [str(x) for x in ions_masks]
            ions_masks = ','.join(ions_masks)
            ions_masks = ', ' + ions_masks
            reswt_str += ions_masks
        if charge2 < 0:
            ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Na+']
            ions_masks = list(ions.loc['residue_number'])
            ions_masks = [str(x) for x in ions_masks]
            ions_masks = ','.join(ions_masks)
            ions_masks = ', ' + ions_masks
            resmut_str += ions_masks
        if charge2 > 0:
            ions = dual_topology.df['ATOM'][dual_topology.df['ATOM']['residue_name'] == 'Cl-']
            ions_masks = list(ions.loc['residue_number'])
            ions_masks = [str(x) for x in ions_masks]
            ions_masks = ','.join(ions_masks)
            ions_masks = ', ' + ions_masks
            resmut_str += ions_masks
        

    replace_dict = {
        "%reswt%": reswt_str,
        "%resmut%": resmut_str
    }
    replace_dict2 = {
        "%r%": str(nres),
        "%rn%": resname,
        "%pairID%": dir
    }
    
    newdir = f'./{dir}/solvation_correction'
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    for x in ['min', 'heat', 'press']:
        copyfile(f'./tmpls/{x}.tmpl', f'{newdir}/{x}.in')
        solvate.replace_in_file(f'{newdir}/{x}.in', replace_dict)
    copyfile('./solvation_correction/run_simulations.slurm', f'{newdir}/run_simulations.slurm')
    solvate.replace_in_file(f'{newdir}/run_simulations.slurm', replace_dict2)

    newdir = f'./{dir}/minimization'
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    copyfile('./tmpls/minimization.tmpl', f'{newdir}/minimization.in')
    solvate.replace_in_file(f'{newdir}/minimization.in', replace_dict)
    copyfile('./minimization/run_minimization.slurm', f'{newdir}/run_minimization.slurm')
    solvate.replace_in_file(f'{newdir}/run_minimization.slurm', replace_dict2)

    return [reswt_str, resmut_str]

def create_free_energy_dir(nres, resname, reswt_str, resmut_str, structure1, structure2):
    """Create directory for free energy simulations and copy templates from free_energy_tmpls"""
    dir = f'{structure1[1]}_{structure2[1]}'
    new_dir = f'./{dir}/free_energy'
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    src = "./free_energy_tmpls/"
    src_files = os.listdir(src)
    for file_name in src_files:
        full_file_name = os.path.join(src, file_name)
        copyfile(full_file_name, f'{new_dir}/{file_name}')
    
    replace_dict1 = {
        "%tmask1%": reswt_str,
        "%tmask2%": reswt_str,
        "%smask1%": resmut_str,
        "%smask2%": resmut_str
    }
    for x in ["heat", "prod"]:
        solvate.replace_in_file(f'{new_dir}/{x}.tmpl', replace_dict1)

    replace_dict2 = {
        "%r%": str(nres),
        "%rn%": resname,
        "%dir%": dir
    }
    for x in ["2_submit_equilibration.slurm", "4_submit_production.slurm"]:
        solvate.replace_in_file(f'{new_dir}/{x}', replace_dict2)
    
    for x in ["1_setup_equilibration.sh", "3_setup_simulations.sh"]:
        solvate.make_executable(f'{new_dir}/{x}')


if __name__ == "__main__":
    structure1 = [PandasPdb().read_pdb('./leap/5kk3.pdb'), "WT"]
    structure2 = [PandasPdb().read_pdb('./leap/mutant_1.pdb'), "mutant_1"]
    nres = 34
    resname = "VAL"
    reswt_str = "23,66,109,152,195,238,281,324,367,410,453,496,539,582,625,668,711,754"
    resmut_str = "43,86,129,172,215,258,301,344,387,430,473,516,559,602,645,688,731,774"
    create_free_energy_dir(nres, resname, reswt_str, resmut_str, structure1, structure2)
