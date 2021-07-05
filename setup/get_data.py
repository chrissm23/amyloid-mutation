from biopandas.pdb import PandasPdb

def get_res_resname(windows_res_file):
    """Get residue number and residue name of the mutant"""
    with open(windows_res_file, 'r') as f:
        l = [x.split('=') for x in f]
        rows = [list(x) for x in zip(*l)]

    nres = int(rows[1][3])
    resname = rows[1][4][:-1]

    return [nres, resname]

def get_amyloid_info(og_structure):
    """Get what residues constitute the amyloid structure and its number of chains"""
    res_min = og_structure.df['ATOM']['residue_number'].loc[og_structure.df['ATOM']['chain_id'] == 'A'].min()
    res_max = og_structure.df['ATOM']['residue_number'].loc[og_structure.df['ATOM']['chain_id'] == 'A'].max()
    nchains = og_structure.df['ATOM']['chain_id'].nunique()

    return [res_min, res_max, nchains]

def get_tiatoms(og_structure, nres, resname):
    """Get number of ti atoms per chain"""
    amino_acid_natoms = {
        "ALA": 13,
        "ARG": 26,
        "ASN": 17,
        "ASP": 16,
        "CYS": 14,
        "GLN": 20,
        "GLU": 19,
        "GLY": 10,
        "HIS": 20,
        "ILE": 22,
        "LEU": 22,
        "LYS": 24,
        "MET": 20,
        "PHE": 23,
        "PRO": 17,
        "SER": 14,
        "THR": 17,
        "TRP": 27,
        "TYR": 24,
        "VAL": 19,
        "HIP": 20
    }
    tiatoms_mt = amino_acid_natoms[resname]
    aa_wt = og_structure.df['ATOM']['residue_name'].loc[og_structure.df['ATOM']['residue_number'] == nres].iloc[0]
    tiatoms_wt = amino_acid_natoms[aa_wt]
    tiatoms = tiatoms_wt + tiatoms_mt

    return tiatoms

def create_mutants(og_structure, n_mutants, nchains, nres, resname):
    """Create PDBs of the mutants separating according to the amount of mutants required"""
    nchains_p_mutant = [int((i+1)*nchains/n_mutants) for i in list(range(n_mutants))]
    chain_ids = og_structure.df['ATOM']['chain_id'].unique()
    mutant = og_structure
    for j in list(range(nchains)):
        mutant.df['ATOM']['residue_name'].loc[
            (mutant.df['ATOM']['chain_id'] == chain_ids[j]) & (mutant.df['ATOM']['residue_number'] == nres)] = resname
        mutant.df['ATOM'] = mutant.df['ATOM'].drop(
            mutant.df['ATOM'][(mutant.df['ATOM']['chain_id'] == chain_ids[j]) & (mutant.df['ATOM']['residue_number'] == nres) &
            (mutant.df['ATOM']['atom_name']!="CA") & (mutant.df['ATOM']['atom_name']!="C") & (mutant.df['ATOM']['atom_name']!="O") &
            (mutant.df['ATOM']['atom_name']!="N")].index)
        if j+1 in nchains_p_mutant:
            mutant_idx = nchains_p_mutant.index(j+1)
            mutant.to_pdb(path=f'./leap/mutant_{mutant_idx + 1}.pdb', records=None, gz=False, append_newline=True)
    
    mutants = [PandasPdb().read_pdb(f'./leap/mutant_{x+1}.pdb') for x in list(range(n_mutants))]
    return mutants

def get_charge(pdb_structure):
    """Get the charge of the structure feeded"""
    charged_aa = {
    "ASP": -1,
    "GLU": -1,
    "CYS": 1,
    "LYS": 1,
    "HIP": 1
    }
    charge = 0
    chain_ids = pdb_structure.df['ATOM']['chain_id'].unique()
    res_nums = pdb_structure.df['ATOM']['residue_number'].unique()
    for x in chain_ids:
        for y in res_nums:
            res_seq = pdb_structure.df['ATOM'].groupby(['chain_id']).get_group(x).groupby(['residue_number']).get_group(y)['residue_name'].iloc[0]
            if res_seq in charged_aa:
                charge += charged_aa[res_seq]
    
    return charge

def get_order(charges):
    order = []
    for i in list(range(len(charges)-1)):
        if charges[i] != charges[i+1]:
            if charges[i]*charges[i+1] >= 0:
                if abs(charges[i]) < abs(charges[i+1]):
                    order.append(i)
                else:
                    order.append(i+1)
            elif charges[i]*charges[i+1] < 0:
                if i == 0:
                    order.append(i)
                else:
                    if charges[order[i-1]]*charges[i] >= 0:
                        order.append(i)
                    elif charges[order[i-1]]*charges[i+1] >= 0:
                        order.append(i+1)
        else:
            order.append(i)
    return order


if __name__ == "__main__":
    [nres, resname] = get_res_resname('../windows_res')
    og_structure = PandasPdb().read_pdb('./leap/5kk3.pdb')
    chains_ids = og_structure.df['ATOM']['chain_id'].unique()
    res_nums = og_structure.df['ATOM']['residue_number'].unique()
    res_seq = og_structure.df['ATOM'].groupby(['chain_id']).get_group(chains_ids[0]).groupby(['residue_number']).get_group(res_nums[0])['residue_name'].iloc[0]
    print(res_seq)