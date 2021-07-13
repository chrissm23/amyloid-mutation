from biopandas.pdb import PandasPdb

def get_res_resname(windows_res_file):
    """Get residue number and residue name of the mutant"""
    with open(windows_res_file, 'r') as f:
        l = [x.split('=') for x in f]
        rows = [list(x) for x in zip(*l)]

    nres = int(rows[1][3])
    resname = rows[1][4][:-1]
    boundary_chains = rows[1][5][:-1]

    return [nres, resname, boundary_chains]

def create_mutants(amyloid_WT):
    """Create PDBs of the mutants separating according to the amount of mutants required"""
    nchains_p_mutant = [int((i+1)*amyloid_WT.nchains/amyloid_WT.n_mutants) for i in list(range(amyloid_WT.n_mutants))]
    chain_ids = amyloid_WT.biopandas_pdb.df['ATOM']['chain_id'].unique()
    mutant = amyloid_WT.biopandas_pdb
    for j in list(range(amyloid_WT.nchains)):
        mutant.df['ATOM']['residue_name'].loc[
            (mutant.df['ATOM']['chain_id'] == chain_ids[j]) & (mutant.df['ATOM']['residue_number'] == amyloid_WT.nres)] = amyloid_WT.resname
        mutant.df['ATOM'] = mutant.df['ATOM'].drop(
            mutant.df['ATOM'][(mutant.df['ATOM']['chain_id'] == chain_ids[j]) & (mutant.df['ATOM']['residue_number'] == amyloid_WT.nres) &
            (mutant.df['ATOM']['atom_name']!="CA") & (mutant.df['ATOM']['atom_name']!="C") & (mutant.df['ATOM']['atom_name']!="O") &
            (mutant.df['ATOM']['atom_name']!="N")].index)
        if j+1 in nchains_p_mutant:
            mutant_idx = nchains_p_mutant.index(j+1)
            mutant.to_pdb(path=f'./leap/mutant_{mutant_idx + 1}.pdb', records=None, gz=False, append_newline=True)
    
    mutants_pdb = [PandasPdb().read_pdb(f'./leap/mutant_{x+1}.pdb') for x in list(range(amyloid_WT.n_mutants))]
    return mutants_pdb, nchains_p_mutant

def num_to_chains(nchains_p_mutant):
    """Change chains per mutant to list of mutated chains"""
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
    
    mutated_chains = [[num_to_alph[j+1] for j in range(x)] for x in nchains_p_mutant]
    return mutated_chains

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
    windows_res_file = '../windows_res'
    [nres, resname, boundary_chains] = get_res_resname(windows_res_file)
    boundary_chains1 = boundary_chains.split(',')
    print(boundary_chains1)