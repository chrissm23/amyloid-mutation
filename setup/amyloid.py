class Amyloid:
    """Class with data of amyloid to be used for analysis of mutations"""
    def __init__(self, biopandas_pdb, nres, resname, boundary_chains) -> None:
        self.biopandas_pdb = biopandas_pdb
        self.nres = nres
        self.resname = resname
        self.boundary_chains = boundary_chains
        self.res_min = None
        self.res_max = None
        self.nchains = None
        self.tiatoms_total = None
        self.n_mutants = None
        self.charge = None

        self.get_amyloid_info()
        self.get_tiatoms()
        self.get_charge()

        self.n_residues = self.res_max - self.res_min + 1
        self.nres_new = self.nres - self.res_min + 1

    def get_amyloid_info(self):
        """Get what residues constitute the amyloid structure and its number of chains"""
        self.res_min = self.biopandas_pdb.df['ATOM']['residue_number'].loc[self.biopandas_pdb.df['ATOM']['chain_id'] == 'A'].min()
        self.res_max = self.biopandas_pdb.df['ATOM']['residue_number'].loc[self.biopandas_pdb.df['ATOM']['chain_id'] == 'A'].max()
        self.nchains = self.biopandas_pdb.df['ATOM']['chain_id'].nunique()

    def get_tiatoms(self):
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
        tiatoms_mt = amino_acid_natoms[self.resname]
        aa_wt = self.biopandas_pdb.df['ATOM']['residue_name'].loc[self.biopandas_pdb.df['ATOM']['residue_number'] == self.nres].iloc[0]
        tiatoms_wt = amino_acid_natoms[aa_wt]
        tiatoms = tiatoms_wt + tiatoms_mt
        self.tiatoms_total = tiatoms*self.nchains
        self.n_mutants = self.tiatoms_total//500 + (self.tiatoms_total % 500 > 0)

    def get_charge(self):
        """Get the charge of the structure feeded"""
        charged_aa = {
        "ASP": -1,
        "GLU": -1,
        "CYS": 1,
        "LYS": 1,
        "HIP": 1
        }
        self.charge = 0
        chain_ids = self.biopandas_pdb.df['ATOM']['chain_id'].unique()
        res_nums = self.biopandas_pdb.df['ATOM']['residue_number'].unique()
        for x in chain_ids:
            for y in res_nums:
                res_seq = self.biopandas_pdb.df['ATOM'].groupby(['chain_id']).get_group(x).groupby(['residue_number']).get_group(y)['residue_name'].iloc[0]
                if res_seq in charged_aa:
                    self.charge += charged_aa[res_seq]

class Mutant(Amyloid):
    """Class with data of mutations of amyloid"""
    def __init__(self, biopandas_pdb, nres, resname, boundary_chains, mutated_chains) -> None:
        super().__init__(biopandas_pdb, nres, resname, boundary_chains)
        self.mutated_chains = mutated_chains