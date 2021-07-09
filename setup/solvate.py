from biopandas.pdb import PandasPdb
import fileinput
from shutil import copyfile
import os
import subprocess

def make_executable(path):
    """Make file executable"""
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2 #copy R bits to X
    os.chmod(path, mode)

def replace_in_file(path, replace_dict):
    """Replace in file changes[0] to changes[1]"""
    with fileinput.FileInput(path, inplace=True, backup='.bak') as f:
        for line in f:
            new_line = line
            for change in replace_dict:
                new_line = new_line.replace(change, replace_dict[change])
            print(new_line, end='')

def solvate(structure1, structure2):
    """Delete hydrogens from structures and solvate dual topology using leap from solvate.tmpl"""
    for structure in [structure1, structure2]:
        structure[0].biopandas_pdb.df['ATOM'] = structure[0].biopandas_pdb.df['ATOM'][structure[0].biopandas_pdb.df['ATOM']['element_symbol']!="H"]
        structure[0].biopandas_pdb.to_pdb(path=f'./leap/{structure[1]}_nohyd.pdb', records=None, gz=False, append_newline=True)

    new_dir = f'{structure1[1]}_{structure2[1]}'
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    
    path_solvate = f'./{new_dir}/solvate.sh'
    copyfile('./solvate/solvate.tmpl', path_solvate)
    replace_dict = {
        "%amyloid1%": f"{structure1[1]}_nohyd",
        "%amyloid2%": f"{structure2[1]}_nohyd",
        "%dir%": new_dir
    }
    replace_in_file(path_solvate, replace_dict)
    print("Solvating dual topology...")
    make_executable(path_solvate)
    subprocess.call(path_solvate)
    print("Solvation finished.")

def strip(structure1, structure2):
    """Separate each amyloid structure and solvation water from dual topology using strip.tmpl"""
    dir = f'{structure1[1]}_{structure2[1]}'
    if not os.path.exists(dir):
        raise FileNotFoundError('Directory with that selection and order of structures does not exist.')

    nres_1 = structure1[0].n_residues*structure1[0].nchains
    nres_2 = structure2[0].n_residues*structure2[0].nchains
    nres_dual = nres_1 + nres_2

    path_strip = f'./{dir}/strip.sh'
    copyfile('./solvate/strip.tmpl', path_strip)
    replace_dict = {
        "%no_residues_dual%": str(nres_dual),
        "%no_residues_dual_pone%": str(nres_dual + 1),
        "%no_residues_single%": str(nres_1),
        "%no_residues_single_pone%": str(nres_1 + 1),
        "%structure_1%": structure1[1],
        "%structure_2%": structure2[1],
        "%dir%": dir
    }
    replace_in_file(path_strip, replace_dict)
    print("Stripping water and structures...")
    make_executable(path_strip)
    subprocess.call(path_strip)
    print("Strip finished.")

def solvate_ions_structure(structure, charge, path, second):
    """Solvate amyloid structure with specified water and final total charge"""
    if second == False:
        water = "water"
        one_or_two = "1"
    else:
        if not os.path.exists(f'./{path}/water_ions_1.pdb'):
            raise FileNotFoundError('water_ions_1.pdb from first solvation with ions not found.')
        water = "water_ions_1"
        one_or_two = "2"
    
    path_solvate_ions = f'./{path}/solvate_ions_{one_or_two}.sh'
    copyfile('./solvate/solvate_ions.tmpl', path_solvate_ions)
    replace_dict = {
        "%dir%": path,
        "%amyloid%": structure[1],
        "%water%": water,
        "%onetwo%": one_or_two,
        "%ch%": str(charge)
    }
    replace_in_file(path_solvate_ions, replace_dict)
    print("Solvating with ions...")
    make_executable(path_solvate_ions)
    subprocess.call(path_solvate_ions)
    print("Solvation with ions finished")

def strip_ions(structure, path, second):
    """Strip water with ions and amyloid with ions for either the first solvation with ions or the second one"""
    if second == False:
        if not os.path.exists(f'./{path}/single_topology_ions_1.pdb'):
            raise FileNotFoundError('single_topology_ions_1.pdb from first solvation with ions not found.')
        one_or_two = "1"
        single_topology = PandasPdb().read_pdb(f'./{path}/single_topology_ions_1.pdb')
        nres_pions = single_topology.df['ATOM'].loc[single_topology.df['ATOM']['residue_name']!='WAT']['residue_number'].max()
    else:
        if not os.path.exists(f'./{path}/single_topology_ions_2.pdb'):
            raise FileNotFoundError('single_topology_ions_2.pdb from second solvation with ions not found.')
        one_or_two = "2"
        single_topology = PandasPdb().read_pdb(f'./{path}/single_topology_ions_2.pdb')
        nres_pions = single_topology.df['ATOM'].loc[single_topology.df['ATOM']['residue_name']!='WAT']['residue_number'].max()
    
    nres = structure[0].n_residues*structure[0].nchains
    
    path_strip_ions = f'./{path}/strip_ions_{one_or_two}.sh'
    copyfile('./solvate/strip_ions.tmpl', path_strip_ions)
    replace_dict = {
        "%onetwo%": one_or_two,
        "%no_residues_single%": str(nres),
        "%no_residues_ions_pone%": str(nres_pions + 1),
        "%dir%": path
    }
    replace_in_file(path_strip_ions, replace_dict)
    print("Stripping with ions...")
    make_executable(path_strip_ions)
    subprocess.call(path_strip_ions)
    print("Strip ions finished")

def solvate_ions_pair(structure1, structure2):
    """Solvate pair of structures already ordered according to charge"""
    dir = f'{structure1[1]}_{structure2[1]}'
    if not os.path.exists(dir):
        raise FileNotFoundError('Directory with that selection and order of structures does not exist.')
    
    solvate_ions_structure(structure1, 0, dir, second=False)
    strip_ions(structure1, dir, second=False)
    charge1 = structure1[0].charge
    charge2 = structure2[0].charge
    if charge1*charge2 > 0:
        charge = 0
    else:
        charge = -charge1
    solvate_ions_structure(structure2, charge, dir, second=True)
    strip_ions(structure2, dir, second=True)