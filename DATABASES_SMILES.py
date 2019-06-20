
""" Some auxillary functions (download files and pdbs, load processed from Drugbank data) 
and functions for jumping between ids:
Uniprot => a/a sequence
Uniprot => PDBs
SMILES => PDBs
PubCHEM => SMILES
Name of ligand => Uniprots of targets
Name of ligand => SMILES

Also function to convert different structural types (pdb, sdf, mol2) between themselves
"""

import os
import requests
import json
import subprocess
from shutil import copyfile

from pathlib import Path
import xml.etree.ElementTree as ET

import pubchempy
import openbabel


def make_dir(dir_path):
    """ Make directory with absolute path dir_name recursively."""
    if not os.path.exists(dir_path):
        os.makedirs(dir_path, exist_ok=True)
        print("Directory " , dir_path ,  " Created ")
    else:
        pass
    
    
def make_dir_from_list(dirList):
    for dirName in dirList:
        if not os.path.exists(dirName):
            os.mkdir(dirName)
            print("Directory " , dirName ,  " Created ")
        else:
            pass
        

def download_url(url, directory=None, name=None, overwrite=False):
    """ Download from url to 'path/name', making path directory, if not existed."""
    # If needs to save
    if directory:
        # If file wasnt' downloaded or need to be overwrited
        if not (Path(directory) / name).is_file() or overwrite:
            r = requests.get(url, allow_redirects=True)
            make_dir(directory)
            open(os.path.join(directory, name), 'wb').write(r.content)
            return r.content.decode('utf-8')
    # If need only to return from function without saving
    else:
        r = requests.get(url, allow_redirects=True)
        return r.content.decode('utf-8')


def download_pdb(pdb, pdb_dir):
    """ Download pdb (e.g. 1A46) and save to pdb_path directory (with name 1A46.pdb).
    INPUT:
        pdb -- identifier in PDB
        pdb_dir -- directory path where to save
    """   
    path = Path(pdb_dir)
    name = pdb + ".pdb"
    url = "https://files.rcsb.org/download/" + pdb + ".pdb"
    full_name = path / name
    # Check if .pdb is already downloaded, if not => download
    if not full_name.is_file():
        download_url(url, str(path), name)

        
def load_info_db_from_namelist(namelist, root):
    """ Load variables listed in namelist from correspondent .txt file 
    in root/'Drugbank_exracted', collected from Drugbank data as json files 
    """
    # All names of files to be loaded from root/Drugbank_extracted with name.txt, where name is from names
    name_full = str(Path(root) / 'Drugbank_extracted')
    for name in namelist:
        with open(str(Path(name_full) / (name + ".txt")), 'r') as f:
            exec('global ' + name + '\n' + name + ' = json.load(f)')
            
            
def convert_single_structure(path_in, path_out):
    """ Convert single structure (or first in multi-structure) with full path path_in 
    and appropriate extension (e.g. 'sdf', 'mol2', 'pdb') to full path path_out with appropriate extension
    """
    obConversion = openbabel.OBConversion()
    # Set up conversion
    # Get extensions of files
    in_format = str(os.path.splitext(path_in)[1])
    out_format = str(os.path.splitext(path_out)[1])
    # Check if formats are the same
    if in_format == out_format:
        if path_in != path_out:
            copyfile(path_in, path_out)
    # Initialize conversion formats
    conversion_initialized = obConversion.SetInAndOutFormats(in_format[1:], out_format[1:])  # Formats without dot
    if not conversion_initialized:
        print('Conversion couldn\'t be done; maybe, file format is/are not appropriate')
        return -1
    mol = openbabel.OBMol()
    # Read file
    readed_normally = obConversion.ReadFile(mol, path_in)
    if not readed_normally:
        print('File wasn\'t read properly; maybe, not existed')
        return -1
    # Add hydrogens if needed
    #mol.AddHydrogens()
    # Write new file
    written_normally = obConversion.WriteFile(mol, path_out)
    if not written_normally:
        print('Writing was done with error')
        return -1
    
    
def get_seq_from_uniprot(uniprot):
    """Return a/a sequence of protein from uniprot"""
    url = 'https://www.uniprot.org/uniprot/?query=id:' + uniprot + '&format=fasta&columns=sequence'
    r = download_url(url)
    if not r:
        print(f'No sequence found, probably invalid uniprot {uniprot}')
        return -1
    else:
        return ''.join(r.split('\n')[1:])

    
def get_pdbs_from_uniprot(uniprot, path_to_save=None):
    """Return list of .pdb which include this uniprot
    INPUT:
        uniprot -- id in UNIPROT
        path_to_save -- directory path where to save with name uniprot_pdbs.txt (if signed by path_to_save)
    """
    url = 'https://www.uniprot.org/uploadlists/?from=ID&to=PDB_ID&format=list&query=' \
        + uniprot
    # If only list is needed
    if path_to_save is None:
        r = download_url(url)
        if not r:
            print(f'No pdbs found, probably invalid uniprot {uniprot}')
        return ''.join(r).split()
    # If need to save the list
    else:
        path_list_of_pdbs = str(Path(path_to_save) / uniprot)
        name = uniprot + '_pdbs.txt'
        full_name = os.path.join(path_list_of_pdbs, name)\
        # If file not exist => download
        if not Path(full_name).is_file():
            r = download_url(url, path_list_of_pdbs, name)
            if not r:
                print(f'No pdbs found, probably invalid uniprot {uniprot}')
            with open(full_name, 'wb') as f:
                f.write(r.encode('utf-8'))
        # If already downloaded => read
        else:
            with open(full_name, 'r') as f:
                r = f.read()
        return ''.join(r).split()
             
                    
def get_pdbs_from_smiles(smiles, root, step_or_exact=-0.05, name='list_of_pdbs'):
    """ Get list of pdbs containing similiar SMILES.
    if step_or_exact < 0 => searching for pdbs decreasing level of similiarity from 1.0 by |step_or_exact|
    if step_or_exact > 0 => searching with this similiarity level (from 0.0 to 1.0)
    if name != None then save list of pdbs in root in name.xml
    Output -- level of similiarity, list of pdb ids
    """
    path = Path(root)# / 'SMILES'
    #make_dir_from_list([str(path)])
    # Trying to find appropriate similarity level to find at least one structure just by descending 
    if step_or_exact < 0:
        print(step_or_exact)
        step = step_or_exact
        sim = 1.0 - step
        pdbs_from_smiles = []
        while not pdbs_from_smiles:
            sim += step
            #print('sim=', sim)
            url = 'http://www.rcsb.org/pdb/rest/smilesQuery?smiles=' + smiles \
            + '&search_type=similarity&similarity=' + str(sim)
            pdbs_from_smiles = [] 
            r = requests.get(url, allow_redirects=True)
            with open(os.path.join(path, str(name) + '.xml'), 'wb') as file:
                file.write(r.content)
            # Parsing the result of request
            tree = ET.parse(os.path.join(path, str(name) + '.xml'))
            root = tree.getroot()
            for child in root:
                for child1 in child:
                    pdbs_from_smiles.append(child1.attrib['structureId'])
            subprocess.check_output(['rm', os.path.join(path, str(name) + '.xml')]) 
    # Exact search
    else:
        sim = step_or_exact
        #print(f'Similarity level {step_or_exact}')
        url = 'http://www.rcsb.org/pdb/rest/smilesQuery?smiles=' + smiles \
            + '&search_type=similarity&similarity=' + str(sim)
        pdbs_from_smiles = []
        r = requests.get(url, allow_redirects=True)
        with open(os.path.join(path, str(name) + '.xml'), 'wb') as file:
            file.write(r.content)
        # Parsing the result of request
        tree = ET.parse(os.path.join(path, str(name) + '.xml'))
        root = tree.getroot()
        for child in root:
            for child1 in child:
                pdbs_from_smiles.append(child1.attrib['structureId'])
        subprocess.check_output(['rm', os.path.join(path, str(name) + '.xml')]) 
    return sim, pdbs_from_smiles


def get_targets_uniprots_from_ligand_name(name_lig, root):
    """ Get list of target's uniprots by name of ligand.
    INPUT - usual name of ligand in Drugbank
    """
    # Load needed info
    namelist = ['ligands_names_and_their_targets_ids', 
                'ligands_names_and_their_targets_resources'
               ]
    load_info_db_from_namelist(namelist, root)
    targets = []
    for index_target, target_resources in enumerate(ligands_names_and_their_targets_resources[name_lig]):
        for index_resource, target_resource in enumerate(target_resources):
            if target_resource == 'UniProtKB':
                targets.append(ligands_names_and_their_targets_ids[name_lig][index_target][index_resource])
    return targets


def get_smiles_from_name_from_pubchem(name, root):
    """ Return SMILES from Pubchem using usual name of drug in Drugbank."""
    # Load needed info
    namelist = ['ligands_ids_by_names', 
                'ligands_resources_by_names'
               ]
    load_info_db_from_namelist(namelist, root)
    # Find indexes of PubChem as a compound
    try:
        ind_compound = ligands_resources_by_names[name].index('PubChem Compound')
        # Find SMILES
        pubchem = ligands_ids_by_names[name][ind_compound]
        # Get smiles from PubCHEM
        c = pubchempy.Compound.from_cid(pubchem)
        smiles = c.isomeric_smiles
        return smiles
    except ValueError:
        print(f'{name} doesn\'t have pubchem id')
        return None


def get_common_pdbs_from_ligand_name_and_target_uniprot(name_lig, uniprot, sim, root):
    """ Get list of pdbs which contain both ligand and target.
    INPUT: 
        name_lig - ligand's usual name in Drugbank
        uniprot - Uniprot ID of target
        sim - level of similarity to search for pdbs from SMILES (see function get_pdbs_from_smiles)
        root - root of the protocol
    OUTPUT -- [(name_lig, uniprot), [list of common pdbs of ligand and target]]
    """
    # Load needed info
    namelist = ['ligands_ids_by_names', 
                'ligands_resources_by_names', 
                'ligands_names_and_their_targets_ids',
                'ligands_names_and_their_targets_resources'
               ]
    load_info_db_from_namelist(namelist, root)
    
    # Get pdbs lists for ligand and target
    smiles_of_ligand = get_smiles_from_name_from_pubchem(name_lig, root)
    pdbs_of_ligand = get_pdbs_from_smiles(smiles_of_ligand, root, sim)
    pdbs_of_target = get_pdbs_from_uniprot(uniprot)
    
    # Get list of common pdbs
    common_pdbs = list(set(pdbs_of_target) & set(pdbs_of_ligand[1]))
    if common_pdbs:
        # Produce list of keys [name_lig, uniprot]
        keys = (name_lig, uniprot)
        # Append list of pdbs to list of keys
        result = []
        result.append(keys)
        result.append(common_pdbs)
        return result
    else:
        return None
        

def get_common_pdbs_with_all_targets_of_ligand(name_lig, sim, root):
    """ Get list of all pdbs connecting ligand and its target.
    INPUT:
        name_lig - ligand's usual name in Drugbank
        sim - level of similarity to search for pdbs from SMILES (see function get_pdbs_from_smiles)
        root - root of the protocol
    OUTPUT -- dictionary {[lig_name, uniprot]:[list of common pdbs]}
    """
    # Load needed info
    namelist = ['ligands_ids_by_names', 'ligands_resources_by_names',
                'ligands_names_and_their_targets_ids', 'ligands_names_and_their_targets_resources']
    load_info_db_from_namelist(namelist, root)
    
    # Get pdbs lists for ligand and target
    smiles_of_ligand = get_smiles_from_name_from_pubchem(name_lig, root)
    pdbs_of_ligand = get_pdbs_from_smiles(smiles_of_ligand, root, sim)
    uniprots_of_targets = get_targets_uniprots_from_ligand_name(name_lig, root)
    
    # List of lists of all common pdbs
    all_common_pdbs = []
    # Keys of the future dict {[lig_name, uniprot]:[list of common pdbs]}
    keys_ligname_uniprot = []
    
    # Get common pdbs and create keys
    for uniprot in uniprots_of_targets:
        pdbs_of_target = get_pdbs_from_uniprot(uniprot)
        common_pdbs = list(set(pdbs_of_target) & set(pdbs_of_ligand[1]))
        all_common_pdbs.append(common_pdbs)
        keys_ligname_uniprot.append((name_lig, uniprot))
                  
    # Final dictionary
    result = dict(zip(keys_ligname_uniprot, all_common_pdbs))
    return result
        
root = '/home/anton_maximov/BACHELOR'
#Examples:
#get_pdbs_from_smiles('CCOC1=CC=C(C=C1)NS(=O)(=O)C2=CC(=NN2)C(=O)NC3=CC(=CC=C3)SC', root, \
#                     step_or_exact=0.7, name='46507011')
#get_pdbs_from_uniprot('P00533')
