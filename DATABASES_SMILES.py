
"""Some functions for jumping between ids:
Uniprot => PDBs
PubCHEM => SMILES
SMILES => PDBs
Name of ligand => SMILES
"""

import os
import requests
import subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
import pubchempy

global ROOT_PATH
ROOT_PATH = '/media/anton/b8150e49-6ff0-467b-ad66-40347e8bb188/anton/BACHELOR'


def make_dir(dirList):
    for dirName in dirList:
        if not os.path.exists(dirName):
            os.mkdir(dirName)
            print("Directory " , dirName ,  " Created ")
        else:
            pass
        

def download_url(url, path=None, name=None):
    """Download from url to 'path/name', making path directory, if not existed"""
    r = requests.get(url, allow_redirects=True)
    if path:
        paths = []
        paths.append(path)
        make_dir(paths)
        open(os.path.join(paths[0], name), 'wb').write(r.content)
    return r.content.decode('utf-8')


def get_pdbs_from_uniprot(uniprot, path_to_save=None):
    """Return list of .pdb which include this uniprot
    uniprot -- id in UNIPROT
    path_to_save -- directory-string where to save with name uniprot_pdbs.txt if needed
    """
    url = 'https://www.uniprot.org/uploadlists/?from=ID&to=PDB_ID&format=list&query=' \
        + uniprot
    # If only list is needed
    if path_to_save is None:
        r = download_url(url)
    # If need to save the list
    else:
        path_list_of_pdbs = str(Path(path_to_save) / uniprot)
        name = uniprot + '_pdbs.txt'
        full_name = os.path.join(path_list_of_pdbs, name)
        if not full_name.is_file():
            r = download_url(url, path_list_of_pdbs, name)
            open(full_name, 'wb').write(r.encode('utf-8'))
    return ''.join(r).split()


def get_all_smiles_from_pubchem_ids(ligands_ids, ligands_resources):
    """ Get all SMILEs of compounds in the ligands_ids list
    https://pubchempy.readthedocs.io/en/latest/guide/install.html
    """
    pubchem_ids = []
    pubchem_smiles = []
    pubchem_numbers = []
    for ind1, lig_res in enumerate(ligands_resources):
        # !! Works too long if many requests
        if ind1 == 5: #-1 for normal work
            pass
            #break
        else:
            for ind2, resource in enumerate(lig_res):
                if resource == 'PubChem Substance':
                    try:
                        # Get correspondent PubCHEM id 
                        pubchem = ligands_ids[ind1][ind2]
                        # Get smiles from this id
                        c = pubchempy.Compound.from_cid(pubchem)
                        pubchem_smiles.append(c.isomeric_smiles)
                        # Add to list places of successfull ligands in whole list
                        pubchem_numbers.append(ind1)  
                        # Add to list of correspondent pubchem ids 
                        pubchem_ids.append(pubchem) 
                    except:
                        # !! What to do if no SMILES?
                        pass            
                    break
    return pubchem_smiles, pubchem_numbers, pubchem_ids              
                    
def get_pdbs_from_smiles(smiles, step_or_exact=-0.05, name='list_of_pdbs', root_path=ROOT_PATH):
    """Get list of pdbs containing similiar SMILES.
    if step_or_exact < 0 => searching for pdbs decreasing level of similiarity from 1.0 by |step_or_exact|
    if step_or_exact > 0 => searching with this similiarity level (from 0.0 to 1.0)
    if name != None then save list of pdbs in ROOT_PATH/SMILES in name.xml
    Output -- level of similiarity, list of pdb ids
    """
    path = Path(root_path) / 'SMILES'
    make_dir([str(path)])
    # Trying to find appropriate similarity level to find at least one structure just by descending 
    # (mb, implement binary search?)
    if step_or_exact < 0:
        print(step_or_exact)
        step = step_or_exact
        sim = 1.0 - step
        pdbs_from_smiles = []
        while not pdbs_from_smiles:
            sim += step
            print('sim=', sim)
            url = 'http://www.rcsb.org/pdb/rest/smilesQuery?smiles=' + smiles \
            + '&search_type=similarity&similarity=' + str(sim)
            pdbs_from_smiles = [] 
            r = requests.get(url, allow_redirects=True)
            a = open(os.path.join(path, str(name) + '.xml'), 'wb').write(r.content)
            # Parsing the result of request
            tree = ET.parse(os.path.join(path, str(name) + '.xml'))
            root = tree.getroot()
            for child in root:
                for child1 in child:
                    pdbs_from_smiles.append(child1.attrib['structureId'])
            #subprocess.check_output(['rm', os.path.join(path, str(name) + '.xml')]) 
    # Exact search
    else:
        print(step_or_exact)
        sim = step_or_exact
        url = 'http://www.rcsb.org/pdb/rest/smilesQuery?smiles=' + smiles \
            + '&search_type=similarity&similarity=' + str(sim)
        pdbs_from_smiles = []
        r = requests.get(url, allow_redirects=True)
        a = open(os.path.join(path, str(name) + '.xml'), 'wb').write(r.content)
        # Parsing the result of request
        tree = ET.parse(os.path.join(path, str(name) + '.xml'))
        root = tree.getroot()
        for child in root:
            for child1 in child:
                pdbs_from_smiles.append(child1.attrib['structureId'])
        #subprocess.check_output(['rm', os.path.join(path, str(name) + '.xml')]) 
                
    return sim, pdbs_from_smiles


def get_targets_uniprots_from_ligand_name(name_lig, ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources):
    """Get list of target's uniprots by name of ligand
    INPUT - name of ligand,
    dictionary name:list of ids of targets of this ligand,
    dictionary name:list of resources of targets of this ligand,
    """
    targets = []
    for index_target, target_resources in enumerate(ligands_names_and_their_targets_resources[name_lig]):
        for index_resource, target_resource in enumerate(target_resources):
            if target_resource == 'UniProtKB':
                targets.append(ligands_names_and_their_targets_ids[name_lig][index_target][index_resource])
    return targets



def get_smiles_from_name(name, ligands_ids_by_names, ligands_resources_by_names):
    """Return SMILES using usual name of drug in Drugbank
    ligands_ids_by_names -- Dictionary name:list of ids, made by Drugbank processing and dump
    ligands_resources_by_names -- Dictionary name:list of resources, made by Drugbank processing and dump 
    """
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
        
        
    

def get_common_pdbs_from_ligand_name_and_target_uniprot(name_lig, uniprot, sim,
                                                       ligands_ids_by_names, ligands_resources_by_names,
                            ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources):
    """Get list of pdbs which are in both lists of pdbs 
    INPUT -- ligand's name name_lig (SMILES searched by sim, see function get_smiles_from_name)
    uniprot of target
    OUTPUT -- [(name_lig, uniprot), [list of common pdbs]]
    """
    # Initializing needed global variables
    #global ligands_ids_by_names, ligands_resources_by_names
    #global ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources
    
    # Get pdbs lists for ligand and target
    smiles_of_ligand = get_smiles_from_name(name_lig, ligands_ids_by_names, ligands_resources_by_names)
    pdbs_of_ligand = get_pdbs_from_smiles(smiles_of_ligand, sim)
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
        

def get_common_pdbs_with_all_targets_of_ligand(name_lig, sim,
                                                       ligands_ids_by_names, ligands_resources_by_names,
                            ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources):
    """Get list of all pdbs connecting ligand and its target
    INPUT -- ligand's name is name_lig (SMILES searched by sim, see function get_smiles_from_name)
    OUTPUT -- dictionary {[lig_name, uniprot]:[list of common pdbs]}
    """
    # Initializing needed global variables
    #global ligands_ids_by_names, ligands_resources_by_names
    #global ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources
    
    # Get pdbs lists for ligand and target
    smiles_of_ligand = get_smiles_from_name(name_lig, ligands_ids_by_names, ligands_resources_by_names)
    pdbs_of_ligand = get_pdbs_from_smiles(smiles_of_ligand, sim)
    uniprots_of_targets = get_targets_uniprots_from_ligand_name(name_lig, ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources)
    
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
        

#Examples:
#get_pdbs_from_smiles('CCOC1=CC=C(C=C1)NS(=O)(=O)C2=CC(=NN2)C(=O)NC3=CC(=CC=C3)SC', \
#                     step_or_exact=0.45, name='46507011')
