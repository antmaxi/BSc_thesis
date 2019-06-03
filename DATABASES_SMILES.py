
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
    """Get list of .pdb which include this uniprot
    uniprot -- id in UNIPROT
    path_to_save -- directory-string where to save with name uniprot_pdbs.txt if needed
    """
    url = 'https://www.uniprot.org/uploadlists/?from=ID&to=PDB_ID&format=list&query=' \
        + uniprot
    # If only list is needed
    if path_to_save is not None:
        r = download_url(url)
    # If need to save the list
    else:
        path_list_of_pdbs = str(Path(path_to_save) / uniprot)
        name = uniprot + '_pdbs.txt'
        r = download_url(url, path_list_of_pdbs, name)
        open(os.path.join(path_list_of_pdbs, name), 'wb').write(r.encode('utf-8'))
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
    print(str(path))
    # Trying to find appropriate similarity level to find at least one structure just by descending 
    # (mb, implement binary search?)
    if step_or_exact < 0:
        step = step_or_exact
        sim = 1.0 + step
        pdbs_from_smiles = []
        while not pdbs_from_smiles:
            url = 'http://www.rcsb.org/pdb/rest/smilesQuery?smiles=' + smiles \
            + '&search_type=similarity&similarity=' + str(sim)
            pdbs_from_smiles = []
            sim -= step
            r = requests.get(url, allow_redirects=True)
            print("Pdbs list downloaded")
            a = open(os.path.join(path, str(name) + '.xml'), 'wb').write(r.content)
            tree = ET.parse(os.path.join(path, str(name) + '.xml'))
            root = tree.getroot()
            for child in root:
                for child1 in child:
                    pdbs_from_smiles.append(child1.attrib['structureId'])
            subprocess.check_output(['rm', os.path.join(path, str(name) + '.xml')]) 
    # Exact search
    else:
        sim = step_or_exact
        url = 'http://www.rcsb.org/pdb/rest/smilesQuery?smiles=' + smiles \
            + '&search_type=similarity&similarity=' + str(sim)
        pdbs_from_smiles = []
        r = requests.get(url, allow_redirects=True)
        a = open(os.path.join(path, str(name) + '.xml'), 'wb').write(r.content)

        tree = ET.parse(os.path.join(path, str(name) + '.xml'))
        root = tree.getroot()
        for child in root:
            for child1 in child:
                pdbs_from_smiles.append(child1.attrib['structureId'])
        subprocess.check_output(['rm', os.path.join(path, str(name) + '.xml')]) 
                
    return sim, pdbs_from_smiles


def get_smiles_from_name(name, ligands_ids_by_names, ligands_resources_by_names):
    """Return SMILES using usual name of drug in Drugbank
    ligands_ids_by_names -- Dictionary name:list of ids, made by Drugbank processing and dump
    ligands_resources_by_names -- Dictionary name:list of resources, made by Drugbank processing and dump 
    """
    ligands_ids_by_names = dict(zip(ligands_names, ligands_ids))
    ligands_resources_by_names = dict(zip(ligands_names, ligands_resources))
    ind = -1
    try:
        # Find index of PubChem id
        ind = ligands_resources_by_names[name].index('PubChem Substance')
        #Find id
        pubchem = ligands_ids_by_names['Cetuximab'][ind]
        # Get smiles from PubCHEM
        c = pubchempy.Compound.from_cid(pubchem)
        smiles = c.isomeric_smiles
        return smiles
    except:
        print(f'{name} doesn\'t have pubchem id')
        return -1
    

def get_targets_uniprots_by_ligand_name(name_lig, ligands_names_and_their_targets_resources,
                                        ligands_names_and_their_targets_ids):
    """Get list of target's uniprots by name of ligand"""
    targets = []
    for target_resources in ligands_names_and_their_targets_resources[name_lig]:
        for index, target_resource in enumerate(target_resources):
            if target_resource == 'UniProtKB':
                targets.append(ligands_names_and_their_targets_ids[name_lig][index])
    return targets

#Examples:
#get_pdbs_from_smiles('CCOC1=CC=C(C=C1)NS(=O)(=O)C2=CC(=NN2)C(=O)NC3=CC(=CC=C3)SC', \
#                     step_or_exact=0.45, name='46507011')
