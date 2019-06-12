"""Get needed info from Drugbank, save it. Possible to load info and update the database"""
import os
import json
import xml.etree.ElementTree as etree
#import xml
from pathlib import Path
import subprocess
import zipfile
import pubchempy

import DATABASES_SMILES as db

# Carefully use it: the whole database is 130 MB packed, 1.3 GB unpacked
def download_drugbank(DRUGBANK_PATH, release='5-1-3'):
    """Download the whole drugbank database, it usually updates from every 2 weeks to 4 months"""
    # Drugbank authentication data
    NAME = 'your name' #'your name'
    PASSWORD = 'your password' #'your password'

    #URL = 'https://www.drugbank.ca/releases/5-1-3/downloads/target-all-polypeptide-ids'
    URL = 'https://www.drugbank.ca/releases/'+ release + '/downloads/all-full-database'

    OUT_FILE = DRUGBANK_PATH / (URL.split('/')[-1] +'.zip')
    # Get file
    subprocess.check_output(['curl', '-Lfv', '-o', str(OUT_FILE), 
                             '-u', NAME + ':' + PASSWORD, URL])
    # Unpack file to the same directory
    with zipfile.ZipFile(str(OUT_FILE), 'r') as zip_ref:
            print("Extracting " + URL.split("/")[-1])
            zip_ref.extractall(str(DRUGBANK_PATH))
    # Delete dewnloaded .zip
    subprocess.check_output(['rm', OUT_FILE])  


def make_dir(dirList):
    for dirName in dirList:
        if not os.path.exists(dirName):
            os.mkdir(dirName)
            print("Directory " , dirName ,  " Created ")
        else:
            pass
        

def dump_info_db(root):
    """Save all collected from Drugbank data to root/Drugbank_exracted"""
    # All names of files to be dumped to root/Drugbank_extracted with name.txt, where name is from names
    names = ['ligands_unii', 'ligands_drugbank_ids', 'ligands_names', 'ligands_ids', 'ligands_resources',
             'ligands_ids_by_names', 'ligands_resources_by_names',
             'ligands_names_and_their_targets_ids', 'ligands_names_and_their_targets_resources',
             'targets_ids', 'targets_resources', 'targets_names',
             'ligands_smiles', 'ligands_names_and_smiles',
            ]
    name_full = str(Path(root) / 'Drugbank_extracted')
    make_dir([name_full])    
    for name in names:
        exec('global ' + name)
    for name in names:
        with open(str(Path(name_full) / (name + ".txt")), 'w') as f:
            exec('json.dump(' + name + ', f, ensure_ascii=False)')

    
def load_info_db(root):
    """Load all collected from Drugbank data from root/Drugbank_exracted"""
    # All names of files to be loaded from root/Drugbank_extracted with name.txt, where name is from names
    names = ['ligands_unii', 'ligands_drugbank_ids', 'ligands_names', 'ligands_ids', 'ligands_resources',
             'ligands_ids_by_names', 'ligands_resources_by_names',
             'ligands_names_and_their_targets_ids', 'ligands_names_and_their_targets_resources',
             'targets_ids', 'targets_resources', 'targets_names',
             'ligands_smiles', 'ligands_names_and_smiles',
            ]
    name_full = str(Path(root) / 'Drugbank_extracted')

    for name in names:
        with open(str(Path(name_full) / (name + ".txt")), 'r') as f:
            exec('global ' + name + '\n' + name + ' = json.load(f)')


def db_tag(element, string):
    """Check that tag of element in Drugbank == needed string"""
    return element.tag.split("{http://www.drugbank.ca}")[1] == string


def process_drugbank(root, name='full database.xml'):
    """Get needed info from the full Drugbank database, placed in root with name
    https://www.drugbank.ca/docs/drugbank.xsd -- scheme of the base
    """
    # Location of the database
    source = str(Path(root) / name)
    # Get an iterable
    context = etree.iterparse(source, events=("start", "end"))

    # Turn it into an iterator
    context = iter(context)

    # Get the root element, for Python 2 here should be: event, root_tree = context.next()
    event, root_tree = next(context)

    # Here go lists with collected information about approved by FDA ligands in Drugbank
    # !!! Maybe use some OOP instead of list of lists
    
    # Initialize variables for ligands
    global ligands_unii, ligands_drugbank_ids, ligands_names, ligands_ids, ligands_resources
    global ligands_ids_by_names, ligands_resources_by_names
    # Forligands and their targets
    global ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources
    global ligands_smiles, ligands_names_and_smiles
    # For targets
    global targets_ids, targets_resources, targets_names
    
    
    ligands_unii = []  # List of ligands' UNII ids
    ligands_drugbank_ids = []  # List of lists of Drugbank ids (one ligand could have several)
    ligands_names = []  # List of usual names
    ligands_ids = []  # List of lists of ids in different DBs
    ligands_resources = []  # List of lists of resources in different DBs
    #ligands_ids_by_names = []  # Dictionary name:list of ids
    #ligands_resources_by_names = []  # Dictionary name:list of resources
    ligands_smiles = []  #list of all SMILES of ligands
    smiles = None
    
    targets_ids = []  # List of lists of lists of ids in different DBs
    targets_resources = []  # List of lists of lists of resources in different DBs
    targets_names = []  # List lists of names of all targets
    unii = ""
    # Flag of FDA approval
    fda_approved = False
    # Flag of having SMILES
    f_smiles = False
    smiles = None
    # Database iteration (DB is too big to parse it directly)
   
    for event, elem in context:
        # Find drugs storage
        if event == "end" and db_tag(elem,"drug"):        
            # Refresh all temporal lists of ids and resources
            # For ligands
            l_ids = []
            l_resources = []
            l_db_ids = []
            # For targets 
            t_id = []  # Id of one target in one database
            t_ids = []  # List of ids of one target in all databases
            t_resource = []  # Resource of one target (name of one database) 
            t_resources = []  # List of resources of one target (names of all databases)
            t_name = []  # Name of one target
            t_names = []  # List of names of all targets

            # Flag of FDA approval
            fda_approved = True
            #first_iter = True
            # Iteration over all fields of one drug
            for item in list(elem):
                # Exit drug loop if it is not approved 
                if not fda_approved:
                    break

                # Get basic info of drug
                if db_tag(item, "unii"):
                    unii = item.text
                if db_tag(item, "drugbank-id"):
                    l_db_ids.append(item.text)
                if db_tag(item, "name"):
                    name = item.text
                    
                # Get SMILES
                if db_tag(item, 'calculated-properties'):
                    for item1 in list(item):
                        for item2 in list(item1):
                            #print(item2.tag.split('}')[1])
                            if db_tag(item2, 'kind'):
                                if item2.text == 'SMILES':
                                    f_smiles = True
                            if db_tag(item2, 'value') and f_smiles:
                                f_smiles = False
                                smiles = item2.text
                            
                # Checking if drug is approved by FDA
                if db_tag(item, "products"):
                    for it1 in list(item):
                        for it in list(it1):
                            if db_tag(it, "approved"):
                                if it.text == "true":
                                    fda_approved = True
                                else:
                                    fda_approved = False

                # Get identifiers of drug and their databases
                # (maybe, it's worth to rewrite the ladder with function)
                if db_tag(item, "external-identifiers"):   
                    for it1 in list(item):
                        for it in list(it1):
                            if db_tag(it, "resource"):
                                l_resources.append(it.text)
                            if db_tag(it, "identifier"):
                                l_ids.append(it.text)

                # Get identifiers of targets and their databases
                # ADD: extract names of targets?
                if db_tag(item, "targets"):
                    for it1 in list(item):
                        for it2 in list(it1):
                            if db_tag(it2, 'name'):
                                t_name = it2.text
                            if db_tag(it2, "polypeptide"):
                                for it3 in list(it2):
                                    if db_tag(it3, "external-identifiers"):
                                        for it4 in list(it3):
                                            for it5 in list(it4):
                                                if db_tag(it5, "resource"):
                                                    t_resource.append(it5.text)
                                                if db_tag(it5, "identifier"):
                                                    t_id.append(it5.text)
                                        # Gather all ids and resorces of one target
                                        t_ids.append(t_id)
                                        t_id = []
                                        t_resources.append(t_resource)
                                        t_resource = []
                                        t_names.append(t_name)
                                        t_name = []

            # In order not to store the whole database in memory    
            root_tree.clear()

        # If approved => add information        
        if unii != "" and fda_approved and l_ids:
            fda_approved = False
            ligands_names.append(name)
            ligands_unii.append(unii)
            ligands_ids.append(l_ids)
            ligands_drugbank_ids.append(l_db_ids)
            ligands_resources.append(l_resources)
            ligands_smiles.append(smiles)
            smiles = None

            targets_ids.append(t_ids)
            targets_resources.append(t_resources)
            targets_names.append(t_names)
    
    # Create some useful dictionaries
    ligands_ids_by_names = dict(zip(ligands_names, ligands_ids))
    ligands_resources_by_names = dict(zip(ligands_names, ligands_resources))
    ligands_names_and_their_targets_ids = dict(zip(ligands_names, targets_ids))
    ligands_names_and_their_targets_resources = dict(zip(ligands_names, targets_resources))
    ligands_names_and_smiles = dict(zip(ligands_names, ligands_smiles))
    
    # Save obtained data        
    dump_info_db(root)
    

def add_smiles_from_pubchem(root, ligands_names_and_smiles):
    """Add SMILES of ligands which don't have one in Drugbank but have it in PubChem (~10 ligands)"""
    # Needed lists with data
    # Iterating over names, finding ones without SMILES and trying to get SMILES from PubChem    
    for name in ligands_names:
        if not ligands_names_and_smiles[name]:
            try:
                ind_compound = ligands_resources_by_names[name].index('PubChem Compound')
                # Find SMILES
                pubchem = ligands_ids_by_names[name][ind_compound]
                # Get smiles from PubCHEM
                c = pubchempy.Compound.from_cid(pubchem)
                smiles = c.isomeric_smiles
                ligands_smiles[ligands_names.index(name)] = smiles
            except ValueError:
                pass
    # Save corrected data
    #global ligands_names_and_smiles
    ligands_names_and_smiles = dict(zip(ligands_names, ligands_smiles))
    name_full = str(Path(root) / 'Drugbank_extracted')
    with open(str(Path(name_full) / ('ligands_smiles' + '.txt')), 'w') as f:
            json.dump(ligands_smiles, f, ensure_ascii=False)
    with open(str(Path(name_full) / ('ligands_names_and_smiles' + ".txt")), 'w') as f:
            json.dump(ligands_names_and_smiles, f, ensure_ascii=False)


if __name__ == "__main__":
    
    # Initialize variables for ligands
    global ligands_unii, ligands_drugbank_ids, ligands_names, ligands_ids, ligands_resources
    global ligands_ids_by_names, ligands_resources_by_names
    # Forligands and their targets
    global ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources
    global ligands_smiles, ligands_names_and_smiles
    # For targets
    global targets_ids, targets_resources, targets_names
    
    # Directory where all data placed
    root = '/media/anton/b8150e49-6ff0-467b-ad66-40347e8bb188/anton/BACHELOR'
    ROOT_PATH = '/media/anton/b8150e49-6ff0-467b-ad66-40347e8bb188/anton/BACHELOR'
    DRUGBANK_PATH = Path(root) / 'Drugbank'

    # If needed to update:
    #update_drugbank(DRUGBANK_PATH)
    
    # Processing if new information needed or wasn't dumped before
    process_drugbank(root)
    # Add SMILES of ligands which don't have one in Drugbank but have it in PubChem
    add_smiles_from_pubchem(root, ligands_names_and_smiles)
    
    # Start of user actions
    # Load data from .txts from root/Drugbank_extracted to program
    load_info_db(root)
    
    # Some checks of work
    #print(ligands_ids_by_names)
    name_lig = 'Acetazolamide'#'Methazolamide'  #'Acetaminophen' 'Acetazolamide' 
    uniprot = 'P00918'
    #a = get_all_smiles_of_ligands(ligands_ids_by_names, ligands_resources_by_names)
    #print(db.get_pdbs_from_smiles(sm, 0.5))
    #print(db.get_common_pdbs_from_ligand_name_and_target_uniprot(name_lig, uniprot, -0.1,
    #                                                  ligands_ids_by_names, ligands_resources_by_names,
    #                       ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources))
    #print(db.get_common_pdbs_with_all_targets_of_ligand(name_lig, -0.1,
    #                                                  ligands_ids_by_names, ligands_resources_by_names,
    #                       ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources))
