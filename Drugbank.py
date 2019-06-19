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
#import RDkit as rd


def download_drugbank(NAME, PASSWORD, DRUGBANK_PATH, release='5-1-3'):
    """Download the whole drugbank database, it usually updates from every 2 weeks to 4 months
    Input -- NAME of the user and PASSWORD
    DRUGBANK_PATH -- where to store
    !! Carefully use it: the whole database is 130 MB packed, 1.3 GB unpacked
    """
    URL = 'https://www.drugbank.ca/releases/'+ release + '/downloads/all-full-database'
    make_dir(DRUGBANK_PATH)
    OUT_FILE = DRUGBANK_PATH / (URL.split('/')[-1] +'.zip')
    # Get file
    subprocess.check_output(['curl', '-Lfv', '-o', str(OUT_FILE), 
                             '-u', NAME + ':' + PASSWORD, URL])
    print('Downloading Drugbank')
    # Unpack file to the same directory
    with zipfile.ZipFile(str(OUT_FILE), 'r') as zip_ref:
            print("Extracting " + str(OUT_FILE))
            zip_ref.extractall(str(DRUGBANK_PATH))
    # Delete downloaded .zip
    subprocess.check_output(['rm', OUT_FILE])


def make_dir(dir_path):
    """Make directory with absolute path dir_name recursively"""
    if not os.path.exists(dir_path):
        os.makedirs(dir_path, exist_ok=True)
        print("Directory " , dir_path ,  " Created ")
    else:
        pass
    
    
def make_dir_from_list(dirList):
    for dirName in dirList:
        if not os.path.exists(dirName):
            os.mkdir(dirName)
            print("Directory ", dirName,  " created ")
        else:
            pass
        

def get_name_seq_from_fasta_lines(fasta):
    """Returns name and a/a sequence from fasta as a string
    INPUT - fasta as a string (with \n)
    OUTPUT - list (name, seq)
    """
    fasta_splitted = fasta.split('\n')
    name = fasta_splitted[0].split('|')[-1]
    seq = ''.join(fasta_splitted[1:])
    return name, seq
        

def dump_info_db(root):
    """Save all collected from Drugbank data as json files to root/Drugbank_exracted"""
    # All names of files to be dumped to root/Drugbank_extracted with name.txt, where name is from names
    names = ['ligands_unii', 'ligands_drugbank_ids', 'ligands_names', 'ligands_ids', 'ligands_resources',
             'ligands_ids_by_names', 'ligands_resources_by_names',
             'ligands_names_and_their_targets_ids', 'ligands_names_and_their_targets_resources',
             'targets_ids', 'targets_resources', 'targets_names',
             'ligands_smiles', 'ligands_names_and_smiles',
             'approved_flags', 'ligands_and_approved_flags',
             'targets_fastas', 'targets_names_and_fastas',
            ]
    name_full = str(Path(root) / 'Drugbank_extracted')
    make_dir_from_list([name_full])    
    for name in names:
        exec('global ' + name)
    for name in names:
        with open(str(Path(name_full) / (name + ".txt")), 'w') as f:
            exec('json.dump(' + name + ', f, ensure_ascii=False)')

    
def load_info_db(root):
    """Load all collected from Drugbank data as json files from root/Drugbank_exracted"""
    # All names of files to be loaded from root/Drugbank_extracted with name.txt, where name is from names
    names = ['ligands_unii', 'ligands_drugbank_ids', 'ligands_names', 'ligands_ids', 'ligands_resources',
             'ligands_ids_by_names', 'ligands_resources_by_names',
             'ligands_names_and_their_targets_ids', 'ligands_names_and_their_targets_resources',
             'targets_ids', 'targets_resources', 'targets_names',
             'ligands_smiles', 'ligands_names_and_smiles',
             'approved_flags', 'ligands_and_approved_flags',
             'targets_fastas', 'targets_names_and_fastas',
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
    source = str(Path(root) / 'Drugbank_extracted' / name)
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
    global approved_flags, ligands_and_approved_flags
    # Forligands and their targets
    global ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources
    global ligands_smiles, ligands_names_and_smiles
    # For targets
    global targets_ids, targets_resources, targets_names
    global targets_fastas, targets_names_and_fastas
    
    
    ligands_unii = []  # List of ligands' UNII ids
    ligands_drugbank_ids = []  # List of lists of Drugbank ids (one ligand could have several)
    ligands_names = []  # List of usual names
    ligands_ids = []  # List of lists of ids in different DBs
    ligands_resources = []  # List of lists of resources in different DBs
    ligands_smiles = []  # List of all SMILES of ligands
    approved_flags = []  # List of True/False of approvance of drug by FDA
    smiles = None
    
    targets_ids = []  # List of lists of lists of ids in different DBs
    targets_resources = []  # List of lists of lists of resources in different DBs
    targets_names = []  # List of lists of names of all targets
    targets_fastas = []  # List of lists of fastas
    unii = ""
    
    # Database iteration (DB is too big to parse it directly)
    for event, elem in context:
        # Iterate over drugs
        if event == "end" and db_tag(elem,"drug"):
            # Will be true in the end if one of products is approved
            # The ligand will be regarded as approved
            fda_approved = False
            # Flag check that
            f_drug_entry = False
            # Flag of having SMILES
            f_smiles = False
            smiles = None
            
            # Refresh all temporal lists of ids and resources
            # for ligands
            l_ids = []
            l_resources = []
            l_db_ids = []
            # for targets 
            t_id = []  # Id of one target in one database
            t_ids = []  # List of ids of one target in all databases
            t_resource = []  # Resource of one target (name of one database) 
            t_resources = []  # List of resources of one target (names of all databases)
            t_name = []  # Name of one target 
            t_names = []  # List of names of all targets
            t_fasta = []  # Fasta of one target
            t_fastas = []  # List of fasta strings of targets

            # Iteration over all fields of one drug
            for item in list(elem):
                # Get basic info of drug
                if db_tag(item, "unii"):
                    unii = item.text
                    # Check that it's not part of other drug entry (mb name 'drug' too as an interacting drug)
                    f_drug_entry = True
                if db_tag(item, "drugbank-id"):
                    l_db_ids.append(item.text)
                if db_tag(item, "name"):
                    name = item.text                    
                    
                # Get SMILES
                if db_tag(item, 'calculated-properties'):
                    # Iterate over properties
                    for item1 in list(item):
                        # Iterate over info about properties
                        for item2 in list(item1):
                            if db_tag(item2, 'kind'):
                                if item2.text == 'SMILES':
                                    f_smiles = True
                            if db_tag(item2, 'value') and f_smiles:
                                f_smiles = False
                                smiles = item2.text
                            
                # Checking if drug is approved by FDA
                if db_tag(item, "products"):
                    # Iterate over products
                    for it1 in list(item):
                        for it in list(it1):
                            # Check whether approved product exists
                            if db_tag(it, "approved"):
                                if it.text == "true":
                                    fda_approved = True

                # Get identifiers of drug and their databases
                if db_tag(item, "external-identifiers"):   
                    for it1 in list(item):
                        for it in list(it1):
                            if db_tag(it, "resource"):
                                l_resources.append(it.text)
                            if db_tag(it, "identifier"):
                                l_ids.append(it.text)

                # Get identifiers of targets and their databases
                if db_tag(item, "targets"):
                    for it1 in list(item):
                        f_polypeptide = False
                        # Iterating over properties of one target
                        for it2 in list(it1):
                            # Get name of the target
                            if db_tag(it2, 'name'):
                                t_name = it2.text
                            # Get info only about polypeptide target
                            if db_tag(it2, "polypeptide"):
                                f_polypeptide = True
                                for it3 in list(it2):
                                    if db_tag(it3, "external-identifiers"):
                                        for it4 in list(it3):
                                            for it5 in list(it4):
                                                if db_tag(it5, "resource"):
                                                    t_resource.append(it5.text)
                                                if db_tag(it5, "identifier"):
                                                    t_id.append(it5.text)
                                    # Get a/a sequence
                                    if db_tag(it3, "amino-acid-sequence"):
                                        t_fasta = it3.text
                                            
                                    # Gather all ids and resources of one polypeptide target
                                if f_polypeptide:
                                    t_ids.append(t_id)
                                    t_id = []
                                    t_resources.append(t_resource)
                                    t_resource = []
                                    t_names.append(t_name)
                                    t_name = []
                                    t_fastas.append(t_fasta)
                                    t_fasta = []


            # Clear in order not to store the whole database in memory    
            root_tree.clear()

            # If it was really drug entry => add information        
            if f_drug_entry:
                approved_flags.append(fda_approved)
                ligands_names.append(name)
                ligands_unii.append(unii)
                ligands_ids.append(l_ids)
                ligands_drugbank_ids.append(l_db_ids)
                ligands_resources.append(l_resources)
                ligands_smiles.append(smiles)

                targets_ids.append(t_ids)
                targets_resources.append(t_resources)
                targets_names.append(t_names)
                targets_fastas.append(t_fastas)
    
    # Create some useful dictionaries
    ligands_ids_by_names = dict(zip(ligands_names, ligands_ids))
    ligands_resources_by_names = dict(zip(ligands_names, ligands_resources))
    ligands_names_and_their_targets_ids = dict(zip(ligands_names, targets_ids))
    ligands_names_and_their_targets_resources = dict(zip(ligands_names, targets_resources))
    ligands_names_and_smiles = dict(zip(ligands_names, ligands_smiles))
    ligands_and_approved_flags = dict(zip(ligands_names, approved_flags))
    # Make dictionary {name of target:a/a sequence} and write to file sequences
    list_names = []
    list_fastas = []
    # Create file where to save fastas
    with open(str(Path(root) / 'Drugbank_extracted' / 'Drugbank_targets.fasta'), "w+") as myfile:
        pass
    for l_targets in targets_fastas:
        for fasta in l_targets:
            name, seq = get_name_seq_from_fasta_lines(fasta)
            if name not in list_names:
                list_names.append(name)
                list_fastas.append(seq)
                with open(str(Path(root) / 'Drugbank_extracted' / 'Drugbank_targets.fasta'), "a+") as myfile:
                    myfile.write(fasta)
                    myfile.write('\n')
    targets_names_and_fastas = dict(zip(list_names, list_fastas))
    
    # Save obtained data        
    dump_info_db(root)
    

def add_smiles_from_pubchem(root, ligands_names_and_smiles):
    """Add SMILES of ligands which don't have SMILES in Drugbank, but have it in PubChem (~10 ligands)"""
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
    global approved_flags, ligands_and_approved_flags
    # For ligands and their targets
    global ligands_smiles, ligands_names_and_smiles
    global ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources
    # For targets
    global targets_ids, targets_resources, targets_names
    global targets_fastas, targets_names_and_fastas
    
    # Directory where all data placed
    #root = '/media/anton/b8150e49-6ff0-467b-ad66-40347e8bb188/anton/BACHELOR'
    root = '/home/anton_maximov/BACHELOR'
    
    # If needed to download new version of Drugbank
    DRUGBANK_PATH = Path(root) / 'Drugbank_extracted'
    #download_drugbank('maksimov.as@phystech.edu', 'drugsandbanks', DRUGBANK_PATH, release='5-1-3')
    
    # Processing if new information needed or wasn't dumped before
    process_drugbank(root)
    # Add SMILES of ligands which don't have one in Drugbank but have it in PubChem
    add_smiles_from_pubchem(root, ligands_names_and_smiles)
    
    # Start of user actions
    # Load data from .txts from root/Drugbank_extracted to program
    #load_info_db(root)
    
    # Some checks of work
    #print(ligands_ids_by_names)
    #name_lig = 'Acetazolamide'#'Methazolamide'  #'Acetaminophen' 'Acetazolamide' 
    #uniprot = 'P00918'
    #a = get_all_smiles_of_ligands(ligands_ids_by_names, ligands_resources_by_names)
    #print(db.get_pdbs_from_smiles(sm, 0.5))
    #print(db.get_common_pdbs_from_ligand_name_and_target_uniprot(name_lig, uniprot, -0.1,
    #                                                  ligands_ids_by_names, ligands_resources_by_names,
    #                       ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources))
    #print(db.get_common_pdbs_with_all_targets_of_ligand(name_lig, -0.1,
    #                                                  ligands_ids_by_names, ligands_resources_by_names,
    #                       ligands_names_and_their_targets_ids, ligands_names_and_their_targets_resources))
