"""Get needed info from Drugbank, save it. Possible to load info and update the database"""

import os
import json
import xml.etree.ElementTree as etree
#import xml
from pathlib import Path
import subprocess
import zipfile


# Carefully use it: the whole database is 130 MB packed, 1.3 GB unpacked
def update_drugbank(DRUGBANK_PATH, release='5-1-3'):
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
    names = ['ligands_names', 'ligands_unii', 'ligands_ids', 'ligands_drugbank_ids', 'ligands_resources',
             'ligands_ids_by_names', 'ligands_resources_by_names',
             'targets_ids', 'targets_resources',
            ]
    name_full = str(Path(root) / 'Drugbank_extracted')
    make_dir([name_full])
    ligands_ids_by_names = dict(zip(ligands_names, ligands_ids))
    ligands_resources_by_names = dict(zip(ligands_names, ligands_resources))
    for name in names:
        with open(str(Path(name_full) / (name + ".txt")), 'w') as f:
            exec('global ' + name + '\n' + 'json.dump(' + name + ', f, ensure_ascii=False)')


def load_info_db(root):
    """Load all collected from Drugbank data from root/Drugbank_exracted"""
    names = ['ligands_names', 'ligands_unii', 'ligands_ids', 'ligands_drugbank_ids', 'ligands_resources',
             'ligands_ids_by_names', 'ligands_resources_by_names',
             'targets_ids', 'targets_resources'
            ]
    name_full = str(Path(root) / 'Drugbank_extracted')

    for name in names:
        with open(str(Path(name_full) / (name + ".txt")), 'r') as f:
            exec('global ' + name + '\n' + name + ' = json.load(f)')


def db_tag(element, string):
    """Check that tag of element in Drugbank == needed string"""
    return element.tag.split("{http://www.drugbank.ca}")[1] == string


def process_drugbank(root, name='full database.xml'):
    """Get needed info from the full Drugbank database, placed in root with name"""
    # !!! Check the first doubled entry
    source = str(Path(root) / name)
    # Get an iterable
    context = etree.iterparse(source, events=("start", "end"))

    # Turn it into an iterator
    context = iter(context)

    # Get the root element, for Python 2 here should be: event, root_tree = context.next()
    event, root_tree = next(context)

    # Here go lists with collected information about approved by FDA ligands in Drugbank
    # !!! Maybe use some OOP instead of list of lists
    global ligands_unii, ligands_drugbank_ids, ligands_names, ligands_ids, ligands_resources
    global ligands_ids_by_names, ligands_resources_by_names
    global targets_ids, targets_resources
    
    ligands_unii = []  # List of ligands' UNII ids
    ligands_drugbank_ids = []  # List of lists of Drugbank ids (one ligand could have several)
    ligands_names = []  # List of usual names
    ligands_ids = []  # List of lists of ids in different DBs
    ligands_resources = []  # List of lists of resources in different DBs
    ligands_ids_by_names = []  # Dictionary name:list of ids
    ligands_resources_by_names = []  # # Dictionary name:list of resources
    
    targets_ids = []  # List of lists of lists of ids in different DBs
    targets_resources = []  # List of lists of lists of resources in different DBs
    unii = ""
    # Flag of FDA approval
    fda_approved = False
    # Database iteration (DB is too big to parse it directly)
    # https://www.drugbank.ca/docs/drugbank.xsd -- scheme of the base
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
                if db_tag(item, "targets"):
                    for it1 in list(item):
                        for it2 in list(it1):
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
                                        t_ids = []
                                        t_resources.append(t_resource)
                                        t_resource = []

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

            targets_ids.append(t_ids)
            targets_resources.append(t_resources)
            
    # Save obtained data        
    dump_info_db(root)


# Directory where all data placed
root = '/media/anton/b8150e49-6ff0-467b-ad66-40347e8bb188/anton/BACHELOR'
ROOT_PATH = Path(root)
DRUGBANK_PATH = Path(root) / 'Drugbank'

process_drugbank(ROOT_PATH)
#update_drugbank(DRUGBANK_PATH)