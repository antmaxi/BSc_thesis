# %load IsoMIF.py

"""Using 2 pdb ids and numbers of two clefts from the biggest gets tanimoto coef between them"""

import os
import subprocess
import requests
from pathlib import Path  # Library for easier paths processing
import mmap

global GSL_PATH, GET_CLEFT_PATH, ISOMIF_PATH, REDUCE_PATH, SYSTEM_NAME
global HIVE_PATH, PDB_PATH, MIF_NAME, ISOMIF_NAME, N_cavities

PDB_PATH = None
N_cavities = 1


def make_dir_from_list(dirList):
    """Make directories with paths from dirList"""
    for dirName in dirList:
        if not os.path.exists(dirName):
            os.mkdir(dirName)
            print("Directory ", dirName,  " Created ")
        else:
            pass
            #print("Directory " , dirName ,  " already exists")


def make_hive(root_isomif):
    """Create all needed for IsoMIF data storage directories"""
    List = []
    List.append(str(Path(root_isomif)))
    List.append(str(Path(root_isomif) / 'hive'))
    for last_name in ('clefts', 'match', 'matchView',
                      'mifs', 'mifView', 'pdb'):
        List.append(str(Path(root_isomif) / 'hive' / last_name))
    make_dir_from_list(List)


def download_url(url, path=None, name=None):
    """Download from url to 'path/name', making path directory, if not existed"""
    r = requests.get(url, allow_redirects=True)
    if path:
        paths = []
        paths.append(path)
        make_dir_from_list(paths)
        open(os.path.join(paths[0], name), 'wb').write(r.content)
    return r.content.decode('utf-8')

def download_pdb(pdb, pdb_path=None):
    """Download pdb (e.g. 1A46) and save to PATH directory (with name 1A46.pdb)
    pdb -- identifier in PDB
    pdb_path -- directory-string where to save
    """
    if pdb_path is None:
        pdb_path = PDB_PATH
        
    path_uniprot = Path(pdb_path)
    name = pdb + ".pdb"
    url = "https://files.rcsb.org/download/" + pdb + ".pdb"
    full_name = path_uniprot / name
    # Check if .pdb is already downloaded, if not => download
    if not full_name.is_file():
        download_url(url, str(path_uniprot), name)


def isomif_init(root, gsl_path=None, reduce_path=None,
                system_name='linux_x86_64', n_cavities=1, 
                overwrite_init=False):
    """ Prepare everything (paths, default settings)  for the IsoMIF work. 
    If only root is in input, then try to find root/init.txt and get settings from it
    
    Input -- strings:
            root - Root for making storage of all needed for Isomif data
            gsl_path - Path to GSL library with 'lib' and 'include' subdirectories
            reduce_path - Path to reduce program
            system_name - string of system name, could be 'linux_x86_64' (default) or 'mac_x86_64'
            
        n_cavities - Default number of cavities to find
        boolean overwrite_init -- True if needed to overwrite old init.txt
    """
    global GSL_PATH, GET_CLEFT_PATH, ISOMIF_PATH, REDUCE_PATH, SYSTEM_NAME
    global HIVE_PATH, PDB_PATH, MIF_NAME, ISOMIF_NAME, N_cavities
    
    path_init = str(Path(root)  / 'Isomif' / 'init.txt')
    # If need to initialize from input and create/overwrite init file 
    if not Path(path_init).is_file() or overwrite_init:
        with open(str(Path(root) / 'Isomif' / 'init.txt'), 'w+') as init_file:
            init_file.write('\n'.join([root, gsl_path, reduce_path, 
                                  system_name, str(n_cavities)]))
    # If init.txt exists and we don't need to overwrite it
    else:
        with open(Path(root) / 'Isomif' / 'init.txt', 'r') as init_file:
            lines = init_file.readlines()
            root = lines[0].rstrip()
            gsl_path = lines[1].rstrip()
            reduce_path = lines[2].rstrip()
            system_name = lines[3].rstrip()
            n_cavities = int(lines[4])

    # Further global paths needed to be defined
    #
    # Root for making storage of all needed for Isomif data    
    root_isomif = str(Path(root) / 'Isomif')
    # Path to GSL library with 'lib' and 'include' subdirectories
    GSL_PATH = gsl_path
    # Path to Get_Cleft-master directory
    GET_CLEFT_PATH = str(Path(root_isomif) / 'Get_Cleft-master')
    # Path to Isomif-master directory
    ISOMIF_PATH = str(Path(root_isomif) / 'IsoMif-master')
    # Path to reduce program
    REDUCE_PATH = reduce_path
    # System name, could be 'linux_x86_64' or 'mac_x86_64'
    SYSTEM_NAME = system_name
    # Default number of cavities to find
    N_cavities = n_cavities

    # Where all data stored
    HIVE_PATH = str(Path(root_isomif) / 'hive')
    PDB_PATH = str(Path(root_isomif) / 'hive' / 'pdb')

    MIF_NAME = 'mif_' + SYSTEM_NAME + '_compiled'
    ISOMIF_NAME = 'isomif_' + SYSTEM_NAME + '_compiled'

    # Creating needed for IsoMIF directories inside ROOT
    make_hive(root_isomif)
        

def isomif_compile():
    """Compile Get_Cleft for this OS. Initialisation by isomif_init function is needed before this"""
    subprocess.check_output(['gcc', str(Path(GET_CLEFT_PATH) / 'Get_Cleft.c'),
                             '-o', str(Path(GET_CLEFT_PATH) / 'Get_Cleft'),
                             '-O3', '-lm',
                             ])
    print('Get_Cleft compiled')
    # Compile MIF and IsoMIF for linux_x86_64, needs path to gsl
    subprocess.check_output(['g++', str(Path(ISOMIF_PATH) / 'mif.cpp'),
                            '-o', str(Path(ISOMIF_PATH) / MIF_NAME),
                            '-O3', '-lm',
                            ])
    print('MIF compiled')
    subprocess.check_output(['g++', str(Path(ISOMIF_PATH) / 'isomif.cpp'),
                            '-o', str(Path(ISOMIF_PATH) / ISOMIF_NAME),
                            '-O3', '-lm', '-lgsl', '-lgslcblas', '-L', str(Path(GSL_PATH) / 'lib'),
                            '-I', str(Path(GSL_PATH) / 'include'),
                            ])
    print('IsoMIF compiled')

    
def get_cavities(pdb, n_cavities=N_cavities, pdb_path=None):
    """Finding n_cavities of the biggest cavities in pdb"""
    if pdb_path is None:
        pdb_path = PDB_PATH
    subprocess.check_output([str(Path(GET_CLEFT_PATH) / 'Get_Cleft'),
                            '-p', str(Path(pdb_path) / (pdb + '.pdb')),
                            '-o', str(Path(HIVE_PATH) / 'clefts' / pdb),
                            '-s', '-t', str(n_cavities),
                            ])


def add_hydrogens_by_reduce(pdb, pdb_path=None):
    """Add hydrogens to .pdb using reduce, save in the same directory adding 'h' to pdb name"""
    if pdb_path is None:
        pdb_path = PDB_PATH
    # Check if hydrogens were added before
    if not Path(str(Path(pdb_path) / (pdb + 'h.pdb'))).is_file():
        print(str(Path(pdb_path) / (pdb + '.pdb')))
        subprocess.check_output(REDUCE_PATH + \
                                ' -p ' + str(Path(pdb_path) / (pdb + '.pdb')) + \
                                ' > ' + str(Path(pdb_path) / (pdb + 'h.pdb')), shell=True
                                )  # Somehow with list of parameters does't work properly
        print(f'Added hydrogens to {pdb}.pdb')
    else:
        print(f'Were before added hydrogens to {pdb}.pdb')


def calc_mif(pdb,  n_cavities=N_cavities, pdb_path=None, res_mif=1, vis=True):
    """Get MIFs for n_cavities of pdb placed in pdb_path. Vis=True to make visualisation
    res -- Resolution
    0 => 2 Ang, 1 => 1.5 Ang, 2 => 1.0 Ang, 3 => 0.5 Ang"""
    if pdb_path is None:
        pdb_path = PDB_PATH
    f_calculated_something = False
    for i in range(1, n_cavities + 1):
        # Check if MIF is already calculated
        if not Path(str(Path(HIVE_PATH) / 'mifs' / (pdb + 'h_' + str(i) + '.mif'))).is_file():
            f_calculated_something = True
            subprocess.check_output([str(Path(ISOMIF_PATH) / MIF_NAME),
                                    '-p', str(Path(PDB_PATH) / (pdb + 'h.pdb')),
                                    '-g', str(Path(HIVE_PATH) / 'clefts' / (pdb + '_sph_' + str(i) + '.pdb')),
                                    '-o', str(Path(HIVE_PATH) / 'mifs'),
                                    'z', str(res_mif),
                                    ])
            # Rename files, as by default makes just smth like 1E8Xh.mif and 1E8Xh_cpy.pdb despite number of cleft
            os.rename(str(Path(HIVE_PATH) / 'mifs' / (pdb + 'h.mif')), 
                      str(Path(HIVE_PATH) / 'mifs' / (pdb + 'h_' + str(i) + '.mif')))
            os.rename(str(Path(HIVE_PATH) / 'mifs' / (pdb + 'h_cpy.pdb')), 
                      str(Path(HIVE_PATH) / 'mifs' / (pdb + 'h_' + str(i) + '_cpy.pdb')))
            # Produce file for visualisation
            if vis:
                subprocess.check_output(['perl', str(Path(ISOMIF_PATH) / 'mifView.pl'),
                                         ' -m ', str(Path(HIVE_PATH) / 'mifs' / (pdb + 'h_' + str(i) + '.mif')),
                                         ' -o ', str(Path(HIVE_PATH) / 'mifView'), #/ (pdb + str(i))
                                    ])
                print('perl' + str(Path(ISOMIF_PATH) / 'mifView.pl')+ \
                                         ' -m ' + str(Path(HIVE_PATH) / 'mifs' / (pdb + 'h_' + str(i) + '.mif')) + \
                                         ' -o ' + str(Path(HIVE_PATH) / 'mifView'))
        else:
            print(f'Was before calculated MIF for {i}th(st, nd)')
    if f_calculated_something:
        print(f'Calculated MIFs of {n_cavities} the biggest cavities in {pdb}h.pdb')
            

def calc_isomif(pdb1, pdb2, i1=1, i2=1, res_isomif=1, node_variab=2.0, vis=False, res_nodes=1):
    """Calculate ISOMIF of files pdb1 + 'h_' + str(i1) + '.mif' and pdb2 + 'h_' + str(i2) + '.mif'
    in HIVE_PATH / 'mifs'
    """
    # SOMETIMES COLLAPSES WHEN NOT ENOUGH MEMORY
    pdbh1 = pdb1 + 'h'
    pdbh2 = pdb2 + 'h'
    # Check if IsoMIF isn't calculated already
    if not Path(str(Path(HIVE_PATH) / 'match' / (pdbh1 + '_' + str(i1) + '_match_' + \
                                                 pdbh2 + '_' + str(i2) + '.isomif'))).is_file():
        print(' '.join([str(Path(ISOMIF_PATH) / ISOMIF_NAME),
                                '-p1', str(Path(HIVE_PATH) / 'mifs' / (pdbh1 + '_' + str(i1) + '.mif')),
                                '-p2', str(Path(HIVE_PATH) / 'mifs' / (pdbh2 + '_' + str(i2) + '.mif')),
                                '-o', os.path.join(str(Path(HIVE_PATH) / 'match'), ''),
                                '-c', str(res_isomif), '-d', str(node_variab),
                                ]))
        subprocess.check_output([str(Path(ISOMIF_PATH) / ISOMIF_NAME),
                                '-p1', str(Path(HIVE_PATH) / 'mifs' / (pdbh1 + '_' + str(i1) + '.mif')),
                                '-p2', str(Path(HIVE_PATH) / 'mifs' / (pdbh2 + '_' + str(i2) + '.mif')),
                                '-o', os.path.join(str(Path(HIVE_PATH) / 'match'), ''),
                                '-c', str(res_isomif), '-d', str(node_variab),
                                ])
        
        if vis:
            subprocess.check_output(['perl', str(Path(ISOMIF_PATH) / 'isoMifView.pl'),
                                     '-m', str(Path(HIVE_PATH) / 'match' / (pdbh1 + '_' + str(i1) + \
                                                                '_match_' + pdbh2 + '_' + str(i2) + '.isomif')),
                                     '-o', os.path.join(str(Path(HIVE_PATH) / 'matchView'), ''),
                                     '-g', str(res_nodes)
                                    ])
        print(f'Calculated IsoMIF for {i1}th(st, nd) and {i2}th(st, nd) biggest cavities in {pdb1}h.pdb and {pdb2}h.pdb')
    else:
        print(f'Was before calculated IsoMIF for {i1}th(st, nd) and {i2}th(st, nd) biggest cavities in {pdb1}h.pdb and {pdb2}h.pdb') 
        
        
def get_tanimoto_from_isomif_file(pdb1, pdb2, i1, i2):
    """Extract tanimoto similiarity from .txt file - result of IsoMIF"""
    filename = pdb1 + 'h_' + str(i1) + '_match' + pdb2 + 'h_' + str(i2) + '.isomif'
    with open(str(HIVE / 'match' / filename), 'rb', 0) as file, \
         mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ) as s:
        if s.find(b'TANIM') != -1:
            pos = s.find(b'TANIM')
            tanim = float(s[pos+6 : pos+12].decode('utf-8'))
        else:
            return -1
    return tanim


def mif_from_pdb_simple(pdb, n_cavities=N_cavities):
    """Produce MIF files from pdb id"""
    download_pdb(pdb, PDB_PATH)
    add_hydrogens_by_reduce(pdb)
    get_cavities(pdb, n_cavities)
    calc_mif(pdb, n_cavities, vis=True)


def calc_mifs_of_whole_drugbank():
    """"""
# Initialisation of paths and default parameters
# (maybe, it's worth to put initialisation parameters into .txt?)
root = '/home/anton_maximov/BACHELOR'
root_isomif = str(Path(root) / 'Isomif')
#root = '/media/anton/b8150e49-6ff0-467b-ad66-40347e8bb188/anton/BACHELOR/Isomif'
isomif_init(
    # Root directory
    root,
    # GSL directory
    '/home/anton_maximov/gsl',
    # Reduce executive
    os.path.join(root_isomif, 'reduce.3.23.130521'),
    # System name, could be 'linux_x86_64' or 'mac_x86_64'
    'linux_x86_64', 
    # Default number of cavities to find
    1,
    )

# Compilation if hasn't been done before
#isomif_compile()

# Example of getting pdbs from uniprot: 
#uniprot1 = 'P00734'  # Prothrombin
#uniprot2 = 'P00736'  # Cetuximab
#pdb1 = get_pdbs_from_uniprot(uniprot1, PDB_PATH)[1]
#pdb2 = get_pdbs_from_uniprot(uniprot2, PDB_PATH)[2]
#mif_from_pdb_simple(pdb1, i1)
#mif_from_pdb_simple(pdb2, i2)
#calc_isomif(pdb1, pdb2, i1, i2) DOESN'T WORK
#return get_tanimoto_from_isomif_file(pdb1, pdb2, i1, i2)
