
""" Here are the functions: 
    to calculate similarity between: 
        two fasta sequences (using BioPython)
        one fasta and the whole database targets (using BioPython)
        SMILES and list of SMILES (using RDkit)
        SMILES and the whole database ligands (using RDkit)
        ligand and the whole database (using pybel fingerprints)
        two protein structure files (TM-score and RMSD, using TM-align) # in process

RDkit: installation https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md
       how to use https://www.rdkit.org/docs/Cookbook.html
TM-align: https://zhanglab.ccmb.med.umich.edu/TM-align/
Open Babel: http://openbabel.org/docs/current/UseTheLibrary/Python_Pybel.html
Biopython: https://biopython.org/
           about alignments http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
"""

import os
import pandas as pd
import subprocess
import json
from pathlib import Path
import ntpath


from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

import openbabel
import pybel

from Bio import pairwise2
import Bio.SubsMat.MatrixInfo
from Bio.SubsMat.MatrixInfo import *

import DATABASES_SMILES as db
#import Drugbank as dr

from Bio import SeqIO

# The same function as in DATABASES_SMILES.py
def load_info_db_from_namelist(namelist, root):
    """Load listed in namelist names.txt collected from Drugbank data as json files from root/Drugbank_extracted"""
    # All names of files to be loaded from root/Drugbank_extracted with name.txt, where name is from names
    name_full = str(Path(root) / 'Drugbank_extracted')
    for name in namelist:
        with open(str(Path(name_full) / (name + ".txt")), 'r') as f:
            exec('global ' + name + '\n' + name + ' = json.load(f)')
            
##################   SEQUENCE SIMILARITY (for targets)   ####################################

# Useful if Bio.SeqIO doesn't work somehow, otherwise useless    
class get_seq_from_fasta_file:
    """ Get name (first line of .fasta) and string from .fasta (Biopython's SeqIO sometimes doesn't work:( ))"""
    def __init__(self, path):
        with open(path, 'r') as f:
            seq = ''
            for line in f.readlines():
                if line[0] == '>':
                    name = line.rstrip()
                else:
                    seq += line.rstrip()
            self.name = name
            self.seq = seq


def get_sequences_similarity(input1, input2, align_matrix='blosum62', verbose=False):
    """ Calculates similarity of two sequences using align_matrix from Biopython (blosum62 by default) and 
    paths to fastas or sequences
    Input -  sequences or paths to single fastas
    Output - integers similiarity and identity
    """
    # Draft for using different substitution matrices
    #print('Available matrices:', Bio.SubsMat.MatrixInfo.available_matrices)
    #print('Which one would you like to use? Type [Enter] to use blosum62')
    #align_matr = input()
    # Process input data
    # If first input is path to fasta
    if input1.split('.')[-1] == 'fasta':
        # If import from Bio import SeqIO works
        seq1 = Bio.SeqIO.read(input1, "fasta")
        # If Bio.SeqIO doesn't work
        #seq1 = get_seq_from_fasta_file(input1)
        str1 = seq1.seq        
    # If first input is string
    else:
        str1 = input1
     # If second input is path to fasta    
    if input2.split('.')[-1] == 'fasta':
        # If import from Bio import SeqIO works
        seq2 = Bio.SeqIO.read(input2, "fasta")
        # If Bio.SeqIO doesn't work
        #seq2 = get_seq_from_fasta_file(input2)     
        str2 = seq2.seq
    # If second input is string
    else:
        str2 = input2
    # Get needed matrix
    exec('matr_bio = Bio.SubsMat.MatrixInfo.' + align_matrix, globals())
    # Make an alignment
    try:
        alignments = pairwise2.align.globalds(str1, str2, matr_bio, -10, -0.5)  
        alignments_id = pairwise2.align.globalms(str1, str2, 1, 0, 0, 0) 
        # Print info
        if verbose:
            print("Matrix " + align_matrix + ", number of alignments = " + str(len(alignments)))
            print(pairwise2.format_alignment(*alignments[0]))
        sim = int(alignments[0][2])
        ident = int(str(alignments_id[0][2]).split('.')[0])
        return sim, ident
    except:
        print('Wasn\'t able to compare with ', input2)
        return -1000, 0


def get_closest_fastas_in_fasta_file_from_fasta_or_seq(fasta, path_to_data_in_fasta, 
                                                       k=0, align_matrix='blosum62', sim_or_ident=True):
    """ Returns k or all (if k == 0) of closest to input fasta molecules from path_to_data_in_fasta multi-fasta.
    OUTPUT -- dataframe: 'query':repeated input fasta, 
                    'position_in_fasta': position in input file (to find later needed info)
                     'similarity':similarity, 'identity': identity, 
                     'sequence': sequence of compared target, 'name':name of compared target
            Also writes the best alignment
    INPUT -- fasta -- input a/a sequence of path to single fasta file,
            path_to_data_in_fasta -- fasta file to compare with, 
            k -- number of the best to find (k == 0 if want to get all), 
            sim_or_ident == True => get closest by similiarity, othrwise by identity
    """
    # Load fastas to compare with
    records = list(SeqIO.parse(path_to_data_in_fasta, "fasta"))
    # Process when input is path to fasta file
    if fasta.split('.')[-1] == '.fasta':
        seq = get_seq_from_fasta_file(fasta).seq
    # If input is sequence string
    else:
        seq = fasta
    # Get similarities and identities for all targets in Drugbank
    similarity_list = [] 
    identity_list = []
    seq_list = []
    name_list = []
    for ind, element in enumerate(records):
        if element.seq == seq and k == 1:
            sim, ident = get_sequences_similarity(seq, element.seq)
            d =  {'query':fasta, 'position_in_fasta':ind,  
                  'similarity':sim, 'identity':ident,
                  'sequence':seq, 'name':element.name,
                 }
            return pd.DataFrame(data=d)
        sim, ident = get_sequences_similarity(seq, element.seq, align_matrix)
        similarity_list.append(sim)
        identity_list.append(ident)
        seq_list.append(element.seq)
        name_list.append(element.name)
    # Create correspondent dictionary and then dataframe
    d = {'query':[fasta]*len(similarity_list), 'position_in_fasta':range(len(similarity_list)), 
         'similarity':similarity_list, 'identity':identity_list,
         'sequence':seq_list, 'name':name_list,
        }
    res = pd.DataFrame(data=d)
    if sim_or_ident:
        res = res.sort_values('similarity', ascending=False)
    else:
        res = res.sort_values('identity', ascending=False)
    if k < res.size and k != 0:
        result = res[0:k]
    else:
        result = res
    get_sequences_similarity(seq, get_element_of_fasta_by_number(path_to_data_in_fasta, 
                                                                 result.iloc[0].loc['position_in_fasta']),
                                                                verbose=True)
    return result
                   

def get_closest_fastas_from_uniprot(uniprot, path_to_data_in_fasta, k=0, sim_or_ident=True):
    """ Returns k of closest to fasta of input molecule (with uniprot Uniprot ID) molecules 
    from path_to_data_in_fasta multi-fasta.
    OUTPUT -- 'query':repeated input fasta, 
                    'position_in_fasta': position in input file (to find later needed info)
                     'similarity':similarity, 'identity': identity, 
                     'sequence': sequence of compared target, 'name':name of compared target
    INPUT -- path to fasta or sequence, k -- number of the best to find (k == 0 if want to get all), 
            sim_or_ident == True => get closest by similiarity, othrwise by identity
    """
    seq = db.get_seq_from_uniprot(uniprot)
    return get_closest_fastas_in_fasta_file_from_fasta_or_seq(seq, path_to_data_in_fasta, k, sim_or_ident)

                   
def get_element_of_fasta_by_number(path_to_data_in_fasta, n):
    """ Returns SeqIO record of n-th order from multi-fasta file.
    OUTPUT -- SeqIO fasta sequence element
    INPUT -- path to fasta file with compared fastas, number of needed element in this fasta."""
    # Load fastas with thich compared           
    records = list(SeqIO.parse(path_to_data_in_fasta, "fasta"))
    return records[n]                           


def print_closest_fastas(path_to_data_in_fasta, list_of_numbers):
    """ Get name and sequence of elements with numbers in 'list_of_numbers' from 'path_to_data_in_fasta'"""
    records = list(SeqIO.parse(path_to_data_in_fasta, "fasta"))
    for n in list_of_numbers:
        elem = records[n]
        print('Name = ', elem.name)
        print('Seq = ', elem.seq)


##################   SMILES SIMILARITY (for ligands)   ####################################

def get_smiles_similiarity(smiles, list_smiles):
    """ Get dataframe 'query':input SMILES, 
                    'target_smiles':list_smiles_cleaned, 
                    'similarity':similarities correspondently
    INPUT -- SMILES (smiles) and list of smiles to compare with (list_smiles)
    """
    # Proof and make a list of SMILES
    c_smiles = []
    # Delete Nones
    list_smiles_cleaned = [i for i in list_smiles if i]
    # List of indices to delete because SMILES are invalid
    del_indices = []
    for ind, ds in enumerate(list_smiles_cleaned):
        try:
            cs = Chem.CanonSmiles(ds)
            c_smiles.append(cs)
        except:
            # Delete smiles if it's invalid
            del_indices.append(ind)
            print('Invalid SMILES, deleted from list:', ds)
    # Delete elements starting from end
    for ind in del_indices[::-1]:
        del list_smiles_cleaned[ind]
    try:
        smiles = Chem.CanonSmiles(smiles)
    except:
        print('Invalid Input SMILES:', ds)
        return -1

    # Make a list of mols
    ms = [Chem.MolFromSmiles(x) for x in c_smiles]

    # Make a list of fingerprints (fp)
    fps = [FingerprintMols.FingerprintMol(x) for x in ms]
    # Input fingerprint
    fp_in = FingerprintMols.FingerprintMol(Chem.MolFromSmiles(smiles))

    # Compare all fps with fp_in
    sim = (DataStructs.BulkTanimotoSimilarity(fp_in, fps[:]))
    #print()

    # Build the dataframe and sort it
    print(len([smiles]*len(sim)), len(list_smiles_cleaned), len(sim))
    d = {'query':[smiles]*len(sim), 'smiles':list_smiles_cleaned, 'similarity':sim}
    df_final = pd.DataFrame(data=d)
    df_final = df_final.sort_values('similarity', ascending=False)
    return df_final#dict(zip(df_final['Similarity'], df_final['target']))


def get_closest_smiles_names(smiles, k=1):
    """ Get k names and smiles of the closest to input smiles, k=1 by default
    INPUT -- SMILES (smiles), k -- number of the best smiles to find (k == 0 if want to get all)
    OUTPUT -- dataframe of similar by smiles ligands: 
            'name' -- names of ligands
            'smiles' -- SMILES of ligands
            'query' -- input SMILES (same for all)
            'similarity' -- level of similarity (1 - identical, 0 - abs. different)
    """
    # Load needed dictionary
    global root
    db.load_info_db_from_namelist(['ligands_names_and_smiles'], root)
    # Delete ligands with None smiles
    dict_cleaned = {k: v for k, v in ligands_names_and_smiles.items() if v is not None}
    # Get dataframe of sorted by descending similarity smiles
    try:
        res = get_smiles_similiarity(smiles, list(dict_cleaned.values()))
    except:
        if res == -1:
            print('Input SMILES is invalid, abort')
            return -1
    # Take only needed amount of smiles
    if k < res.size and k != 0:
        result = res[0:k]
    else:
        result = res
    # Get names of correspondent ligands
    names = []
    for sm in result['smiles']:
        for name in ligands_names_and_smiles.keys():
            if ligands_names_and_smiles[name] == sm:
                names.append(name)
    result['name'] = names
    return result


##################   Structure-Fingerprints SIMILARITY (for ligands)   ####################################

def change_extension(path, new_ext):
    """ Returns path with changed extension
    INPUT - path and new extension
    OUTPUT -- path with changed to new_ext extension, if path had extension. If not returns -1"""
    f, f_ext = os.path.splitext(path)
    if f_ext == '':
        print('Path without extension')
        return -1
    return f + new_ext


def extract_approved_sdf(path_to_sdf_from_drugbank, root, overwrite=False):
    """ Extract approved ligands taken from 'ligands_drugbank_ids' to new multi-sdf file
    with changed name as added by _approved. Overwrite - flag of overwriting this file
    """
    # Set path of file with approved sdfs
    path_to_approved_sdf = path_to_sdf_from_drugbank.split('.sdf')[0] + '_approved.sdf'
    # If file with approved structured doesn't exist or should be overwrited
    if not Path(path_to_approved_sdf).is_file() or overwrite:
        # Load list of ids of approved ligands
        load_info_db_from_namelist(['ligands_drugbank_ids'], root)
        sdf_approved = pybel.Outputfile("sdf", path_to_approved_sdf, overwrite=True)
        for mol in pybel.readfile('sdf', path_to_sdf_from_drugbank):
            mol_id = mol.data['DATABASE_ID']
            # Check if ligand is in approved list
            f_approved = False
            for lig in ligands_drugbank_ids:
                if mol_id in lig:
                    f_approved = True
                    break
            if f_approved:
                sdf_approved.write(mol)
        sdf_approved.close()
    return path_to_approved_sdf


def get_closest_ligands_from_3d_structure(path_to_structure, path_to_sdf_approved, root, fptype='fp2', number_to_print=1):
    """ Get sorted descending list of tanimoto coeff from fingerprints and correspondent DB IDs list
    More about fingerprints http://openbabel.org/docs/current/UseTheLibrary/Python_Pybel.html
    Their formats: http://openbabel.org/docs/current/Fingerprints/fingerprints.html#fingerprint-format-details
    INPUT:
        path_to_structure -- path to single .sdf, .pdb or .mol2 structure of molecule
        path_to_sdf_approved -- path to multi-sdf file to compare with
        fptype - type of fingerptint, such as 'ftp2', 'maccs', 'ecfp0' etc.
        list of all available can be taken by 'pybel.fps'
        number_to_print - how many to print
    OUTPUT:
        dataframe: 'Name' : name of compared ligand,'Tanimoto coeff' : corresp similarity,
        'Drugbank ID' : Drugbank ID of compared ligand, 'Fingerprint type' : name of used fingerprint      
    """
    # Fingerprints of ligands
    fps = []
    # Drugbank IDs of ligands
    ids = []
    # Tanimoto coefficients between fingerprints
    tanim = []
    for mymol in pybel.readfile('sdf', path_to_structure):
        fp_mol = mymol.calcfp(fptype)
    for mol in pybel.readfile('sdf', path_to_sdf_approved):
        # Get correspondent fingerprint
        fp = mol.calcfp(fptype)
        fps.append(fp)
        tanim.append(fp_mol | fp)
        # Get ID in Drugbank
        ids.append(mol.data['DATABASE_ID'])
    tanim, ids= zip(*sorted(zip(tanim, ids)))
    # Make them descending
    tanim = tanim[::-1]
    ids = ids[::-1]
    # Get names of ligands
    load_info_db_from_namelist(['ligands_names', 'ligands_drugbank_ids'], root)
    ligands_db_ids_by_names = dict(zip(ligands_names, ligands_drugbank_ids))
    names = []
    for lig_id in ids:
        for name in ligands_db_ids_by_names.keys():
            if lig_id in ligands_db_ids_by_names[name]:
                names.append(name)
    # Make dataframe from obtained data
    data_tuples = list(zip(names, tanim, ids, [fptype]*len(ids)))
    df = pd.DataFrame(data_tuples, columns=['Name','Tanimoto coeff', 'Drugbank ID', 'Fingerprint type'])
    # Print the best k
    print(df[0:number_to_print])
    return df

##################   Structure TM-score and RMSD SIMILARITY (for targets)   ####################################
def get_TMscore_and_RMSD(struct1_path, struct2_path):
    """ Get TM-score and RMSD of two protein structure files
    INPUT -- paths of two protein structure files (.sdf, mol2, pdb)
    OUTPUT -- (TM-score, RMSD) of the files
        if 0.0 < TM-score < 0.17, then random structural similarity 
        if 0.5 < TM-score < 1.00, then in about the same fold 
        Returns 0.0 if no common residues were found
    """
    # Convert files to .pdb if needed, saving in the same directory and changing extension
    pdb1_path = change_extension(struct1_path, '.pdb')
    pdb2_path = change_extension(struct2_path, '.pdb')
    db.convert_single_structure(struct1_path, pdb1_path)
    db.convert_single_structure(struct1_path, pdb2_path)
    # ?? Delete new structures after calculation or not??
    res = subprocess.check_output(['TMscore', pdb1_path, pdb2_path])
    text = res.decode('utf-8')
    #print(text)
    if text.find('TM-score') == -1:
        return(0.0)
    else:
        tm_score = float(text.split('TM-score')[3].split()[1])
        rmsd = float(text.split('RMSD')[1].split()[4])
        return tm_score, rmsd
    
    
def get_closest_targets_from_3d_structure(path_to_structure, path_to_sdf_from_drugbank, k=0):
    """"""
    # Load file with all sdfs to compare
    sdf_drugbank = pybel.readfile('sdf', path_to_sdf_from_drugbank)
    # New path to converted to pdb structure
    new_path = change_extension(path_to_structure, '.pdb')
    # Convert input structure to PDB
    db.convert_structure(path_to_structure, new_path)
    x.calcfp()
    tm_list = []
    rmsd_list = []
    # Make directory where to save all pdbs
    dir_name = str(Path(ntpath.dirname(path_to_structure)) / 'pdbs_ligands')
    db.make_dir(dir_name)


# Tests
#print(get_TMscore_and_RMSD('hive/pdb/1A3B.pdb', 'hive/pdb/1RDQ.pdb')) #1MD7 1RDQ
#print(get_sequences_similarity('P08069.fasta', 'P00533.fasta', verbose=True))
#print(get_sequences_similarity('ABB', 'ABC'))
#a = get_SMILES_similiarity('ClCCNC(=O)N(CCCl)N=O', 
#                       ('CN1CCC[C@@H]1CCO[C@](C)(C1=CC=CC=C1)C1=CC=C(Cl)C=C1', 'ClCCNC(=O)N(CCCl)N=O'))
#print(a)

# Test of seq search
#root = '/media/anton/b8150e49-6ff0-467b-ad66-40347e8bb188/anton/BACHELOR'
root = '/home/anton_maximov/BACHELOR'
uniprot = 'P00533'
#print(db.get_seq_from_uniprot(uniprot))
fasta = '/home/anton_maximov/BACHELOR/P08069.fasta'
path_to_data_in_fasta = '/home/anton_maximov/BACHELOR/Drugbank_extracted/Drugbank_targets.fasta'
#df = get_closest_fastas_in_fasta_file_from_fasta_or_seq(fasta, path_to_data_in_fasta, k=3, sim_or_ident=True)
#print(df['position_in_fasta'], df['similarity'])
#print(df)
path_to_structure = str(Path(root) / 'Drugbank_extracted' / 'SDF_ideal.sdf')
path_to_sdf_from_drugbank =  str(Path(root) / 'Drugbank_extracted' / 'structures.sdf') #str(Path(root) / 'Drugbank_extracted' / 'structures_approved_by_db.sdf')#
#path_to_sdf_approved = extract_approved_sdf(path_to_sdf_from_drugbank, root, overwrite=True)

#get_closest_ligands_from_3d_structure(path_to_structure, path_to_sdf_approved, root,
#                                                          fptype='maccs', number_to_print=5)
