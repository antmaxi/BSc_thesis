"""
About RDkit: https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md
Here is a function to  calculate similarity between two fasta sequences and in the whole database with BioPython
.....
"""
import pandas as pd
import subprocess
import json
from pathlib import Path

# More about RDkit https://www.rdkit.org/docs/Cookbook.html
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

import openbabel

from Bio import pairwise2
import Bio.SubsMat.MatrixInfo
from Bio.SubsMat.MatrixInfo import *

import DATABASES_SMILES as db
#import Drugbank as dr

from Bio import SeqIO  # Doesn't work

def load_info_db_from_namelist(namelist):
    """Load listed in namelist names.txt collected from Drugbank data as json files from root/Drugbank_exracted"""
    # All names of files to be loaded from root/Drugbank_extracted with name.txt, where name is from names
    global root
    name_full = str(Path(root) / 'Drugbank_extracted')
    for name in namelist:
        with open(str(Path(name_full) / (name + ".txt")), 'r') as f:
            exec('global ' + name + '\n' + name + ' = json.load(f)')


def convert_structure(in_format, path_in, out_format, path_out, name):
    """Convert structure in in_format (e.g. 'sdf', 'mol2') 
    from path_in directory with name (e.g. A3551) to out_format (e.g. 'pdb')
    path to input structure file is therefore path_in/name.in_format, to out file analogously
    """
    obConversion = openbabel.OBConversion()
    # Set up conversion
    conversion_initialized = obConversion.SetInAndOutFormats(in_format, out_format)
    if not conversion_initialized:
        print('Conversion couldn\'t be done; maybe, file format is/are not appropriate')
        return -1
    mol = openbabel.OBMol()
    # Read file
    readed_normally = obConversion.ReadFile(mol, str(Path(path_in) / (name + '.' + in_format)))
    if not readed_normally:
        print('File wasn\'t read properly; maybe, not existed')
        return -1
    # Add hydrogens if needed
    #mol.AddHydrogens()
    # Write new file
    written_normally = obConversion.WriteFile(mol, str(Path(path_out) / (name + '.' + out_format)))
    if not written_normally:
        print('Writing was done with error')
        return -1
    
# Useful if Bio.SeqIO doesn't work somehow, otherwise useless    
class get_seq_from_fasta_file:
    """Get name (first line of .fasta) and string from .fasta (biopython's SeqIO sometimes doesn't work:( ))"""
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
    """Calculates similarity of two sequences using align_matrix from Biopython (blosum62 by default) and 
    paths to fastas or sequences
    Input -  sequences or paths to single fastas
    Output - integers similiarity and identity
    """
    # Draft for using different sucstitution matrices
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


        
def get_closest_fastas_in_fasta_file_from_fasta_or_seq(fasta, path_to_data_in_fasta, k=0, sim_or_ident=True):
    """ OUTPUT -- dataframe: 'query':repeated input fasta, 
                            'position_in_fasta': position in input file (to find later needed info)
                     'similarity':similarity, 'identity': identity, 'target_name':name of compared target
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
    for ind, element in enumerate(records):
        if element.seq == seq and k == 1:
            sim, ident = get_sequences_similarity(seq, element.seq)
            return pd.DataFrame(data={d: 'query':fasta, 'position_in_fasta':ind, 
         'similarity':sim, 'identity':ident})
        sim, ident = get_sequences_similarity(seq, element.seq)
        similarity_list.append(sim)
        identity_list.append(ident)
        print(sim, ident)
    # Create correspondent dictionary and then dataframe
    d = {'query':[fasta]*len(similarity_list), 'position_in_fasta':range(len(similarity_list)), 
         'similarity':similarity_list, 'identity':identity_list}
    res = pd.DataFrame(data=d)
    if sim_or_ident:
        res = res.sort_values('similarity', ascending=False)
    else:
        res = res.sort_values('identity', ascending=False)
    if k < res.size and k != 0:
        result = res[0:k]
    else:
        result = res
    return result
                   

def get_closest_fastas_from_uniprot(uniprot, path_to_data_in_fasta, k=0, sim_or_ident=True):
    """ OUTPUT -- dataframe: 'query':repeated input fasta, 
                            'position_in_fasta': position in input file (to find later needed info)
                     'similarity':similarity, 'identity': identity, 'target_name':name of compared target
    INPUT -- path to fasta or sequence, k -- number of the best to find (k == 0 if want to get all), 
            sim_or_ident == True => get closest by similiarity, othrwise by identity
    """
    seq = db.get_seq_from_uniprot(uniprot)
    return get_closest_fastas_in_fasta_file_from_fasta_or_seq(seq, path_to_data_in_fasta, k, sim_or_ident)

                   
def get_element_of_fasta_by_number(path_to_data_in_fasta, n):
    """OUTPUT -- SeqIO fasta sequence element
    INPUT -- path to fasta file with compared fastas, number of needed element"""
    # Load fastas with thich compared           
    records = list(SeqIO.parse(path_to_data_in_fasta, "fasta"))
    return records[n]                           


def get_TMscore_and_RMSD(pdb1_path, pdb2_path):
    """Returns (TM-score, RMSD) of two pdbs, whose paths are in INPUT
    if 0.0 < TM-score < 0.17, then random structural similarity 
    if 0.5 < TM-score < 1.00, then in about the same fold 
    Returns 0.0 if no common residues were found
    """
    res = subprocess.check_output(['TMscore', pdb1_path, pdb2_path])
    text = res.decode('utf-8')
    if text.find('TM-score') == -1:
        return(0.0)
    else:
        tm_score = float(text.split('TM-score')[3].split()[1])
        rmsd = float(text.split('RMSD')[1].split()[4])
        return tm_score, rmsd
    
    
def get_closest_target_from_pdb(pdb):
    """"""
    


def get_smiles_similiarity(smiles, list_smiles):
    """Get dataframe 'query':input SMILES, 
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
    """Get k names and smiles of the closest to input smiles, k=1 by default
    INPUT -- SMILES (smiles), k -- number of the best smiles to find (k == 0 if want to get all)
    OUTPUT -- dataframe of similar by smiles ligands: 
            'name' -- names of ligands
            'smiles' -- SMILES of ligands
            'query' -- input SMILES (same for all)
            'similarity' -- level of similarity (1 - identical, 0 - abs. different)
    """
    # Load needed dictionary
    global root
    load_info_db_from_namelist(['ligands_names_and_smiles'])
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

# Tests
#print(get_TMscore_and_RMSD('hive/pdb/1A3B.pdb', 'hive/pdb/1RDQ.pdb')) #1MD7 1RDQ
#get_SMILES_sim()
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
df = get_closest_fastas_in_fasta_file_from_fasta_or_seq(fasta, path_to_data_in_fasta, k=3, sim_or_ident=True)
print(df['position_in_fasta'], df['similarity'])
print(df)
