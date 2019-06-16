"""From .csv file of SMILES makes .csv file with sorted correspondent pairwise tanimoto similarities 
(alongside checks validity of SMILES)
About RDkit: https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md
Also here is a function to  calculate similiarity between two fasta sequences with BioPython
"""
import pandas as pd
import subprocess
import json

# More about RDkit https://www.rdkit.org/docs/Cookbook.html
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

from Bio import pairwise2
import Bio.SubsMat.MatrixInfo
from Bio.SubsMat.MatrixInfo import *

from pathlib import Path

#import DATABASES_SMILES as db
#import Drugbank as dr
# Doesn't work
#from Bio import SeqIO


class get_seq_from_fasta:
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


def get_sequences_similiarity(input1, input2):
    """Calculates similiarity of two sequences using blosum62 matrix and 
    paths to fastas or sequences
    Input - paths to fasta sequences or these sequences
    Output - integers similiarity and identity
    """
    # Draft for using different sucstitution matrices
    #print('Available matrices:', Bio.SubsMat.MatrixInfo.available_matrices)
    #print('Which one would you like to use? Type [Enter] to use blosum62')
    #align_matr = input()
    #if align_matr == "":
    align_matr = "blosum62"
    print(f'Using {align_matr}')
    # If input is paths to fasta
    if input1.split('.')[-1] == 'fasta':
        # If import from Bio import SeqIO works
        #seq1 = Bio.SeqIO.read(name1, "fasta")
        #seq2 = Bio.SeqIO.read(name2, "fasta")
        seq1 = get_seq_from_fasta(input1)
        str1 = seq1.seq
        seq2 = get_seq_from_fasta(input2)
        str2 = seq2.seq
    # If input is strings
    else:
        str1 = input1
        str2 = input2
    alignments = pairwise2.align.globalds(str1, str2, blosum62, -10, -0.5)
    alignments_id = pairwise2.align.globalms(str1, str2, 1, 0, 0, 0)
    print("Number of alignments = " + str(len(alignments)))
    print(pairwise2.format_alignment(*alignments[0]))
    sim = int(alignments[0][2])
    ident = int(str(alignments_id[0][2]).split('.')[0])
    return sim, ident


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
    

def get_closest_smiles_names(smiles, ligands_names_and_smiles, k=1):
    """Get k names and smiles of the closest to input smiles, k=1 by default
    INPUT -- SMILES (smiles) and dictionary names:SMILES (ligands_names_and_smiles)
            k -- number of smiles to find
    OUTPUT -- dataframe of similar by smiles ligands: 
            'name' -- names of ligands
            'smiles' -- SMILES of ligands
            'query' -- input SMILES (same for all)
            'similarity' -- level of similarity (1 - identical, 0 - abs. different)
    """
    # Delete ligands with None smiles
    dict_cleaned = {k: v for k, v in ligands_names_and_smiles.items() if v is not None}
    # Get dataframe of sorted by descending similarity smiles
    try:
        res = get_smiles_similiarity(smiles, list(dict_cleaned.values()))
    except:
        if res == -1:
            print('Input SMILES is invalid, abort')
    # Take only needed amount of smiles
    if k < res.size:
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
#print(get_sequences_similiarity('RDkit/O14733.fasta', 'RDkit/P00533.fasta'))
#print(get_sequences_similiarity('ABB', 'ABC'))
#a = get_SMILES_similiarity('ClCCNC(=O)N(CCCl)N=O', 
#                       ('CN1CCC[C@@H]1CCO[C@](C)(C1=CC=CC=C1)C1=CC=C(Cl)C=C1', 'ClCCNC(=O)N(CCCl)N=O'))
#print(a)
# Test of SMILES search
name = 'ligands_names_and_smiles'
root = '/media/anton/b8150e49-6ff0-467b-ad66-40347e8bb188/anton/BACHELOR'
name_full = str(Path(root) / 'Drugbank_extracted')
with open(str(Path(name_full) / (name + ".txt")), 'r') as f:
    exec('global ' + name + '\n' + name + ' = json.load(f)')
#df = get_closest_smiles_names('ClC1=CC=CC=C1CN1CCCC2=C(C1)C=CS2', ligands_names_and_smiles, 3)
#print(df, df['name'], df['smiles'],df['query'], df['similarity'])
