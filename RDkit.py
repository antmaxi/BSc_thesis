"""From .csv file of SMILES makes .csv file with sorted correspondent pairwise tanimoto similarities 
(alongside checks validity of SMILES)
About RDkit: https://github.com/rdkit/rdkit/blob/master/Docs/Book/Install.md
Also here is a function to  calculate similiarity between two fasta sequences with BioPython
"""
import pandas as pd
import subprocess

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

from Bio import pairwise2
import Bio.SubsMat.MatrixInfo
from Bio.SubsMat.MatrixInfo import *
# Doesn't work
#from Bio import SeqIO


class get_seq_from_fasta:
    """Get name and string from .fasta (biopython's SeqIO sometimes doesn't work:( ))"""
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
        #print(text.split('TM-score')[3])
        tm_score = float(text.split('TM-score')[3].split()[1])
        rmsd = float(text.split('RMSD')[1].split()[4])
        return tm_score, rmsd
    #print(text[ind_tm:ind_tm+40].split()[2])
    

def get_smiles_similiarity(smiles, list_smiles):
    """Get dictionary sims:smiles_from_list for smiles comparing list_smiles"""
    # Proof and make a list of SMILES
    c_smiles = []
    for ds in list_smiles:
        print(ds)
        try:
            cs = Chem.CanonSmiles(ds)
            c_smiles.append(cs)
        except:
            print('Invalid SMILES:', ds)
    try:
        smiles = Chem.CanonSmiles(smiles)
    except:
        print('Invalid SMILES:', ds)    

    # Make a list of mols
    ms = [Chem.MolFromSmiles(x) for x in c_smiles]

    # Make a list of fingerprints (fp)
    fps = [FingerprintMols.FingerprintMol(x) for x in ms]
    # Input fingerprint
    fp_in = FingerprintMols.FingerprintMol(Chem.MolFromSmiles(smiles))

    # Compare all fps with fp_in
    sim = (DataStructs.BulkTanimotoSimilarity(fp_in, fps[:])) # +1 compare with the next to the last fp
    #print()

    # Build the dataframe and sort it
    d = {'query':smiles, 'target_smiles':list_smiles, 'similarity':sim}
    df_final = pd.DataFrame(data=d)
    df_final = df_final.sort_values('similarity', ascending=False)
    return df_final#dict(zip(df_final['Similarity'], df_final['target']))
    

def get_closest_smiles_name(smiles, ligands_names_and_smiles, k=1):
    """Get k names and smiles of the closest to input smiles, k=1 by default"""
    res = get_SMILES_similiarity(smiles, list(ligands_names_and_smiles.values()))
    result = res[0:k]
    names = []
    for sm in result['target_smiles']:
        for name in ligands_names_and_smiles.keys():
            if ligands_names_and_smiles[name] == sm:
                names.append(name)
    result['Name'] = names
    return result

# Tests
#print(get_TMscore_and_RMSD('hive/pdb/1A3B.pdb', 'hive/pdb/1RDQ.pdb')) #1MD7 1RDQ
#get_SMILES_sim()
#print(get_sequences_similiarity('RDkit/O14733.fasta', 'RDkit/P00533.fasta'))
#print(get_sequences_similiarity('ABB', 'ABC'))
#a = get_SMILES_similiarity('ClCCNC(=O)N(CCCl)N=O', 
#                       ('CN1CCC[C@@H]1CCO[C@](C)(C1=CC=CC=C1)C1=CC=C(Cl)C=C1', 'ClCCNC(=O)N(CCCl)N=O'))
#print(a)
