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
from Bio import SeqIO
#from Bio.SubsMat.MatrixInfo import *
import Bio.SubsMat.MatrixInfo

def get_fasta_similiarity(name1, name2):
    """Calculates similiarity of two sequences in fasta using blosum62 matrix and paths to fastas
    Input - paths to fasta sequences
    Output - similiarity and identity
    """
    #print('Available matrices:', Bio.SubsMat.MatrixInfo.available_matrices)
    #print('Which one would you like to use? Type [Enter] to use blosum62')
    #align_matr = input()
    #if align_matr == "":
    align_matr = "blosum62"
    print(align_matr)
    seq1 = SeqIO.read(name1, "fasta")
    seq2 = SeqIO.read(name2, "fasta")
    alignments = pairwise2.align.globalds(seq1.seq, seq2.seq, blosum62, -10, -0.5)
    alignments_id = pairwise2.align.globalms(seq1.seq, seq2.seq, 1, 0, 0, 0)
    print("Number of alignments = " + str(len(alignments)))
    print(pairwise2.format_alignment(*alignments[0]))
    sim = alignments[0][2]
    ident = str(alignments_id[0][2]).split('.')[0]
    return sim, ident


def get_TMscore_and_RMSD(pdb1_path, pdb2_path):
    """Returns (TM-score, RMSD) of two pdbs, whose paths are in INPUT
    0.0 < TM-score < 0.17, random structural similarity                 *
    0.5 < TM-score < 1.00, in about the same fold 
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
    

def get_SMILES_sim():
    """!!! fix the function"""
    # Read and Concaonate the csv's
    df = pd.read_csv('first.csv')

    # Proof and make a list of SMILES
    df_smiles = df['smiles']
    c_smiles = []
    for ds in df_smiles:
        try:
            cs = Chem.CanonSmiles(ds)
            c_smiles.append(cs)
        except:
            print('Invalid SMILES:', ds)
    print()

    # Make a list of mols
    ms = [Chem.MolFromSmiles(x) for x in c_smiles]

    # Make a list of fingerprints (fp)
    fps = [FingerprintMols.FingerprintMol(x) for x in ms]

    # The list for the dataframe
    qu, ta, sim = [], [], []

    # Compare all fp pairwise without duplicates
    for n in range(len(fps)-1): # -1 so the last fp will not be used
        s = DataStructs.BulkTanimotoSimilarity(fps[n], fps[n+1:]) # +1 compare with the next to the last fp
        print(c_smiles[n], c_smiles[n+1:]) # witch mol is compared with what group
        # Collect the SMILES and values
        for m in range(len(s)):
            qu.append(c_smiles[n])
            ta.append(c_smiles[n+1:][m])
            sim.append(s[m])
    print()

    # Build the dataframe and sort it
    d = {'query':qu, 'target':ta, 'Similarity':sim}
    df_final = pd.DataFrame(data=d)
    df_final = df_final.sort_values('Similarity', ascending=False)
    print(df_final)

    # Save as csv
    df_final.to_csv('result.csv', index=False, sep=',')
    
#print(get_TMscore_and_RMSD('hive/pdb/1A3B.pdb', 'hive/pdb/1RDQ.pdb')) #1MD7 1RDQ
