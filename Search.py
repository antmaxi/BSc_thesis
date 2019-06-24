
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
import subprocess  # To execute like as from cmd
import json
from pathlib import Path  # To process paths
import ntpath
import pickle


#from rdkit import Chem
#from rdkit import DataStructs
#from rdkit.Chem.Fingerprints import FingerprintMols

import openbabel
import pybel

from Bio import pairwise2  # To make sequence alignments
import Bio.SubsMat.MatrixInfo  # To get info about available distance matrices
from Bio.SubsMat.MatrixInfo import *
from Bio import PDB  # To parse PDB files

import Auxiliary as aux  # Needed for work auxiliary functions
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
    """ Get name (first line of .fasta) and string from .fasta 
    (Biopython's SeqIO sometimes doesn't work due to improper installation:( ))"""
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

            
def get_seq_from_fasta_uniprot_or_seq(input1):
    """ Return sequence from input as sequnce, path to fasta or Uniprot ID."""
    # If input is fasta
    if input1.split('.')[-1] == 'fasta':
        # If import from Bio import SeqIO works
        seq1 = Bio.SeqIO.read(input1, "fasta")
        # If Bio.SeqIO doesn't work
        #seq1 = get_seq_from_fasta_file(input1)
        str1 = seq1.seq        
    # If input is seq or uniprot
    else:
        # Uniprot ID have length of 6, so checking if it is possibly ID or not
        if len(input1) <= 10:
            # Flag of being Uniprot ID
            f_uniprot = False
            for i in range(len(input1)):
                if input1[i].isdigit():
                    f_uniprot = True
            if f_uniprot:
                str1 = aux.get_seq_from_uniprot(input1)
            else:
                print("Is it really such a short protein sequence or invalid Uniprot ID")
                str1 = input1
        else:
            str1 = input1
    return str1


def get_sequences_similarity(input1, input2, align_matrix='blosum62', verbose=False):
    """ Calculates similarity of two inputs (could be raw seq, path to fasta or uniprot ID) 
    using align_matrix from Biopython (blosum62 by default, list of all by Bio.SubsMat.MatrixInfo.available_matrices)
    Input -  sequences, paths to single fastas or Uniprot IDs of proteins
    Output - float similiarity and integer identity
    """
    # Draft for using different substitution matrices
    #print('Available matrices:', Bio.SubsMat.MatrixInfo.available_matrices)
    #print('Which one would you like to use? Type [Enter] to use blosum62')
    #align_matr = input()
    # Process input data
    seq1 = get_seq_from_fasta_uniprot_or_seq(input1)
    seq2 = get_seq_from_fasta_uniprot_or_seq(input2)
    # Get needed matrix
    exec('matr_bio = Bio.SubsMat.MatrixInfo.' + align_matrix, globals())
    # Make an alignment
    try:
        alignments = pairwise2.align.globalds(seq1, seq2, matr_bio, -10, -0.5)  
        alignments_id = pairwise2.align.globalms(seq1, seq2, 1, 0, 0, 0) 
        # Print info
        sim = float(alignments[0][2])
        ident = int(str(alignments_id[0][2]).split('.')[0])
        print(f'Similarity={sim}, identity={ident}')
        if verbose:
            print("Matrix " + align_matrix + ", number of alignments = " + str(len(alignments)))
            print(pairwise2.format_alignment(*alignments[0]))
        return sim, ident
    except:
        print('Smth went wrong with comparison to ', input2)
        return -1000, 0


def get_closest_fastas_in_fasta_file_from_fasta_uniprot_or_seq(input1, path_to_data_in_fasta, 
                                                       k=0, align_matrix='blosum62', sort_by='s'):
    """ Returns k or all (if k == 0) of closest to input fasta molecules from path_to_data_in_fasta multi-fasta.
    OUTPUT -- dataframe: 'query':repeated input fasta, 
                    'position_in_fasta': position in input file (to find later needed info)
                     'similarity':similarity, 'identity': identity, 
                     'sequence': sequence of compared target, 'name':name of compared target
            Also writes the best alignment
    INPUT -- input1 -- input a/a sequence, path to single fasta file or Uniprot ID of protein,
            path_to_data_in_fasta -- fasta file to compare with, 
            k -- number of the best to find (k == 0 if want to get all), 
            sort_by == 's' => sort descending by seimilarity. == 'i' => by identity
    """
    # Load fastas to compare with
    records = list(SeqIO.parse(path_to_data_in_fasta, "fasta"))
    # Process when input is path to fasta file
    seq = get_sequence_from_fasta_uniprot_or_sequence(input1)
    # Get similarities and identities for all targets in Drugbank
    similarity_list = [] 
    identity_list = []
    seq_list = []
    name_list = []
    for ind, element in enumerate(records):
        if element.seq == seq and k == 1:
            sim, ident = get_sequences_similarity(seq, element.seq, align_matrix)
            d =  {'query':fasta, 'position_in_fasta':ind,  
                  'similarity':sim, 'identity':ident,
                  'sequence':seq, 'name':element.name,
                 }
            return pd.DataFrame(data=d)
        sim, ident = get_sequences_similarity(seq, element.seq, align_matrix)
        similarity_list.append(sim)
        identity_list.append(ident)
        seq_list.append(element.seq)
        name_list.append(element.description)
    # Create correspondent dictionary and then dataframe
    d = {'query':[fasta]*len(similarity_list), 'position_in_fasta':range(len(similarity_list)), 
         'similarity':similarity_list, 'identity':identity_list,
         'sequence':seq_list, 'name':name_list,
        }
    res = pd.DataFrame(data=d)
    if sort_by == 's':
        res = res.sort_values('similarity', ascending=False)
    else:
        res = res.sort_values('identity', ascending=False)
        if sort_by != 'i':
            print('Sorted descending by identity. If you need by similarity, corresp. key should be "s"')
    if k < res.size and k != 0:
        result = res[0:k]
    else:
        result = res
    return result

                   
def get_element_of_fasta_by_number(path_to_data_in_fasta, n):
    """ Returns SeqIO record of n-th order from multi-fasta file.
    OUTPUT -- SeqIO fasta sequence element
    INPUT -- path to fasta file with compared fastas, number of needed element in this fasta."""
    # Load records           
    records = list(SeqIO.parse(path_to_data_in_fasta, "fasta"))
    elem = records[n]
    print('Name = ', elem.description)
    print('Seq = ', elem.seq)
    return elem                          


def print_closest_fastas_data(input1, path_to_data_in_fasta, k=0, align_matrix='blosum62', sim_or_ident=True):
    """ Get k or all (if k == 0) closest proteins in file by fasta from seq/uniprot/path to fasta.
    Print alignments of k or 5 (if  k > 5) with input
    """
    # Get seq of input
    seq1 = get_seq_from_fasta_uniprot_or_seq(input1)
    # Get dataframe with sorted by similarity or identity fastas
    df = get_closest_fastas_in_fasta_file_from_fasta_uniprot_or_seq(seq1, path_to_data_in_fasta, 
                                                       k=0, align_matrix='blosum62', sim_or_ident=True)
    print(df)
    # How many to align, from 1 to 5
    if k:
        k1 = min(k, 5)
    else:
        k1 = 5
    res = df[0:k1]
    # Load records
    records = list(SeqIO.parse(path_to_data_in_fasta, "fasta"))
    for row in res.iterrows():
        # Position of this target in whole fasta file
        n = row['position_in_fasta']
        seq2 = row['sequence']
        print('Name = ', records[n].description)
        # Print alignment, sim and indent coeffs
        get_sequences_similarity(seq1, seq2, align_matrix='blosum62', verbose=True)
    return df


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
            print('Invalid SMILES, deleted from list to compare with:', ds)
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


def get_closest_smiles_names(smiles, root, k=1):
    """ Get k names and smiles of the closest to input smiles, k=1 by default
    INPUT -- SMILES (smiles), 
            k -- number of the best smiles to find (k == 0 if want to get all)
            root - where all protocol is located (ligands_names_and_smiles.txt in root/Drugbank_extracted)
    OUTPUT -- dataframe of similar by smiles ligands: 
            'name' -- names of ligands
            'smiles' -- SMILES of ligands
            'query' -- input SMILES (same for all)
            'similarity' -- level of similarity (1 - identical, 0 - abs. different)
    """
    # Load needed dictionary
    load_info_db_from_namelist(['ligands_names_and_smiles'], root)
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


def get_closest_ligands_from_3d_structure(path_to_structure, path_to_sdf_approved, root, 
                                          fptype='fp2', number_to_print=1):
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


##################   Structure TM-score and RMSD SIMILARITY (for targets and complexes)   ####################################
def get_TMscore_and_RMSD_of_proteins_or_complexes(input1, input2, root,
                                                  compare_type='p', verbose=False):
    """ Get TM-score and RMSD of two protein or complex structure files
    INPUT -- paths of two protein structure files (.sdf, mol2, pdb)
            root - save to root/'pdb', if PDB ID is in input
            compare_type == 'p' => comparing proteins
            compare_type == 'c' => comparing complexes, else also try as complexes, but with warning
            
    OUTPUT -- (TM-score, RMSD) of the files
        if 0.0 < TM-score < 0.17, then random structural similarity 
        if 0.5 < TM-score < 1.00, then in about the same fold 
        Returns 0.0 if no common residues were found
    """
    # Convert files to .pdb if needed, saving in the same directory and changing extension
    # Or download PDB file, if PDB ID is as input
    pdb1_path = aux.get_path_to_pdb_from_pdb_id_or_path_to_structure(input1, root)
    pdb2_path = aux.get_path_to_pdb_from_pdb_id_or_path_to_structure(input2, root)
    # ?? Delete new structures after calculation or not??
    
    # Get result of TM-align work in protein and complex comparison types
    if compare_type == 'p':
        res = subprocess.check_output(['TMscore', pdb1_path, pdb2_path])
    else:
        res = subprocess.check_output(['TMscore', '-c', pdb1_path, pdb2_path])
        # Process invalid key
        if compare_type != 'c':
            print(f"Calculated as for complexes, but key was '{compare_type}' and in ['p', 'c']")
            compare_type = 'c'
    text = res.decode('utf-8')        
    # Find needed results in the output of TM-align
    if text.find('TM-score') == -1:
        print(f'Something went wrong when TM-align compared {pdb1} and {pdb2}')
        return(0.0)
    else:
        if text.find('Warning') != -1:
            print(f'When comparing {input1} and {input2}:')
            print(text.split('*')[0])
        if compare_type == 'p':
            tm_score = float(text.split('TM-score')[3].split()[1])
            rmsd = float(text.split('RMSD')[1].split()[4])
            common_res = int(text.split('common=')[1].split()[0])
            
        if compare_type == 'c':
            tm_score = float(text.split('TM-score')[3].split()[1])
            rmsd = float(text.split('RMSD')[1].split()[4])
            common_res = int(text.split('common=')[1].split()[0])
        if verbose:
            print(f'TM-score = {tm_score}, RMSD = {rmsd}, number of common residues in alignment = {common_res}')
            print(text)
        return tm_score, rmsd, common_res


def download_proteomes(root, overwrite=False):
    """Download reviewed proteomes of human, rat and mouse from Uniprot and save them to root/'Uniprot_proteomes'
    with names organism_proteome.fasta
    """
    # Create directory (if needed) where to save 
    uniprot_dir = str(Path(root) / 'Uniprot_proteomes')
    aux.make_dir(uniprot_dir)
    # Ids of species in Uniprot
    list_ids = ['9606', '10116', '10090']
    # Correspondent names
    list_names = ['human', 'rat', 'mouse']
    organism_dict = dict(zip(list_ids, list_names))
    # Download reviewed proteomes with ids from organism_dict.keys
    for organism_id in list_ids: #human, rat, mouse
        url = 'https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:' + organism_id + '&format=fasta'
        aux.download_url(url, uniprot_dir, organism_dict[organism_id] + '_proteome.fasta', overwrite)

        
def get_fasta_from_pdb(pdb, directory_to_save=None):
    """ Download fasta file for PDB structure from PDB ID and save it to directory, creating it if not existed"""
    url = 'https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=' \
                + pdb + '&compressionType=uncompressed'
    # If needed to save
    if directory_to_save:
        aux.make_dir(directory_to_save)
        r = aux.download_url(url, str(Path(directory_to_save)), (pdb + '.fasta'))
    else:
        r = aux.download_url(url)
    return r

def get_best_pdb_of_target(uniprot, root, verbose=False):
    """ Get pdb of target which has the most biggest similarity to its sequence.
    Input:
        uniprot -- Uniprot ID of target
        root -- root of the protocol
    Output:
        
        """
    # Get list of pdbs where this uniprot is mentioned
    pdbs = aux.get_pdbs_from_uniprot(uniprot)
    # Get sequence of target
    seq = aux.get_seq_from_uniprot(uniprot)
    for pdb in pdbs:
        print()
        print(pdb)
        print()
        # Get fastas attached to pdbs
        pdb_dir = str(Path(root) / 'pdb')
        get_fasta_from_pdb(pdb, pdb_dir)
        for record in SeqIO.parse(str(Path(pdb_dir) / (pdb + '.fasta')), 'fasta'):
            print(record.seq)
            get_sequences_similarity(seq, record.seq, verbose=True)
    
    
def get_pdbs_of_drugbank_targets(path_to_sdf_from_drugbank, root):
    """"""
    
def get_closest_targets_from_3d_structure(path_to_structure, type='DB', k=0):
    """ type -- DB, human, rat, mouse
    """
    # Load file with all sdfs to compare
    sdf_drugbank = pybel.readfile('sdf', path_to_sdf_from_drugbank)
    # New path to converted to pdb structure
    new_path = aux.change_extension(path_to_structure, '.pdb')
    # Convert input structure to PDB
    aux.convert_structure(path_to_structure, new_path)
    x.calcfp()
    tm_list = []
    rmsd_list = []
    # Make directory where to save all pdbs
    dir_name = str(Path(ntpath.dirname(path_to_structure)) / 'pdbs_ligands')
    aux.make_dir(dir_name)

    
##################   Search by TM-score SIMILARITY (for complexes)   ####################################

def get_closest_complexes(input1, sim, sim_min, root, k_print=1):
    """ Get info about closest complexes to input pdb, print about k_print
    Input:
        input1 - path to structure in .pdb, .sdf or .mol2, or PDB ID
        sim - level of similarity of SMILES of ligands to search for their pdbs
        root -- root of the protocol
    Output:
        dictionary with keys (ligand_name, target_uniprot, number) 
                        values (TM-score, RMSD, number of residues in alignment)
    """
    # Load dict of pdbs by name of ligand
    filename, path = produce_name_and_path_of_file_with_sims('all_pdbs_of_all_connections_', sim, sim_min, root)
    # Load dict (name of ligand, uniprot, sim_of_SMILES, sim_min) : [pdbs where is ligand and target]
    with open(path, 'rb') as f:
        connect_dict = pickle.load(f)
    sim_values = []
    connection_keys = []
    # Produce dictionary
    # keys --  (ligand_name, target_uniprot, sim, sim_min, input1) and
    # values -- (TM-score, RMSD, number of residues in alignment)
    # Iterating over pairs (ligand, target)
    for name_uniprot_sim, pdbs in connect_dict.items():     
        # Iteration over pdbs correspondent to one (ligand, target)
        for pdb in pdbs:
            # Get similarity measures for this pdb with input
            tm, rmsd, common_res = get_TMscore_and_RMSD_of_proteins_or_complexes(input1, pdb, root, compare_type='c')
            sim_value = (tm, rmsd, common_res)
            sim_values.append(sim_value)
            connection_key = (name_uniprot_sim[0], name_uniprot_sim[1], 
                              name_uniprot_sim[2], name_uniprot_sim[3],
                             input1)
            connection_keys.append(connection_key)
    res = dict(zip(connection_keys, sim_values))
    # Sort by tm-score
    res_sort = sorted(res.items(), key=lambda e: e[1][1])
    path = str(Path(root) / 'pdb' / '1.txt')
    with open(path, 'wb') as f:
        pickle.dump(res_sorted, f, pickle.HIGHEST_PROTOCOL)
    print(res_sorted[0:k_print])
    return path


# Tests
#root = '/media/anton/b8150e49-6ff0-467b-ad66-40347e8bb188/anton/BACHELOR'
root = '/home/anton_maximov/BACHELOR'

# Test of SMILES search
#df = get_closest_smiles_names('ClCCNC(=O)N(CCCl)N=O', root)
#print(df)

# Test of seq search
uniprot = 'P00533'
uniprot = 'O43451'
#print(aux.get_seq_from_uniprot(uniprot))
fasta = '/home/anton_maximov/BACHELOR/P08069.fasta'
path_to_data_in_fasta = '/home/anton_maximov/BACHELOR/Drugbank_extracted/Drugbank_targets.fasta'
#df = get_closest_fastas_in_fasta_file_from_fasta_uniprot_or_seq(fasta, path_to_data_in_fasta, k=3, sim_or_ident=True)
#print(df['position_in_fasta'], df['similarity'])
#print(df)

# Test of fingerprint search
path_to_structure = str(Path(root) / 'Drugbank_extracted' / 'SDF_ideal.sdf')
path_to_sdf_from_drugbank =  str(Path(root) / 'Drugbank_extracted' / 'structures.sdf') #str(Path(root) / 'Drugbank_extracted' / 'structures_approved_by_db.sdf')#
path_to_sdf_approved = extract_approved_sdf(path_to_sdf_from_drugbank, root, overwrite=True)

get_closest_ligands_from_3d_structure(path_to_structure, path_to_sdf_approved, root,
                                                          fptype='maccs', number_to_print=5)

pdb_dir = Path(root) / 'pdb'
#pdb1 = '1AZM'
#pdb2 = '3W6H'
pdb1 = '3L4Y'
pdb2 = '3L4Z'
aux.download_pdb(pdb1, pdb_dir)
aux.download_pdb(pdb2, pdb_dir)
struct1_path = str(pdb_dir / (pdb1 + '.pdb'))
struct2_path = str(pdb_dir / (pdb2 + '.pdb')) 
# Test of target structure search 
#get_best_pdb_of_target(uniprot, root, path_to_save=None)
#get_best_pdb_of_target(uniprot, root)
# Test of complex structure search by TM-align
#a = get_TMscore_and_RMSD_of_proteins_or_complexes('/home/anton_maximov/BACHELOR/pdb/3W2O.pdb', 
 #                                                 struct2_path, compare_type='c', verbose=True)
#print(a)

# Other auxiliary functions
#get_seq_from_fasta_uniprot_or_seq('P08100')
#get_sequences_similarity('P08100', 'P32238', align_matrix='blosum62', verbose=True)
#get_element_of_fasta_by_number(path_to_data_in_fasta, 1)
