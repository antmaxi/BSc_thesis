# BSc_thesis
A screening protocol to find possible off-targets and search for polypharmacology using different similarity measures and comparison to approved by FDA drugs and correspendent targets or just some proteomes.

See *Installation* for the list of requirements.

# Files description

**Search.py** - main file, there are the functions to calculate similarity between: 
        1. two fasta sequences (using BioPython)
        2. one fasta and the whole database targets (using BioPython)
        3. SMILES and list of SMILES (using RDkit)
        4. SMILES and the whole database ligands (using RDkit)
        5. ligand and the whole database (using pybel fingerprints)
        6. two protein structure files (TM-score and RMSD, using TM-align) 
        
**Auxiliaray.py** - Some auxillary functions (download files and pdbs, load processed from Drugbank data) 
and functions for jumping between ids:
Uniprot => a/a sequence
Uniprot => PDBs
SMILES => PDBs
PubCHEM => SMILES
Name of ligand => Uniprots of targets
Name of ligand => SMILES
Also function to convert different structural types (pdb, sdf, mol2) between themselves

**Drugbank.py** - Functions to get needed info from Drugbank, save it. Possible to load info and update the database

**IsoMIF.py** - Using 2 pdb ids and numbers of two clefts from the biggest gets tanimoto coefficient between them
