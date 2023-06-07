# NKA-Dimerization-Master-Thesis
This repository contains all the code written by me for my master thesis "Co-evolutional anaylsis of the Na+,K+-ATPase’s β-subunit dimerization"

It contains 4 python scripts:

1.) global_function.py: This is the main mapping and filtering code. It takes as an input the contact map matrix, the MSA, the MSA indices used for the contact map and different structure/sequence files of the various used structures.

The required structure/sequence files for one used protein are:
  1. .pdb file
  2. .fasta file of the gene in UniProt
  3. .fasta file of the sequence generated from the pdb file
  4. .fasta file of an MSA between the UniProt- and the pdb-sequence file
  5. .fasta file of an MSA between the UniProt- and the original MSA-sequence file

2.) mapping_library.py: This file contains all the written functions, which are used in "global_function.py".

3.) mdtraj_contacts.py: This script contains the function to compute the intra contacts of a protein, which is used to differentiate between intra- and inter-contacts.

4.) pymol_visualization.py: With this file, a pymol script is genereated, that visualizes all predicted inter contacts on two beta sububits of a given pdb.
