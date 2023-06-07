import os
import mdtraj_contacts as mdc
import numpy as np
import mapping_library as ml
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# specify the folder path
folder_path = 'enter_your_folder_path'

W = np.load("enter_your_contact_map_matrix") #co-evolution score matrix
V_idx = np.load("enter_your_array_of_used_indices_in_the_MSA") #indices of used columns of the MSA
msa_pfam = "enter_your_MSA_file" #MSA

#Define Threshold and shortest distance between neighbors
thresh = 3
dist = 3
sum = 0
cutoff = 0.5

#clear and filter W
cont_coev = ml.filter_coevolution(W, thresh, dist, 'yes')

#transfer to V, which corresponds to contacts in the PFAM MSA
cont_pfam = ml.create_cont_pairs(ml.coev_to_pfam(cont_coev, V_idx))
contacts_df = ml.create_cont_dict(cont_pfam)

# iterate through all folders in the folder
for folder in os.listdir(folder_path):

    folder_abs_path = os.path.join(folder_path, folder)
    
    if os.path.isdir(folder_abs_path) and "Uniprot" not in folder_abs_path and "Contacts" not in folder_abs_path:
        #print(folder_abs_path)
        # do something with the folder

        pdb_structure_file = folder_abs_path + "/" + folder_abs_path[-4:] + ".pdb"
        pdb_fasta_file = ml.extract_seq(folder_abs_path + "/Sequences/beta/rcsb_pdb_" + folder_abs_path[-4:].upper() + ".fasta")
        msa_pfam_pdb_file = ml.extract_seq(folder_abs_path + "/Sequences/beta/msa_pfam_pdb_" + folder_abs_path[-4:].upper() + ".fasta")
        ml.create_seq_from_pdb_structure(pdb_structure_file, "B")
        msa_pdbf_pdbs_file = ml.extract_seq(folder_abs_path + '/Sequences/beta/msa_pdbf_pdbs_' + folder_abs_path[-4:].upper() + ".fasta")
        
        gene_of_interest = ml.find_gof_pfam(msa_pfam, pdb_structure_file)

        print("PDB structure: " + str(folder_abs_path[-4:]))
        intra_cont = mdc.mdtraj_contact_calc(pdb_structure_file, cutoff = cutoff, scheme = "closest")
        
        for seq in msa_pdbf_pdbs_file:
            if 'structure' in seq.id:
                pdbs = seq
            else:
                pdbf = seq
        
        msa_pdbs_conts = ml.get_idx_msa(intra_cont, pdbs)
        msa_pdbf_conts = ml.transfer_other_seq_msa(msa_pdbs_conts, pdbf)

        for seq in msa_pfam_pdb_file:
            if gene_of_interest.id in seq.id:
                pfam = seq
            else:
                pdb = seq

        msa_pfam_pdb_conts = ml.get_idx_msa(msa_pdbf_conts, pdb)

        gof_conts = ml.transfer_other_seq_msa(msa_pfam_pdb_conts, pfam)
        pfam_conts = ml.get_idx_msa(gof_conts, gene_of_interest)
            
        for i, cont in enumerate(contacts_df["contact point"]):

            if np.any(np.all(pfam_conts == cont, axis=1)):
                
                nth_row = i
                df = contacts_df.iloc[nth_row]
                
                if df['orientation'] != 'INTRA':
                    sum += 1
                    contacts_df.iloc[nth_row]['orientation'] = 'INTRA'
        print("Current amount of found INTRA contact points: " + str(sum) + "\n")
        
print(contacts_df)
print(str(sum) + " INTRA contacts were found. This leaves " + str(len(contacts_df) - sum) + " possible inter contacts with a cutoff distance of " + str(cutoff) + " nm and a coevoultion threshold of " + str(thresh) + ".")

#co.make_dist_hist(contacts_df)     

contacts_df.to_csv('Structures/Contacts/contacts.csv', index=False)

#########################################
#########################################
'FIND CONTACTS IN ONE SPECIFIC STRUCTURE'
#########################################
#########################################

contact = "7e1z"
contact_upper = contact.upper()

msa_pfam = "/Users/sebi/Desktop/Master_Thesis/Sequences/NaK_ATPase_PFAM.fasta"
pdb = f"Structures/{contact}/{contact}.pdb"
pdbs = ml.extract_seq(f"Structures/{contact}/Sequences/beta/pdb_structure_{contact_upper}.fasta")
pdbf = ml.extract_seq(f"Structures/{contact}/Sequences/beta/rcsb_pdb_{contact_upper}.fasta")[1]
gof = ml.find_gof_pfam(msa_pfam, pdb)
msa_uniprot = ml.extract_seq(f'/Users/sebi/Desktop/Master_Thesis/Structures/{contact}/Sequences/beta/msa_pfam_pdb_{contact_upper}.fasta')

pfam = msa_uniprot[0]
pdb = msa_uniprot[1]

msa_pdb = ml.extract_seq(f"Structures/{contact}/Sequences/beta/msa_pdbf_pdbs_{contact_upper}.fasta")

for seq in msa_pdb:
    #print(seq)
    #print("\n")
    if 'Chain' in seq.id:
        pdbf_msa = seq
    else:
        pdbs_msa = seq

###### PRINT CONTACTS IN BROAD PFAM SEQ #####
starting_contacts = ml.add_aa(contacts_df, gof)
print("STARTING CONT")
print(starting_contacts)

uniprot_contacts = ml.transfer_to_orig_seq(starting_contacts, gof)
uniprot_pdb_msa_contacts = ml.transfer_to_msa(uniprot_contacts, pfam)
pdb_uniprot_msa_contacts = ml.filter_errors(uniprot_pdb_msa_contacts, pdb)
a = ml.transfer_to_orig_seq(pdb_uniprot_msa_contacts, pdbf)
b = ml.transfer_to_msa(a, pdbf_msa)
c = ml.filter_errors(b, pdbs_msa)
d = ml.transfer_to_orig_seq(c, pdbs_msa)     
e = d[d['orientation'] == 'inter']
inter_contacts = e

inter_contacts_shift = ml.add_seq_shift_of_pdb(inter_contacts, pdb, contact)
print("INTER CONTACTS SHIFT TO UNIPROT NUMBERS")
print(inter_contacts_shift)

inter_contacts_shift.to_csv(f'{contact}_intercontacts.csv', index=False) #outputfile
