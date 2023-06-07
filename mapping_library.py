#import CoEv_To_PDB_func as co
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import itertools
import pandas as pd
import mdtraj as md
from Bio import SeqIO

def split_seq(seq):
    aminoacids = []
    for aa in seq:
        aminoacids.append(aa)
    
    aminoacids_array = np.array(aminoacids)
    return aminoacids_array

def split_cont(cont):
    return cont[0], cont[1]

def coev_to_pfam(prev_cont, V):
    
    cont_row, cont_col = split_cont(prev_cont)
    
    cont_row_pfam = np.zeros(len(cont_row))
    cont_col_pfam = np.zeros(len(cont_col))
    
    for i in range(len(cont_row)):
        cont_row_pfam[i] = V[cont_row[i]]
        cont_col_pfam[i] = V[cont_col[i]]
    
    return cont_row_pfam, cont_col_pfam

def create_seq_from_pdb_structure(file, chainID):
    # Load PDB file
    traj = md.load(file)

    top = traj.top

    atoms, bonds = top.to_dataframe()
    if chainID == "A":
        beta_mask = atoms["chainID"] == 0
    elif chainID == 'B':
        beta_mask = atoms["chainID"] == 1

    atoms_beta = atoms[beta_mask]

    ca_mask = atoms_beta['name'] == 'CA'

    atoms_beta_ca = atoms_beta[ca_mask]

    aa_seq = atoms_beta_ca["resName"]

    aa_map = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y', 'HSE':'H', 'HSP':'H'}

    seq_str = ''

    for aa in aa_seq:
        seq_str += aa_map[aa]

    name_str = file[-8:-4] + "_structure"

    seq = Seq(seq_str)

    # Create a SeqRecord object to store the sequence and its metadata
    record = SeqRecord(seq, id=name_str, name=name_str, description="")

    # Write the SeqRecord object to a fasta file using SeqIO.write()
    #/Users/sebi/Desktop/Master_Thesis/Structures/7e1z/Sequences/rcsb_pdb_7E1Z.fasta
    SeqIO.write(record, file[:-9] + "/Sequences/beta/pdb_structure_" + file[-8:-4].upper() + ".fasta", "fasta")
    
    return None

def create_cont_pairs(cont):
    
    pairs = [(a, b) for a, b in itertools.zip_longest(cont[0], cont[1]) if a != -1 and b != -1]
    
    return np.array(pairs)

def filter_coevolution(matrix, threshold, neighbor_dist, display = 'no'):
    
    #filters the coev matrix, via neglecting close neighbors, only looking at one side
    #and applying a threshold for individual contact values
    
    N = matrix.shape[0]
    filtered_matrix = np.zeros((N, N))
    for i in range(N):
        for j in range(i, N):
            if abs(i-j) >= neighbor_dist and matrix[i, j] > threshold:
                filtered_matrix[i, j] = matrix[i, j]
                filtered_matrix[j, i] = matrix[j, i]
                
                
    n_row, n_col = matrix.shape
    il1 = np.tril_indices(n_row)

    filtered_matrix[il1] = 0
    
    indices = np.where(filtered_matrix > 0.0)
    
    if display == 'yes':
        plot_coevolution(matrix)
        plot_coevolution(filtered_matrix)

    return indices

def plot_coevolution(W):
    
    #Plots the coevolution as a HeatMap
    
    plt.imshow(W, cmap='hot', vmin= np.min(W), vmax= np.max(W))
    plt.colorbar()
    plt.xlabel("contact map residue index")
    plt.ylabel("contact map residue index")

    #plt.title('Coevolution Heatmap')
    plt.show()
    return None


def create_mask(seq):
    msa = split_seq(seq.seq)
    mask = np.array(msa != '-')
    cont_msa_short = np.where(mask)[0]
    
    return cont_msa_short

def get_idx_msa(prev_conts, seq):
    cont_msa_short = create_mask(seq)
    new_conts = []

    for cont in prev_conts:

        new_cont = np.array([cont_msa_short[cont[0]], cont_msa_short[cont[1]]])
        new_conts.append(new_cont)

            
    new_conts = np.array(new_conts)
    return new_conts

def transfer_other_seq_msa(conts, seq):
    mask = create_mask(seq)
    new_conts = []
    
    for cont in conts:
        if cont[0] in mask and cont[1] in mask:
            new_cont = np.array([int(np.where(cont[0] == mask)[0]), int(np.where(cont[1] == mask)[0])])
            new_conts.append(new_cont)
    
    new_conts = np.array(new_conts)
    return new_conts

def test(df, seq):
    sum = 0
    
    df_new = create_empty_df()
    
    contacts = df["contact point"]
    mask = create_mask(seq)
    for i, cont in contacts.items():
        if cont[0] in mask and cont[1] in mask:
            cont_point = np.array([int(np.where(mask == cont[0])[0]), int(np.where(mask == cont[1])[0])])
            df_new = add_element_to_df(df_new, cont_point, '-', df['orientation'][i])
        else:
            sum += 1

    return df_new

def add_aa(df, seq):
  
    #for i, cont in df["contact point"].items():
    i = 0
    for cont in df["contact point"]:
        #print(cont)
        nth_row = i

        if seq.seq[int(cont[0])] != '-' and seq.seq[int(cont[1])] != '-':
            
            aa_pair = np.array([seq.seq[int(cont[0])], seq.seq[int(cont[1])]])
            df.iloc[nth_row]["contact aa"] = aa_pair
        else:

            df.iloc[nth_row]["contact aa"] = 'ERROR'
        i += 1
        

    return df

def transfer_to_orig_seq(prev_cont_df, long_seq):
    #print(prev_cont_df)
    
    
    new_cont_df = create_empty_df()
    
    long_seq_split = split_seq(long_seq.seq)
    mask = np.array(long_seq_split != '-')
    long_seq_aa_idx = np.where(mask)[0]
    short_seq = long_seq_split[mask]

    
    # Convert sequence data to string
    seq_str = "".join(short_seq)
    
    # Create a Seq object from the sequence string
    seq = Seq(seq_str)
    
    # Create a SeqRecord object with the Seq object and ID
    seq_record = SeqRecord(seq, id=long_seq.id)
    
    sum = 0
    i = 0

    for i, prev_cont in prev_cont_df["contact point"].items():
        sum += 1
        #print(sum)

        if isinstance(prev_cont_df["contact aa"][i], np.ndarray):

            #print(prev_cont)
            #print(np.where(long_seq_aa_idx == prev_cont[0])[0])
            #print(type(np.where(long_seq_aa_idx == prev_cont[0])[0]))
            
            
            new_cont = np.array([int(np.where(long_seq_aa_idx == prev_cont[0])[0]), int(np.where(long_seq_aa_idx == prev_cont[1])[0])])
            new_cont_df = add_element_to_df(new_cont_df, new_cont, '-', prev_cont_df["orientation"][i])
            #print(new_cont_df)
            #print("\n")
    new_cont_df = add_aa(new_cont_df, seq_record)
    
    #print(new_cont_df)
    return new_cont_df

def filter_errors(prev_cont_df, transfer_seq):
    
    new_cont_df = create_empty_df()
    
    for i, prev_cont in prev_cont_df["contact point"].items():
        
        if transfer_seq[prev_cont[0]] == '-' or transfer_seq[prev_cont[1]] == '-':
            new_cont_df = add_element_to_df(new_cont_df, prev_cont, 'ERROR', 'ERROR')
        else:
            new_cont_df = add_element_to_df(new_cont_df, prev_cont, '-', prev_cont_df["orientation"][i])
    
    new_cont_df = add_aa(new_cont_df, transfer_seq)
    #print(new_cont_df)
    return new_cont_df

def transfer_to_msa(prev_cont_df, input_seq):
    
    new_cont_df = create_empty_df()
    
    seq_split = split_seq(input_seq.seq)
    mask = np.array(seq_split != '-')
    msa_seq_aa_idx = np.where(mask)[0]
    #print(msa_seq_aa_idx)
    
    for i, prev_cont in prev_cont_df["contact point"].items():
        
        new_cont = np.array([msa_seq_aa_idx[int(prev_cont[0])], msa_seq_aa_idx[int(prev_cont[1])]])
        new_cont_df = add_element_to_df(new_cont_df, new_cont, '-', prev_cont_df["orientation"][i])
    
    new_cont_df = add_aa(new_cont_df, input_seq)
    #print(new_cont_df)
    return new_cont_df

def create_empty_df():
    return pd.DataFrame({
        "contact point": [],
        "contact aa": [],
        "orientation": []
    })

def add_element_to_df(df, new_contact_point, new_contact_aa, new_orientation):
    new_row = pd.DataFrame({
        "contact point": [new_contact_point],
        "contact aa": [new_contact_aa],
        "orientation": [new_orientation]
    })
    #ret = df.append(new_row)
    return pd.concat([df, new_row], ignore_index=True)


def create_cont_dict(conts):
    
    cont_dict = create_empty_df()
    
    for cont in conts:
        
        #cont_aa = get_aa(seq, cont)
        
        cont_dict = add_element_to_df(cont_dict, np.array(cont), '-', 'inter')
        
    return cont_dict

def extract_seq(file_dir):
    
    records = list(SeqIO.parse(file_dir, "fasta"))
    
    return records

def find_gof_pdb(pdb):
    #print(pdb)
    with open(pdb, 'r') as file:
        gof = ''
        for line in file:
            
            # Check if the line starts with 'DBREF'
            if line.startswith('DBREF'):
                # Do something with the line
                #print(line)
                # Split the line into words
                words = line.split()
                #print(words)
                # Check if the 3rd word is 'B'
                if len(words) > 2 and words[2] == 'B':
                    # Extract the 7th word

                    gof = words[6]
                    # Do something with the seventh word
                    #print(gof)
    return gof

def find_gof_pfam(pfam, pdb):
    # Define the input file name and the search term
    input_file = pfam
    #print(pdb)
    search_term = find_gof_pdb(pdb)
    #print(search_term)

    # Iterate over the records in the input file
    abc = 0
    for record in SeqIO.parse(input_file, 'fasta'):
        abc += 1
        # Check if the search term appears in the record ID
        if search_term in record.id:
            # Do something with the record
            #print(record.id)
            #print(record.seq)

            pfam_seq = record
    

    return pfam_seq



def add_seq_shift_of_pdb(prev_cont_df, seq, contact):
    with open(f"Structures/{contact}/{contact}.pdb", 'r') as f:
        shift = None  # default value
        for line in f:
            if line.startswith('ATOM') and line.split()[4] == 'B':
                shift = int(line.split()[5])
                break  # stop searching after finding the first match
    
    new_cont_df = create_empty_df()
    #print(shift)
    for i, cont in prev_cont_df["contact point"].items():
        new_cont = np.array([cont[0]+shift, cont[1] + shift])
        new_cont_df = add_element_to_df(new_cont_df, new_cont, prev_cont_df["contact aa"][i], prev_cont_df["orientation"][i])

    
    return new_cont_df
