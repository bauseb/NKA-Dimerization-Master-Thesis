
import mdtraj as md
import numpy as np

def mdtraj_contact_calc(file,  cutoff = 0.5, scheme='closest-heavy'):
    #file = 'Structures/7e1z.pdb'
    
    traj = md.load(file)
    top = traj.topology
    
    chain_1_atoms = traj.topology.select('chainid == 1')    
    chain_1_traj = traj.atom_slice(chain_1_atoms)
    
    
    cont = md.compute_contacts(chain_1_traj, scheme = scheme)
    
    cont_distances = cont[0][0]
    
    
    mask_dist = cont_distances < cutoff
    #print(mask_dist)
    ids = np.where(mask_dist)
    cont_idx = []
    for i in ids:
        cont_idx.append(cont[1][i])
        
    cont_idx = np.array(cont_idx)[0]
    #print(cont_idx)
    
    return cont_idx
