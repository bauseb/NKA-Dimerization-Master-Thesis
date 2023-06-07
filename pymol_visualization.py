import pandas as pd
#import pymol as py
import numpy as np
import ast
import mapping_library as ml
import random


structure = '7e1z'

df = pd.read_csv(f'/Users/sebi/Desktop/Master_Thesis/Structures/Contacts/{structure}_intercontacts_beta_shift_oldW.csv')
df['contact point'] = df['contact point'].apply(lambda x: np.fromstring(x[1:-1], sep=' '))

import subprocess

pymol_path = '/opt/homebrew/Cellar/pymol/2.5.0/bin/pymol'  # Replace this with the path to your PyMOL executable
pdb_path = f'/Users/sebi/Desktop/Master_Thesis/Structures/{structure}/{structure}.pdb' # Replace this with the path to your PDB file

#pymol_command = f'load {pdb_file}' # Create the PyMOL command to load the PDB file

    
    # Command to launch PyMOL and load the PDB file
cmd = [pymol_path, pdb_path]



# Launch PyMOL and load the PDB file
#subprocess.Popen(cmd)

new_df = ml.create_empty_df()
loc1 = []
loc0 = []

def determiner(cont):
    if cont <= 34:
        return "intra cellular"
    elif cont <= 62:
        return "transmembrane"
    else:
        return "extra cellular"
    
    
for cont in df["contact point"]:
    loc0.append(determiner(cont[0]))
    loc1.append(determiner(cont[1]))
    
df["loc0"] = loc0
df["loc1"] = loc1 
location = []
contact_points = []
contact_aas = []
orientation = []
for i, loc in df["loc0"].items():
    if loc == df["loc1"][i]:
        contact_points.append(df["contact point"][i])
        contact_aas.append(df["contact aa"][i])
        location.append(loc)
        orientation.append("inter")

new_df["contact point"] = contact_points
new_df["contact aa"] = contact_aas
new_df["orientation"] = orientation
new_df["location"] = location
     
    
print(new_df)
    
file = open(f"Structures/{structure}/contact_predictions_old.pml", "w")

# Define the RGB color range (0-255)
color_range = range(256)

# Loop through the selections

def rgb_to_hex(rgb):
    r, g, b = rgb
    return f"#{r:02x}{g:02x}{b:02x}"

def divide_tuple_by_256(tup):
    return tuple(num/256 for num in tup)
# Example usage:

colors = ['green', 'blue', 'pink', 'purple', 'chocolate', 'brown', 'carbon', 'orange']


# write some text to the file

file.write("load '/Users/sebi/Desktop/Master_Thesis/Structures/7e1z/7e1z_beta.pdb', beta_1")
file.write("\n")
'''file.write("load '/Users/sebi/Desktop/Master_Thesis/Structures/7e1z/7e1z_beta.pdb', beta_2")
file.write("\n")'''

file.write("create beta_2, beta_1")
file.write("\n")
file.write("translate [0, 0, -70], beta_1")
file.write("\n")
file.write("rotate x, 180, beta_2")
file.write("\n")
file.write("color palegreen, beta_1")
file.write("\n")
file.write("color palecyan, beta_2")
file.write("\n")

for i, cont in new_df["contact point"].items():
    
    # Generate a random RGB color tuple
    color = random.choice(colors)

    file.write(f"select cont_{i}.{int(cont[0])}, beta_1 and resid {int(cont[0])}")
    file.write("\n")
    file.write(f"show licorice, beta_1 and cont_{i}.{int(cont[0])}")
    file.write("\n")
    file.write(f"color {color}, beta_1 and cont_{i}.{int(cont[0])}")
    file.write("\n")
    file.write(f"select cont_{i}.{int(cont[1])}, beta_2 and resid {int(cont[1])}")
    file.write("\n")
    file.write(f"show licorice, beta_2 and cont_{i}.{int(cont[1])}")
    file.write("\n")
    file.write(f"color {color}, beta_2 and cont_{i}.{int(cont[1])}")
    file.write("\n")
    file.write(f"distance line_{int(cont[0])}_{int(cont[1])}, cont_{i}.{int(cont[0])} and name ca, cont_{i}.{int(cont[1])} and name ca")
    file.write("\n")
    file.write(f"hide labels, line_{int(cont[0])}_{int(cont[1])}")
    file.write("\n")
    file.write(f"color {color}, line_{int(cont[0])}_{int(cont[1])}")
    file.write("\n")


# close the file
file.close()
