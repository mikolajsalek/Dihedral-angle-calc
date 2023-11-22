import Bio
import numpy as np
from Bio.PDB import *
#import matplotlib.pyplot as plt
import sys
import matplotlib.pyplot as plt
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from matplotlib.cm import get_cmap

def get_secondary_structure(structure):
    dssp = DSSP(structure[0], pdb_filename, file_type="pdb")
    secondary_structure = [dssp[(chain.id, res.id)][2] for model in structure for chain in model for res in chain if res.id[0] == " "]
    return secondary_structure

pdb_list = PDBList()
print(pdb_list)

arg1 = sys.argv[1] #nazwa pliku

#pobieranie pliku
pdb_id = arg1

pdb_filename = pdb_list.retrieve_pdb_file(pdb_id, pdir="data/PDB_files", file_format="pdb")
structure = Bio.PDB.PDBParser().get_structure(pdb_id, pdb_filename)

residue_list = Selection.unfold_entities(structure, "R")
atm_list = []
res_list = []
phi_angles = []
psi_angles = []

secondary_structure = get_secondary_structure(structure)
print(secondary_structure)


phi = ["N", "CA", "C"]
phi_prev = ["CA"]
psi = ["N", "CA", "C"]
psi_next = ["N"]

phi_selected = []
psi_selected = []

matrix_row = []

for r in residue_list:
    tags = r.get_full_id()
    if tags[3][0] == " ":
        #print(r)
        res_list.append(r)


for i in range(0, len(res_list)):
    for atm in res_list[i]:
        #print(atm)
        atm_list.append(atm)
    #print("----------------")

    if i != 0:
        r = res_list[i - 1]
        for x in r:
            if x.get_id() in phi_prev:
                phi_selected.append(x.get_vector())

    for atom in atm_list:
        if atom.get_id() in phi:
            phi_selected.append(atom.get_vector())
        if atom.get_id() in psi:
            psi_selected.append(atom.get_vector())

    if i != len(res_list) - 1:
        r = res_list[i + 1]
        for x in r:
            if x.get_id() in psi_next:
                psi_selected.append(x.get_vector())

    if len(phi_selected) == 4:
        angle_phi = calc_dihedral(phi_selected[0], phi_selected[1], phi_selected[2], phi_selected[3])
        angle_phi_degree = angle_phi * (180.0 / np.pi)
        matrix_row.append(angle_phi_degree)
        phi_angles.append(angle_phi_degree)

    else:
        matrix_row.append("brak")

    if len(psi_selected) == 4:
        angle_psi = calc_dihedral(psi_selected[0], psi_selected[1], psi_selected[2], psi_selected[3])
        angle_psi_degree = angle_psi * (180.0 / np.pi)
        matrix_row.append(angle_psi_degree)
        psi_angles.append(angle_psi_degree)

    else:
        matrix_row.append("brak")

    print(secondary_structure[i])
    #print(res_list[i].get_resname())
    #print(matrix_row)
    matrix_row.clear()
    phi_selected.clear()
    psi_selected.clear()
    atm_list.clear()


#trzeba dokonczyc sb
plt.figure(figsize=(8, 6))

#color_dict = {'H': 'red', 'G': 'blue', 'T': 'green', 'S': 'purple', 'None': 'gray',}

legend_labels = {}  # Dictionary to store unique legend labels
legend_handles = {} # Dictionary to store unique legend handles

for i in range(len(secondary_structure)-1):
    if secondary_structure[i] == 'H':
        label = 'Alpha helix 4-12'
        color = 'red'
    elif secondary_structure[i] == 'G':
        label = '3-10 Helix'
        color = 'blue'
    elif secondary_structure[i] == 'T':
        label = 'Turn'
        color = 'green'
    elif secondary_structure[i] == 'S':
        label = 'Bend'
        color = 'purple'
    elif secondary_structure[i] == 'E':
        label = 'Strand'
        color = 'purple'
    elif secondary_structure[i] == 'B':
        label = 'Isolated beta-bridge residue'
        color = 'purple'
    elif secondary_structure[i] == 'I':
        label = 'Sheet'
        color = 'purple'
    else:
        label = 'None'
        color = 'gray'

    # Use label and color in scatter and store handles for legend
    handle = plt.scatter(phi_angles[i], psi_angles[i], color=color, s=10)
    legend_labels[label] = None
    legend_handles[label] = handle


plt.legend(legend_handles.values(), legend_labels.keys())
plt.title('Ramachandran Plot')
plt.xlabel('Phi (degrees)')
plt.ylabel('Psi (degrees)')
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.grid(color='gray', linestyle='--', linewidth=0.5)



output_file = 'ramachandran_plot.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
