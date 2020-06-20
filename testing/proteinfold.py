"""
Minor programmeren

Team proti

Case Protein Pow(d)er

Just testing some fold begin stages random algo
"""
import numpy as np
import matplotlib.pyplot as plt
from random import randint

def protein():
    """
    Returns the protein as a string
    Later to be updated to read in all proteins
    """
    protein = 'HHPHHHPH'

    return protein

def folding(protein):
    """
    Folds protein? figure it out
    """
    
    length = len(protein)
   
    # start in the center (can be any random point but middle is just convention)
    start_position = (length, length)
    start_atom = protein[0]
    # will be a list of the coordinates of each atom in the form of tuples       
    atom_coords = [start_position]
    
    first_step = True
    
    for i in range(len(protein) - 1):
        
        next_atom = protein[i + 1]
        pre_coords = atom_coords[i]        
        step = atom_step()

        potential_move = check_step(step, pre_coords, atom_coords)

        while potential_move in atom_coords:
            new_step = atom_step()
            potential_move = check_step(new_step, pre_coords, atom_coords)
            
        if first_step:                      
            
            next_coords = potential_move
            atom_coords.append(next_coords)

            first_step = False
            continue

        next_coords = potential_move
        atom_coords.append(next_coords)

    print(atom_coords)
    return atom_coords
                      
def atom_step():

    directions = randint(-2, 2)

    while directions == 0:
        directions = randint(-2, 2)

    return directions 


def check_step(step, coords, coords_list):
    
    if step == -1 or step == 1:
        potential_coords = (coords[0] + step, coords[1])

    else:
        potential_coords = (coords[0], coords[1] + step)

    # if potential_coords in coords_list:
    #     print("Protein folded on itself")
    #     return

    return potential_coords

def plot_fold(coords):
    plt.scatter(*zip(*coords))
    plt.show()

if __name__ == "__main__":
    protein = protein()
    fold = folding(protein)
    plot_fold(fold)


