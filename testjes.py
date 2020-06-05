"""
File om testjes te runnen voordat we aanpassingen doen in prot.py
"""

"""
Deze functie zouden we kunnen grebruiken om de proteins uit een file te lezen 
Waardoor we meerdere proteins in 1 keer door de code kunnen sturen (van prot.py)
Ik ga er later mee door om het in prot.py te verwerken
"""

import numpy as np 
from random import randint
import sys
import matplotlib.pyplot as plt 
import csv

def get_proteins():

    """Saves the proteins to analyse in a dictionary"""

    proteins = {}
    filename = "test.txt"

    with open(filename, "r") as file:
        lines = file.readlines()

        for i in range(len(lines)):
            protein = lines[i].strip("\n") 
            length = len(protein)
            proteins[i] = [protein, length] 

    return proteins 

def make_grid(proteins):
    """
    Creates a grid of specified size.
    Saves the grid for all proteins in the dictionary.
    """

    for key in proteins:

        grid = []
        size = 2*proteins[key][1]+1
        for i in range(size):
            row = []
            for j in range(size):
                row.append(0)
            grid.append(row)

        proteins[key].append(grid)

    return proteins

def atom_positions(proteins):
    
    for key in proteins:

        protein = proteins[key][0]
        length = proteins[key][1]
        grid = proteins[key][2]

        current_position = [length, length, 0]
        log_x = []
        log_y = []

        total_score = 0
        first_loop = True

        for atom in protein:
           
            if first_loop:
                grid[length][length] = atom
                log_x.append(length)
                log_y.append(length)
                first_loop = False

                continue

            moves = possible_moves(current_position[0], current_position[0], atom, grid)
            
            if len(moves) == 0:
                print("Protein is stuck.")
                # would like to change this so if 1 protein stuck it does continue for the other proteins in the file
                sys.exit()
                

            current_position = moves[randint(0, len(moves) - 1)]
            total_score += current_position[2]
            
            grid[current_position[0]][current_position[1]] = atom
            log_x.append(current_position[0])
            log_y.append(current_position[1])
      
        proteins[key][2] = grid
        proteins[key].append(total_score)       

    return proteins

def possible_moves(x, y, _type, grid):
    """Determines which moves are possible from a given position. 
    If possible function will also determine the points gained by making that move."""

    # list to hold all valid moves
    moves = []
   
    # check if the target tile is still empty, and whether that tile doesn't have 4 filled neighbours already, which would result in the string getting wrappped into itself
    if grid[x-1][y] == 0 and neighbour_count(x-1, y, grid) < 4:
        score = 0
        # if the protein to be placed is hydrophobe, determine what the score would be if placed in this position
        if _type == 'H':
            score = points(x-1, y, x, y, grid)
        # add to list of valid moves
        moves.append([x-1, y, score])
    
    if grid[x][y+1] == 0 and neighbour_count(x, y+1, grid) < 4:
        score = 0
        if _type == 'H':
            score = points(x, y+1, x, y, grid)
        moves.append([x, y+1, score])

    if grid[x+1][y] == 0 and neighbour_count(x+1, y, grid) < 4:
        score = 0
        if _type == 'H':
            score = points(x+1, y, x, y, grid)
        moves.append([x+1, y, score])

    if grid[x][y-1] == 0 and neighbour_count(x, y-1, grid) < 4:
        score = 0
        if _type == 'H':
            score = points(x, y-1, x, y, grid)
        moves.append([x, y-1, score])

    return moves

def neighbour_count(x, y, grid):
    """Counts how many non-empty neighbours a given position has."""
    counter = 0

    # goes to all 4 neighbouring tiles, if not empty, add 1 to counter
    if grid[x-1][y] != 0:
        counter += 1
    
    if grid[x][y+1] != 0:
        counter += 1

    if grid[x+1][y] != 0:
        counter += 1

    if grid[x][y-1] != 0:
        counter += 1

    return counter

def points(x, y, old_x, old_y, grid):
    """Determines how many points are gained by making a move to x, y, coming from old_x, old_y."""
    
    counter = 0

    # add 1 point if the target position is a hydrophobe and is not the previous protein in the string
    if grid[x-1][y] == 'H' and not (x-1 == old_x and y == old_y):
        counter += 1
    
    if grid[x][y+1] == 'H' and not (x == old_x and y+1 == old_y):
        counter += 1

    if grid[x+1][y] == 'H' and not (x+1 == old_x and y == old_y):
        counter += 1

    if grid[x][y-1] == 'H' and not (x == old_x and y-1 == old_y):
        counter += 1

    return counter

def display(proteins):
    """Gives a graphical display of the protein string after folding."""

    for key in proteins:

        # make two coordinate list for both tpyes of protein
        red_dots_x = []
        red_dots_y = []
        blue_dots_x = []
        blue_dots_y = []

        grid_size = 2 * proteins[key][1] + 1
        grid = proteins[key][2]
        score = proteins[key][3]

        # iterate through the grid and add coordinates of the proteins found
        for i in range(grid_size):
            for j in range(grid_size):

                if grid[i][j] == 'H':
                    red_dots_x.append(i)
                    red_dots_y.append(j)
                
                if grid[i][j] == 'P':
                    blue_dots_x.append(i)
                    blue_dots_y.append(j)

        # make a colorcoded graph
        plt.plot(red_dots_x,red_dots_y,'or',blue_dots_x,blue_dots_y,'ob')
        plt.title('Folded protein string')
        plt.text(1, 1, f'Score: {score}')
        plt.xlim(0, grid_size)
        plt.ylim(0, grid_size)
        plt.show()

        
if __name__ == "__main__":
    p = get_proteins()
    g = make_grid(p)
    f = atom_positions(g)
    display(f)
    