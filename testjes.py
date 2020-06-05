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
                # sys.exit()
                break
                
            best_moves = lowest_score(moves)
            current_position = best_moves[randint(0, len(best_moves) - 1)]
            total_score += current_position[2]
            
            grid[current_position[0]][current_position[1]] = atom
            log_x.append(current_position[0])
            log_y.append(current_position[1])
      
        proteins[key][2] = grid
        proteins[key].append(total_score) 
        proteins[key].append(log_x)  
        proteins[key].append(log_y)

    return proteins

def possible_moves(x, y, _type, grid):
    """Determines which moves are possible from a given position. 
    If possible function will also determine the points gained by making that move."""

    moves = []
    directions = [[x-1,y], [x,y+1], [x+1,y], [x,y-1]]

    # check if the target tile is still empty, and whether that tile doesn't have 4 filled neighbours already, which would result in the string getting wrappped into itself
    for direction in directions:
        if grid[direction[0]][direction[1]] == 0 and neighbour_count(direction[0], direction[1], grid) < 4:
            # calculate points gained by making this move
            score = points(direction[0], direction[1], x, y, _type, grid)
            # place move and associated score in the moves list
            moves.append([direction[0],direction[1], score])

    return moves

def neighbour_count(x, y, grid):
    """Counts how many non-empty neighbours a given position has."""
    counter = 0

    directions = [[x-1,y], [x,y+1], [x+1,y], [x,y-1]]

    for direction in directions:
        # check whether the neighbouring tile is occupied, if so add one to counter
        if grid[direction[0]][direction[1]] != 0:
            counter += 1

    return counter

def points(x, y, old_x, old_y, _type, grid):
    """Determines how many points are gained by making a move to x, y, coming from old_x, old_y."""
    
    counter = 0

    HH_HC_score = -1
    CC_score = -5
    directions = [[x-1,y], [x,y+1], [x+1,y], [x,y-1]]

    #  add to the score if the target position is occupied by an H or C atom
    if _type == 'H':

        for direction in directions:
            if (grid[direction[0]][direction[1]] == 'H' or grid[direction[0]][direction[1]] == 'C') and not (direction[0] == old_x and direction[1] == old_y):
                counter += HH_HC_score            

    # for C atom differentiate between the cases of neighbouring with a C or H atom
    if _type == 'C':

        for direction in directions:
            if grid[direction[0]][direction[1]] == 'H' and not (direction[0] == old_x and direction[1] == old_y):
                counter += HH_HC_score   

        for direction in directions:
            if grid[direction[0]][direction[1]] == 'C' and not (direction[0] == old_x and direction[1] == old_y):
                counter += CC_score  

    return counter


def lowest_score(moves):
    """Determines which move has te highest score."""

    best_moves = []
    best_score = 0

    for move in moves:
        
        # take the score assiciated with each move
        score = move[2]
    
        # if score is equal to the best add it to the list of best moves
        if score == best_score:
            best_moves.append(move)
        # when a lower score is found, empty the list and place the new score in there
        elif score < best_score:
            best_moves = [move]
            best_score = score 
        
    return best_moves


def display(proteins):
    """Gives a graphical display of the protein string after folding."""

    for key in proteins:

        # make two coordinate list for both tpyes of protein
        red_dots_x = []
        red_dots_y = []
        blue_dots_x = []
        blue_dots_y = []
        yellow_dots_x = []
        yellow_dots_y = []

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

                if grid[i][j] == "C":
                    yellow_dots_x.append(i)
                    yellow_dots_y.append(j)

        # make a colorcoded graph
        plt.plot(red_dots_x, red_dots_y, 'or', label='Hydrophobe')
        plt.plot(blue_dots_x, blue_dots_y,'ob', label='Polar')
        
        if len(yellow_dots_x) != 0:
            plt.plot(yellow_dots_x, yellow_dots_y, 'oy', label='Cysteine')
        plt.plot(proteins[key][4], proteins[key][5],'--')
        plt.title(f'Folded protein string, score: {score}')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        plt.xlim(0, grid_size)
        plt.ylim(0, grid_size)
        plt.show()

        
if __name__ == "__main__":
    p = get_proteins()
    g = make_grid(p)
    f = atom_positions(g)
    display(f)
    