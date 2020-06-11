"""
Completely random configuration without score rewritten and tested
"""

import matplotlib.pyplot as plt
import numpy as np
from random import randint

PROTEIN = 'HHPHHHPH'
LENGTH = len(PROTEIN)
DIMENSION = 2
GRID_SIZE = 2 * LENGTH + 1
DIRECTIONS = [[1, 0], [-1, 0], [0, -1], [0, 1]]

def init_grid():
    grid = [[0 for i in range(GRID_SIZE)] for i in range(GRID_SIZE)]
    
    return grid

def display_grid(grid):    
    matrix = "\n".join([" ".join(['{} '.format(element) for element in row]) for row in grid])
    
    return matrix

def protein_positions():
    grid = init_grid()

    first_amino = True
    pos_x = [LENGTH]
    pos_y = [LENGTH]

    for amino in range(LENGTH):
        
        if first_amino:
            grid[pos_x[amino]][pos_y[amino]] = PROTEIN[amino]
            first_amino = False
            continue

        x = pos_x[amino - 1]
        y = pos_y[amino - 1]

        new_x, new_y = protein_fold(x, y)
        while not check_grid_space(new_x, new_y, grid):
            new_x, new_y = protein_fold(x, y)
        
        grid[new_x][new_y] = PROTEIN[amino]
        pos_x.append(new_x)
        pos_y.append(new_y)

    return grid, pos_x, pos_y

def protein_fold(x, y):
    ran_direct = randint(0, len(DIRECTIONS) - 1)
    direction = DIRECTIONS[ran_direct]
    
    new_x = x - direction[1]
    new_y = y + direction[0]

    return new_x, new_y

def check_grid_space(new_x, new_y, grid):
    potential_position = grid[new_x][new_y]

    if potential_position != 0:
        return False

    return True 


def plot_protein_config(x_list, y_list):
    H_x = []
    H_y = []
    P_x = []
    P_y = []

    for i in range(LENGTH):
        amino = PROTEIN[i]
        x = x_list[i]
        y = y_list[i]

        if amino == 'H':
            H_x.append(x)
            H_y.append(y)
        else:
            P_x.append(x)
            P_y.append(y)

    plt.plot(H_x, H_y, "ro", label="hydrofoob")
    plt.plot(P_x, P_y, "bo", label="polair")
    plt.plot(x_list, y_list, "k-")
    plt.legend()
    plt.show()

        

if __name__ == "__main__":
    protein_config = protein_positions()
    grid = display_grid(protein_config[0])
    plot_protein_config(protein_config[1], protein_config[2])
    print(grid)
