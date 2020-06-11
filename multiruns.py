"""
Run semi random algorithm N times and pick the configuration with the highest score
From this I think we can modify to a hillclimber
"""

import matplotlib.pyplot as plt
import numpy as np
from random import randint
import sys
import copy

#PROTEIN = 'HHPHHHPH'
PROTEIN ='PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP'
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

    total_score = 0
    for amino in range(LENGTH):
        
        element = PROTEIN[amino]
        if first_amino:
            grid[pos_x[amino]][pos_y[amino]] = element
            first_amino = False
            continue

        x = pos_x[amino - 1]
        y = pos_y[amino - 1]

        step = score_valid_folds(x, y, grid, element)

        if not step:
            # print("Skipping stuck protein")
            return False

        new_x = step[0]
        new_y = step[1]
        score = step[2]

        grid[new_x][new_y] = element
        pos_x.append(new_x)
        pos_y.append(new_y)
        total_score += score

    return grid, pos_x, pos_y, total_score

def all_possible_folds(x, y):

    possible_folds = []

    for direction in DIRECTIONS:
        new_x = x - direction[1]
        new_y = y + direction[0]
        possible_folds.append([new_x, new_y])
    
    return possible_folds

def check_valid_folds(x, y, grid):

    possible_folds = all_possible_folds(x, y)
    valid_folds = []
    for fold in possible_folds:
        new_x = fold[0]
        new_y = fold[1]
        value = grid[new_x][new_y]
        if value == 0:
            valid_folds.append(fold)
            

    return valid_folds

def score_valid_folds(x, y, grid, amino):
    valid_folds = check_valid_folds(x, y, grid)

    if len(valid_folds) == 0:   
        # print("Protein got stuck")     
        return False

    fold_neighbours = []
    
    for fold in valid_folds:
        score = 0
        new_x = fold[0]
        new_y = fold[1]
        local_score = 0
        neighbours = []
        if [new_x - 1, new_y] != [x, y]:
            grid_value = grid[new_x - 1][new_y]
            neighbours.append([new_x - 1, new_y, grid_value])

        if [new_x + 1, new_y] != [x, y]:
            grid_value = grid[new_x + 1][new_y]
            neighbours.append([new_x + 1, new_y, grid_value])

        if [new_x, new_y - 1] != [x, y]:
            grid_value = grid[new_x][new_y - 1]
            neighbours.append([new_x, new_y - 1, grid_value])

        if [new_x, new_y + 1] != [x, y]:
            grid_value = grid[new_x][new_y + 1]
            neighbours.append([new_x, new_y + 1, grid_value])

        for neighbour in neighbours:
            score = 0
            grid_value = neighbour[2]
            
            if grid_value == 'H' and (amino == 'H' or amino == 'C'):
                score = -1
            if grid_value == 'C' and amino == 'H':
                score = -1
            if grid_value == 'C' and amino == "C":
                score = -5

            neighbour.append(score)

        fold_neighbours.append([fold, neighbours])

    lowest_score = 0
    first_best_move = []

    for fold_info in fold_neighbours:
        fold = fold_info[0]
        neighbours = fold_info[1]
                
        for neighbour in neighbours:
            score = neighbour[3]
            local_score += score
            if score < lowest_score:
                lowest_score = score
                first_best_move = fold
                new_x = first_best_move[0]
                new_y = first_best_move[1]

    if len(first_best_move) == 0:
        random_choice = randint(0, len(valid_folds) - 1)
        random_move = valid_folds[random_choice]
        new_x = random_move[0]
        new_y = random_move[1]
        
    return new_x, new_y, local_score

def best_score_configuration():
    best_score = 0
    N = 100000

    for i in range(N):
        current_config = protein_positions()

        while not current_config:
            current_config = protein_positions()

        score = current_config[-1]
        
        if best_score > score:
            best_score = score
            best_configuration = copy.deepcopy(current_config)

    grid, x, y, score = best_configuration
    return grid, x, y, score



def plot_protein_config(x_list, y_list, total_score):
    H_x = []
    H_y = []
    P_x = []
    P_y = []
    C_x = []
    C_y = []

    for i in range(LENGTH):
        amino = PROTEIN[i]
        x = x_list[i]
        y = y_list[i]

        if amino == 'H':
            H_x.append(x)
            H_y.append(y)
        elif amino == 'C':
            C_x.append(x)
            C_y.append(y)
        else:
            P_x.append(x)
            P_y.append(y)

    plt.plot(x_list, y_list, "k--")
    plt.plot(H_x, H_y, "ro", label="hydrofoob")
    plt.plot(P_x, P_y, "bo", label="polair")

    if len(C_x) > 0:
        plt.plot(C_x, C_y, "go", label="cysteine")
    
    plt.title("Length of protein: {}. Total score: {}".format(LENGTH, total_score))
    plt.legend()
    plt.show()

        
if __name__ == "__main__":
    grid, x, y, score = best_score_configuration()
    plot_protein_config(x, y, score)
    
