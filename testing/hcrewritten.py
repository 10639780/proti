"""
ACO try from https://arxiv.org/pdf/1309.7690.pdf
"""

import numpy as np 
import matplotlib.pyplot as plt 
import random
import math
import copy
# from helpers import * 
import timeit
# PROTEIN = 'PPHPPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHPHH'
PROTEIN = 'CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC'
LENGTH = len(PROTEIN)
ITERATIONS = 1000
DIRECTIONS = [[-1, 0], [0, 1], [1, 0], [0, -1]]

def aco():

    ant_amount = 0
    ants = []

    while ant_amount < ITERATIONS:
        ant = protein_configuration()
        ants.append(ant)
        ant_amount += 1

    lowest_score = 0
    for a in ants:
        score = a[2]

        if score < lowest_score:
            best_ant_x = a[0]
            best_ant_y = a[1]
            lowest_score = score
            best_ant_scores = a[3]

    return best_ant_x, best_ant_y, lowest_score, best_ant_scores




def protein_configuration():
    x = [i for i in range(LENGTH)]
    y = [0] * LENGTH
    
    configurations = []

    rotations = 0
    lowest_score = 0
    best_x = []
    best_y = []
    scores = []

    score = 0
    while rotations < ITERATIONS:

        backup_x = copy.deepcopy(x)
        backup_y = copy.deepcopy(y)

        # old_score = copy.deepcopy(score)

        rotating_amino = random.randint(0, LENGTH - 1)
        fold_protein(x, y, rotating_amino)

        if double(x, y):
            x = backup_x
            y = backup_y
            continue

        configurations.append([x, y])
        
        rotations += 1

    lowest_score = 0
    for c in configurations:
        x = c[0]
        y = c[1]

        score = score_func(x, y)
        scores.append(score)
        if score < lowest_score:
            best_x = x
            best_y = y
            lowest_score = score   

    return best_x, best_y, lowest_score, scores


def double(x, y):
    coordinates = []

    for amino_x, amino_y in zip(x, y):
        if [amino_x, amino_y] in coordinates:
            return True
        coordinates.append([amino_x, amino_y])

    return False


def fold_protein(x, y, amino_pos):

    rotating_amino_x = x[amino_pos]
    rotating_amino_y = y[amino_pos]

    direction = random.random()

    score = 0
    for i in range(amino_pos + 1, LENGTH):
        relative_x = x[i] - rotating_amino_x
        relative_y = y[i] - rotating_amino_y
        
        if direction > 0.5:
            # rotate left
            x[i] = rotating_amino_x - relative_y
            y[i] = rotating_amino_y + relative_x
            
        else:
            # rotate right
            x[i] = rotating_amino_x + relative_y
            y[i] = rotating_amino_y - relative_x

        # if i == amino_pos + 1 and PROTEIN[i] != "P":
        #     score = check_score(x, y, rotating_amino_x, rotating_amino_y, x[i], y[i], PROTEIN[i]) 

    # return score

def score_func(list_x, list_y):
    """Given the coordinates of a protein string, calculate the score of the shape."""

    # list to place the 'already scored' atoms into
    coordinates = []
    directions = [[-1,0],[0,1],[1,0],[0,-1]]
    score = 0
    # list_x, list_y = direction_to_xy(string)
    length = len(list_x)
    protein = PROTEIN


    for i in range(length):

        # P's dont interact so can skip those cases
        if not protein[i] == 'P':

            # for every atom look around in all 4 directions
            for d in directions:

                # check whether one of the previously placed atoms is in the vicinity and determine the score of the interaction with it and the current atom
                for j in range(len(coordinates)):

                    if [list_x[i] + d[0], list_y[i] + d[1]] == coordinates[j] and not(list_x[i] + d[0] == list_x[i-1] and list_y[i] + d[1] == list_y[i-1]):
                        
                        if protein[i] == 'H':
                            if protein[j] == 'H' or protein[j] == 'C':
                                score += -1

                        if protein[i] == 'C':
                            if protein[j] == 'C':
                                score += -5
                            if protein[j] == 'H':
                                score += -1 

        # place in the list with coordinates
        coordinates.append([list_x[i], list_y[i]])
    
    return score


def visualization(x, y, score, scores):

    H_x = []
    H_y = []
    P_x = []
    P_y = []
    C_x = []
    C_y = []

    iterations = [i for i in range(ITERATIONS)]

    for i in range(LENGTH):
        amino = PROTEIN[i]
        amino_x = x[i]
        amino_y = y[i]

        if amino == 'H':
            H_x.append(amino_x)
            H_y.append(amino_y)
        if amino == 'P':
            P_x.append(amino_x)
            P_y.append(amino_y)
        if amino == 'C':
            C_x.append(amino_x)
            C_y.append(amino_y)


    plt.figure()

    plt.subplot(211)
    plt.plot(x, y, "k--")
    plt.plot(H_x, H_y, "ro", label = "hydrofoob")
    plt.plot(P_x, P_y, "bo", label = "polair")
    plt.plot(C_x, C_y, "go", label = "cysteine")
    plt.xlabel("step")
    plt.ylabel("step")
    plt.legend()
    plt.title("Total score: {}".format(score))

    plt.subplot(212)
    plt.plot(iterations, scores)
    plt.title("Statistics")
    plt.xlabel("iterations")
    plt.ylabel("scores")
    plt.show()


if __name__ == "__main__":
    # protein_configuration()
    # best_x, best_y, lowest_score, scores = protein_configuration()
    # visualization(best_x, best_y, lowest_score, scores)
    # print(aco())
    best_x, best_y, lowest_score, scores = aco()
    visualization(best_x, best_y, lowest_score, scores)