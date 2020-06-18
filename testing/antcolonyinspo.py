"""
Inspired by ant colony
Try, not sure if it will work

So from reading the ACO stuff 
I couldn't exactly figure out how to implement te optimization parameters
But I thought up another method inspired by ACO

So what I thought to do was pick a few aminos within the protein.
Fold the protein randomly for all the aminos chosen starting from the current configuration
For each position where you fold check the score 
Pick the lowest score and fold at that amino 
"""

import numpy as np 
import matplotlib.pyplot as plt 
import random
import math
import copy

# PROTEIN = 'PPHPPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHPHH'
PROTEIN = 'CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC'
# PROTEIN = 'HHPHHHPH'
LENGTH = len(PROTEIN)
ITERATIONS = 100
DIRECTIONS = [[-1, 0], [0, 1], [1, 0], [0, -1]]


def protein_configuration():
    x = [i for i in range(LENGTH)]
    y = [0] * LENGTH
    

    rotations = 0
    lowest_score = 0
    best_x = []
    best_y = []
    scores = []

    score = 0
    while rotations < ITERATIONS:

        backup_x = copy.deepcopy(x)
        backup_y = copy.deepcopy(y)

        old_score = copy.deepcopy(score)

        rotating_aminos = ants()
        x, y, score = fold_protein(x, y, rotating_aminos)

        if double(x, y):
            x = backup_x
            y = backup_y
            continue

        if score > old_score:
            x = backup_x
            y = backup_y

        scores.append(old_score)

        if score < lowest_score:
            best_x = copy.deepcopy(x)
            best_y = copy.deepcopy(y)
            lowest_score = copy.deepcopy(score)

        rotations += 1

    return best_x, best_y, lowest_score, scores

def ants():
    """
    choose multiple point on protein to turn 
    choose like half of the amount of points available
    like the "ants"
    """

    # arbitratry integer just to test
    ant_amount = int(LENGTH/2)
    ants = []

    first_ant = True
    for i in range(ant_amount):

        ant = random.randint(0, LENGTH - 1)
        if first_ant:
            ants.append(ant)
            first_ant = False
            continue

        while ant in ants:
            ant = random.randint(0, LENGTH - 1)

        ants.append(ant)

    return ants        


def double(x, y):
    coordinates = []

    for amino_x, amino_y in zip(x, y):
        if [amino_x, amino_y] in coordinates:
            return True
        coordinates.append([amino_x, amino_y])

    return False


def fold_protein(x, y, rot_aminos):

    found_fold = False

    x_folds = []
    y_folds = []
    all_scores = []
    lowest_score = 0

    start_x = copy.deepcopy(x)
    start_y = copy.deepcopy(y)

    for amino_pos in rot_aminos:

        rotating_amino_x = start_x[amino_pos]
        rotating_amino_y = start_y[amino_pos]

        direction = random.random()

        score = 0        
        first_amino = True
        for i in range(amino_pos + 1, LENGTH):
            relative_x = start_x[i] - rotating_amino_x
            relative_y = start_y[i] - rotating_amino_y

            fold_x = start_x
            fold_y = start_y

            if direction > 0.5:
                # rotate left
                fold_x[i] = rotating_amino_x - relative_y
                fold_y[i] = rotating_amino_y + relative_x                
            else:
                # rotate right
                fold_x[i] = rotating_amino_x + relative_y
                fold_y[i] = rotating_amino_y - relative_x

            if first_amino:
                score = check_score(fold_x, fold_y, rotating_amino_x, rotating_amino_y, fold_x[i], fold_y[i], PROTEIN[i]) 
                first_amino = False

            x = fold_x
            y = fold_y

        x_folds.append(x)
        y_folds.append(y)
        all_scores.append(score)

        if score < lowest_score:
            remember_x = copy.deepcopy(x)
            remember_y = copy.deepcopy(y)
            found_fold = True
            lowest_score = copy.deepcopy(score)

    if not found_fold:
        random_fold = random.randint(0, len(x_folds) - 1)
        remember_x = x_folds[random_fold]
        remember_y = y_folds[random_fold]
        lowest_score = all_scores[random_fold]

    return remember_x, remember_y, lowest_score


def check_score(x, y, previous_x, previous_y, amino_x, amino_y, amino):
    protein_coords = list(zip(x, y))   
    score = 0

    for direction in DIRECTIONS:

        dir_x = direction[0]
        dir_y = direction[1]
        check_x = amino_x + dir_x
        check_y = amino_y + dir_y
        neighbour_coords = (check_x, check_y)

        if neighbour_coords in protein_coords and not (check_x == previous_x and check_y == previous_y):
            # print(check_x, check_y)
            # print(protein_coords.index(neighbour_coords))
            neighbour_index = protein_coords.index(neighbour_coords)
            
            neighbour = PROTEIN[neighbour_index]

            if amino == "H":
                if neighbour == "H" or neighbour == "C":
                    score += -1

            if amino == "C":
                if neighbour == "H":
                    score += -1
                if neighbour == "C":
                    score += -5

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
    # plt.savefig("hillclimber.png")
    plt.show()


if __name__ == "__main__":

    best_x, best_y, lowest_score, scores = protein_configuration()
    visualization(best_x, best_y, lowest_score, scores)