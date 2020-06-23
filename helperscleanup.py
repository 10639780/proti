"""
helpers.py

Minor Programmeren 
Team Proti

Algorithms are implemented by calling the functions 
in this script.
"""

# import modules
import matplotlib.pyplot as plt
import random
from classes import protein
from mpl_toolkits import mplot3d
from anytree import Node, RenderTree, Walker, PreOrderIter
import operator
import copy
import math
from classes.amino import *
from classes.protein import *
from helpersrepeated import * 

def plot(proti, score, list_x, list_y, list_z=False, scores=False):

    H_coords, P_coords, C_coords = plot_coords(proti, list_x, list_y, list_z)

    if len(H_coords[-1]) == 0:

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9))

        ax1.plot(H_coords[0], H_coords[1], 'or')
        ax1.plot(P_coords[0], P_coords[1], 'ob')
        ax1.plot(C_coords[0], C_coords[1], 'oy')
        ax1.set_title(f'Folded protein, score: {score}')
        ax1.axis('equal')

        ax2.plot(scores)
        ax2.set_title('Scores of the configurations after each rotation')
        ax2.set(xlabel='Rotation', ylabel='Score')

    else:
        fig = plt.figure()

        ax = fig.add_subplot(2,1, 1, projection='3d')
        ax.plot(H_coords[0], H_coords[1], H_coords[2], 'or')
        ax.plot(P_coords[0], P_coords[1], P_coords[2], 'ob')
        ax.plot(C_coords[0], C_coords[1], C_coords[2], 'oy')
        ax.set_title(f'Folded protein, score: {score}')
        
        ax = fig.add_subplot(2, 1, 2)
        ax.plot(scores)
        ax.set_title('Scores of the configurations after each rotation')
        ax.set(xlabel='Rotation', ylabel='Score')

    plt.show()

def score(proti, list_x, list_y, list_z=False):

    if not list_z:
        list_z = [0] * proti.length

    coordinates = []
    directions = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
    score = 0 

    for i in range(proti.length):
        amino_type = proti.listed[i]
        # P's dont interact so can skip those cases
        if not amino_type == 'P':
            current_x = list_x[i]
            current_y = list_y[i]
            current_z = list_z[i]

            prev_x = list_x[i - 1]
            prev_y = list_y[i - 1]
            prev_z = list_z[i - 1]
            
            # for every atom look around in all 6 directions
            for d in directions:

                neighbour_x = current_x + d[0]
                neighbour_y = current_y + d[1]
                neighbour_z = current_z + d[2]

                # check whether one of the previously placed atoms is in the vicinity 
                # and determine the score of the interaction with it and the current atom
                for j in range(len(coordinates)):
                    neighbour_amino = proti.listed[j]

                    if [neighbour_x, neighbour_y, neighbour_z] == coordinates[j] \
                        and not (neighbour_x == prev_x and neighbour_y == prev_y\
                             and neighbour_z == prev_z):

                        if amino_type == 'H':
                            if neighbour_amino == 'H' or neighbour_amino == 'C':
                                score += -1

                        if amino_type == 'C':
                            if neighbour_amino == 'C':
                                score += -5
                            if neighbour_amino == 'H':
                                score += -1

        # place in the list with coordinates
        coordinates.append([current_x, current_y, current_z])

    return score


def double(list_x, list_y, list_z=False):
    

    coordinates = []

    if not list_z:
        list_z = [0] * len(list_x)

    for x, y, z in zip(list_x, list_y, list_z):
        if [x, y, z] in coordinates:
            return True
        
        coordinates.append([x, y, z])

    return False

    
def output(proti, score, list_x, list_y, list_z=False):
    """
    Prints the folded string to a csv file in the Bas Terwijn style.
    """

    numbers = []

    for i in range(proti.length - 1):

        # new position is compared to the old
        delta_x = list_x[i + 1] - list_x[i]
        delta_y = list_y[i + 1] - list_y[i]
        delta_z = list_z[i + 1] - list_z[i]

        # conversion between coordinates  and bas terwijn numbers
        if delta_x == 1:
            number = -2
        elif delta_x == -1:
            number = 2
        elif delta_y == 1:
            number = 1
        elif delta_y == -1:
            number = -1

        if list_z:
            delta_z = list_z[i + 1] - list_z[i]
            if delta_z == 1:
                number = 3
            elif delta_z == -1:
                number = -3

        numbers.append(number)

    # add 0 to signal end of protein
    numbers.append(0)

    # write the list to a file
    f = open('output.csv', 'w')
    f.write('amino,fold\n')
    for p, n in zip(proti.listed, numbers):
        f.write(f'{p}, {n}\n')
    f.write(f'score,{score}')
    f.close()
        



if __name__ == "__main__":
    x = [i for i in range(20)]
    y = [i for i in range(20)]
    z = [i for i in range(20)]
    scores = [0] * 20
    bool_score = 0

    proteinstring = 'HPHPPHHPHPPHPHHPPHPH'
    length = len(proteinstring)
    proti = Protein(proteinstring, length)

    # test = plot(proti=proti, score=bool_score, list_x=x, list_y=y)

    # string = 'LRSLRSLRS'
    # test1 = plot(proti=proti, score=bool_score, scores=scores, list_x = x, list_y=y, list_z=z)

    # test_score = score(proti=proti, list_x=x, list_y=y)
    # print(test_score)


