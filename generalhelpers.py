"""
generalhelpers.py

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


def plot(proti, score, list_x, list_y, list_z=False, scores=False):
    """
    Plots a 2D or 3D fold of the protein in space.
    Also displays the scores of the fold for every iteration.
    """

    # splits up the coordinates of aminos types for plot
    H_coords, P_coords, C_coords = plot_coords(proti, list_x, list_y, list_z)

    # no coordinates in z direction so plot 2D
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

    # otherwise plot 3D
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


def plot_coords(proti, x_list, y_list, z_list=False):

    # seperate amino type coordinates in lists
    H_x = []
    H_y = []
    H_z = []

    P_x = []
    P_y = []
    P_z = []

    C_x = []
    C_y = []
    C_z = []

    H_coords = []
    P_coords = []
    C_coords = []

    # check for z direction
    if z_list:
        for x, y, z, p in zip(x_list, y_list, z_list, proti.listed):

            if p == 'H':
                H_x.append(x)
                H_y.append(y)
                H_z.append(z)
            if p == 'P':
                P_x.append(x)
                P_y.append(y)
                P_z.append(z)
            if p == 'C':
                C_x.append(x)
                C_y.append(y)
                C_z.append(z)
    else:
        for x, y, p in zip(x_list, y_list, proti.listed):
            if p == 'H':
                H_x.append(x)
                H_y.append(y)                
            if p == 'P':
                P_x.append(x)
                P_y.append(y)                
            if p == 'C':
                C_x.append(x)
                C_y.append(y)

    # seperate all the coordinate lists per amino type
    H_coords.append(H_x)
    H_coords.append(H_y)
    H_coords.append(H_z)

    
    P_coords.append(P_x)
    P_coords.append(P_y)
    P_coords.append(P_z)

    C_coords.append(C_x)
    C_coords.append(C_y)
    C_coords.append(C_z)

    return H_coords, P_coords, C_coords


def score_it(proti, list_x, list_y, list_z=False):
    """
    Determines the score for iterative folding algorithms
    """

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


def score_con(pos_x, pos_y, proti):
    """
    Determines the score of the protein for constructive folding algorithms
    """

    score = 0
    coordinates = []
    check_coordinates = []

    for x, y in zip(pos_x, pos_y):
        coordinates.append([x, y])

    delta_x = pos_x[-1] - pos_x[-2]
    delta_y = pos_y[-1] - pos_y[-2]

    left_x = pos_x[-1] - delta_y
    left_y = pos_y[-1] + delta_x
    check_coordinates.append([left_x, left_y])

    right_x = pos_x[-1] + delta_y
    right_y = pos_y[-1] - delta_x
    check_coordinates.append([right_x, right_y])

    straight_x = pos_x[-1] + delta_x
    straight_y = pos_y[-1] + delta_y
    check_coordinates.append([straight_x, straight_y])

    for c in check_coordinates:

        if c in coordinates:

            index = coordinates.index(c)

            if proti.listed[len(pos_x) - 1] == 'H' and \
                (proti.listed[index] == 'H' or proti.listed[index] == 'C'):
                score += -1
            elif proti.listed[len(pos_x) - 1] == 'C' and proti.listed[index] == 'H':
                score += -1
            elif proti.listed[len(pos_x) - 1] == 'C' and proti.listed[index] == 'C':
                score += -5

    return score

def directions(string):
    """
    Converts a series of string with directions like ['L', 'R', 'U', 'D'] to
    lists with xyz positions.
    """

    pos_x = [0, 1]
    pos_y = [0, 0]
    pos_z = [0, 0]

    # go over every node
    for s in string:

        # previous direction is determined
        delta_x = pos_x[-1] - pos_x[-2]
        delta_y = pos_y[-1] - pos_y[-2]
        delta_z = pos_z[-1] - pos_z[-2]

        # rotation matrices used to turn into the desired direction
        if s == 'S':
            pos_x.append(pos_x[-1] + delta_x)
            pos_y.append(pos_y[-1] + delta_y)
            pos_z.append(pos_z[-1] + delta_z)

        elif s == 'L':
            pos_x.append(pos_x[-1] - delta_y)
            pos_y.append(pos_y[-1] + delta_x)
            pos_z.append(pos_z[-1] + delta_z)

        elif s == 'R':
            pos_x.append(pos_x[-1] + delta_y)
            pos_y.append(pos_y[-1] - delta_x)
            pos_z.append(pos_z[-1] + delta_z)

        elif s == 'U':
            pos_x.append(pos_x[-1] - delta_z)
            pos_y.append(pos_y[-1] + delta_y)
            pos_z.append(pos_z[-1] + delta_x)

        elif s == 'D':
            pos_x.append(pos_x[-1] + delta_z)
            pos_y.append(pos_y[-1] + delta_y)
            pos_z.append(pos_z[-1] - delta_x)

    return pos_x, pos_y, pos_z


def double(list_x, list_y, list_z=False):
    """
    Checks if a protein will fold in on itself
    """

    coordinates = []

    if not list_z:
        list_z = [0] * len(list_x)

    for x, y, z in zip(list_x, list_y, list_z):
        if [x, y, z] in coordinates:
            return True
        
        coordinates.append([x, y, z])

    return False


def random_rotation_xy(list_x, list_y, n, proti):
    """
    Rotates the string 90 degrees to the left or right from the nth atom onwards.
    """

    rotation_point_x = list_x[n]
    rotation_point_y = list_y[n]

    # left or right rotation are randomly chosen
    p = random.random()

    # calculates the new positions for the remainder of the string using the
    #  equations from a 2D rotation matrix
    for i in range(n + 1, proti.length):

        relative_x = list_x[i] - rotation_point_x
        relative_y = list_y[i] - rotation_point_y

        if p > 0.5:
            # rotate left
            list_x[i] = rotation_point_x - relative_y
            list_y[i] = rotation_point_y + relative_x
        else:
            # rotate right
            list_x[i] = rotation_point_x + relative_y
            list_y[i] = rotation_point_y - relative_x

def random_rotation_xyz(list_x, list_y, list_z, n, proti):
    """
    Rotates the string 90 degrees to the left, right, up or down from the nth
    atom onwards.
    """

    rotation_point_x = list_x[n]
    rotation_point_y = list_y[n]
    rotation_point_z = list_z[n]

    # rotation direction is chosen at random
    p = random.random()

    # calculates the new positions for the remainder of the string using the
    #  equations from a 3D rotation matrix
    for i in range(n + 1, proti.length):

        relative_x = list_x[i] - rotation_point_x
        relative_y = list_y[i] - rotation_point_y
        relative_z = list_z[i] - rotation_point_z

        if p < 0.25:
            # rotate left
            list_x[i] = rotation_point_x - relative_y
            list_y[i] = rotation_point_y + relative_x
            list_z[i] = rotation_point_z + relative_z
        elif p > 0.25 and p < 0.5:
            # rotate right
            list_x[i] = rotation_point_x + relative_y
            list_y[i] = rotation_point_y - relative_x
            list_z[i] = rotation_point_z + relative_z
        elif p > 0.5 and p < 0.75:
            # rotate up
            list_x[i] = rotation_point_x - relative_z
            list_y[i] = rotation_point_y + relative_y
            list_z[i] = rotation_point_z + relative_x
        elif p > 0.75:
            # rotate down
            list_x[i] = rotation_point_x + relative_z
            list_y[i] = rotation_point_y + relative_y
            list_z[i] = rotation_point_z - relative_x

    
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



