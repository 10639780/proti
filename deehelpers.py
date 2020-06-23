"""
deehelpers.py

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
from generalhelpers import double, score_it


def possible_score_func_dee(nodes_to_visit, partial_score, lowest_known_score, proti):
    """
    Calculates the best score the remaining bit of the protien can acquire.
    """

    possible_score = 0

    # an H atom can get at most -2 and a C atom at best -10
    for i, n in enumerate(nodes_to_visit):

        if i == len(nodes_to_visit) - 1:
            if n == 'H':
                possible_score += -3
            if n == 'C':
                possible_score += -15
            continue

        if n == 'H':
            possible_score += -2
        if n == 'C':
            possible_score += -10

    if possible_score + partial_score < proti.min_score:
        possible_score = proti.min_score - partial_score
        if lowest_known_score == proti.min_score:
            possible_score += -1

    return possible_score

def partial_score_func(nodes_visited, proti):
    """
    Calculates the socre obtained by the nodes visit so far.
    """

    pos_x, pos_y = nodes_to_xy(nodes_visited)

    # weed out the proteins folded into itself
    if double(list_x=pos_x, list_y=pos_y):
        return 1000

    partial_score = score_it(proti=proti, list_x=pos_x, list_y=pos_y)

    return partial_score

def nodes_to_xy(nodes_visited):
    """
    Converts a series of string with directions like ['left', 'right'] 
    to lists with xy positions.
    """

    pos_x = [0, 1]
    pos_y = [0, 0]

    # go over every node
    for i, n in enumerate(nodes_visited):

        # first two nodes are already placed to confine solutions to one quadrant
        if i == 0 or i == 1:
            continue

        # previous direction is determined
        delta_x = pos_x[-1] - pos_x[-2]
        delta_y = pos_y[-1] - pos_y[-2]

        # rotation matrices used to turn into the desired direction
        if n == 'S':
            pos_x.append(pos_x[-1] + delta_x)
            pos_y.append(pos_y[-1] + delta_y)
        elif n == 'L':
            pos_x.append(pos_x[-1] - delta_y)
            pos_y.append(pos_y[-1] + delta_x)
        elif n == 'R':
            pos_x.append(pos_x[-1] + delta_y)
            pos_y.append(pos_y[-1] - delta_x)

    return pos_x, pos_y

def dee_plot(list_x, list_y, score, loop_time, total_time, proti):
    """
    Makes a graph of two lists list_x, list_y.
    """

    # differentiate between types of atom
    red_dots_x = []
    red_dots_y = []
    blue_dots_x = []
    blue_dots_y = []
    yellow_dots_x = []
    yellow_dots_y = []

    # search through protein and place each atom in the appropiate list
    for x, y, p in zip(list_x, list_y, proti.listed):

        if p == 'H':
            red_dots_x.append(x)
            red_dots_y.append(y)
        if p == 'P':
            blue_dots_x.append(x)
            blue_dots_y.append(y)
        if p == 'C':
            yellow_dots_x.append(x)
            yellow_dots_y.append(y)

    # create graphs with colors
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9))

    ax1.plot(list_x, list_y, '--', color='darkgrey')
    ax1.plot(red_dots_x, red_dots_y, 'or', markersize=15)
    ax1.plot(blue_dots_x, blue_dots_y, 'ob', markersize=15)
    ax1.plot(yellow_dots_x, yellow_dots_y, 'oy', markersize=15)
    ax1.axis('equal')
    ax1.set_title(f'Folded protein of length {proti.length}, score: {score}')

    ax2.plot(loop_time)
    ax2.set_title(f'Time per atom, {round(total_time, 2)} seconds total')
    ax2.set(xlabel='Atom', ylabel='Time')

    plt.show()
