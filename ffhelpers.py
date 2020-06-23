"""
ffhelpers.py

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


def similarities(list_x1, list_y1, list_x2, list_y2):
    """
    Check similarities between lists
    """

    d = 0

    for x1, y1, x2, y2 in zip(list_x1, list_y1, list_x2, list_y2):
        d += math.sqrt((x2-x1)**2 +(y2-y1)**2)

    return d

def ff_plot(list_x, list_y, score, proti):
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
    fig, ax1 = plt.subplots(1, 1, figsize=(7, 7))

    ax1.plot(list_x, list_y, '--', color='darkgrey')
    ax1.plot(red_dots_x, red_dots_y, 'or', markersize=17)
    ax1.plot(blue_dots_x, blue_dots_y, 'ob', markersize=17)
    ax1.plot(yellow_dots_x, yellow_dots_y, 'oy', markersize=17)
    ax1.set_title(f'Folded protein of length {proti.length}, score: {score}')

    plt.show()
