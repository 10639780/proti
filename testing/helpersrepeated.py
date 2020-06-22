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

def plot_coords(proti, x_list, y_list, z_list=False):
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




    
    
                