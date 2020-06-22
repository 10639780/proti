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

def plot(proti, score, scores, list_x=False, list_y=False, list_z=False, string=False):

    if string:
        list_x, list_y, list_z = directions(string)

    if all(elem == 0 for elem in list_z):
        print("Not 3D Folded")
        list_z = False

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


if __name__ == "__main__":
    x = [i for i in range(10)]
    y = [i for i in range(10)]
    z = [i for i in range(10)]
    scores = [0] * 10
    score = 0

    proteinstring = 'HPHPPHHPHPPHPHHPPHPH'
    length = len(proteinstring)
    proti = Protein(proteinstring, length)

    string = 'LRSLRSLRS'
    test1 = plot(proti=proti, score=score, scores=scores, string=string)


