"""
treehelpers.py

Minor Programmeren 
Team Proti

Algorithms are implemented by calling the functions 
in this script.
"""

# import modules
import matplotlib.pyplot as plt
from anytree import PreOrderIter
import copy
from classes.amino import *


def create_tree(proti):
    """Create a tree structure, each atom can branch off into 4 directions, creating 4^n possible routes."""
    directions = [[1, 0], [-1, 0], [0, 1], [0, -1]]
    node_counter = 0
    parent_counter = 0
    node_list = []
    first = True
    second = True

    for p, i in zip(proti.listed, [i for i in range(proti.length)]):

        # first atom has no parents and is located at 0,0
        if first:
            root = Atom(p, [0, 0])
            node_list.append(root)
            first = False
            continue
        
        # second atom
        if second:
            node = Atom(p, [1, 0], parent=root)
            node_list.append(node)
            second = False
            continue

        for k in range((len(directions) - 1) ** (i - 2)):

            # make sure the protein doesn't walk back into itself
            previous_direction = copy.deepcopy(node_list[parent_counter + 1].direction)
            previous_direction[0] = -previous_direction[0]
            previous_direction[1] = -previous_direction[1]

            temp_direction = copy.deepcopy(directions)

            if previous_direction in temp_direction: 
                temp_direction.remove(previous_direction)

            # fill the tree and keep track of the parent child relations
            for j in range(len(directions) - 1):
            
                j = Atom(p, temp_direction[j], parent=node_list[parent_counter + 1])
                node_list.append(j)

                node_counter += 1
                if node_counter % (len(directions) - 1) == 0:
                    parent_counter += 1

    return root


def create_routes(root):
    """Creates a list that hold all the possible configurations a protein could have."""

    route_directions = []
    
    # list of all possible paths in the tree
    # source:
    # https://stackoverflow.com/questions/59917058/how-to-get-all-possible-branch-with-python-anytree
    paths = [list(leaf.path) for leaf in PreOrderIter(root, filter_=lambda no_de: no_de.is_leaf)]

    # add those paths to a list
    for p in paths:
        direction = []

        for node in p:
            direction.append(node.direction)

        route_directions.append(direction)
    
    return route_directions


def tree_plot(list_x, list_y, score, proti):
    """Makes a graph of two lists list_x, list_y."""
    
    # differentiate between types of atom
    red_dots_x = []
    red_dots_y = []
    blue_dots_x = []
    blue_dots_y = []
    yellow_dots_x = []
    yellow_dots_y = []

    # search through protein and place each atom in the appropriate list
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

    plt.plot(list_x, list_y, '--', color='darkgrey')
    plt.plot(red_dots_x, red_dots_y, 'or')
    plt.plot(blue_dots_x, blue_dots_y, 'ob')
    plt.plot(yellow_dots_x, yellow_dots_y, 'oy')
    plt.title(f'Folded protein, score: {score}')

    plt.show()
