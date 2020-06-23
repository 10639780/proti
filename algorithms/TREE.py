"""
TREE.py
Minor Programmeren
Team Proti
Determines the highest possible score a protein fold can get.
"""
from anytree import Node, RenderTree, NodeMixin, LevelOrderGroupIter, PreOrderIter
import copy
import matplotlib.pyplot as plt 
from helpers import *
import timeit


def run(proti):

    start = timeit.default_timer()

    # create a tree with all possible folds
    root = create_tree(proti)

    # extract all possible protein configurations from the tree
    route_directions = create_routes(root)

    # variables to hold the optimal solution
    best_x = []
    best_y = []
    lowest_score = 0
    same_score_counter = 0
    invalid_solutions = 0

    # loop over each of the possible configurations
    for r in route_directions:

        pos_x = []
        pos_y = []
        x = 0
        y = 0

        # lists all filled with coordinates of the atoms
        for d in r:

            x += d[0]
            y += d[1]

            pos_x.append(x)
            pos_y.append(y)

        # discard all the structures that fold into themselves
        if double_xy(pos_x, pos_y):
            invalid_solutions += 1
            continue
        
        # calculate te score of the current structure
        current_score = score(pos_x, pos_y, proti)

        # save the structure with the lowest score
        if current_score < lowest_score:
            best_x = copy.deepcopy(pos_x)
            best_y = copy.deepcopy(pos_y)
            lowest_score = copy.deepcopy(current_score)
            same_score_counter = 0
        if current_score == lowest_score:
            same_score_counter += 1
    
    # print(f'Invalid solutions: {invalid_solutions}')
    # print(f'Same score solutions: {same_score_counter}')
    # print(f'Lowest score: {lowest_score}')


    # print(f'Length')
    # # plot the structure with the best score
    # tree_plot(best_x, best_y, lowest_score, proti)

    stop = timeit.default_timer()
    total_time = stop - start
    print(f'Length: {proti.length}')
    print(f'Score: {lowest_score}')
    print(f'Total runtime: {total_time}')
    print(f'Configuration:\nx:{best_x}\ny:{best_y}')
    tree_plot(best_x, best_y, lowest_score, proti)

    return total_time, lowest_score, [best_x, best_y]