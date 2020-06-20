"""
DEE.py
Minor Programmeren
Team Proti
Uses a Dead End Elimination like algorithm described by Okke to create a tree of all routes that are possibly lower than a preset score.
"""
from anytree import Node, RenderTree, Walker, PreOrderIter
import copy
import matplotlib.pyplot as plt
import timeit
from progress.bar import Bar
from statistics import mean
import random
from helpers import * 


def run(proti):
    length = proti.length
    protein = proti.listed

    lowest_known_score = 0

    # initiate a progress bar
    bar = Bar('Progress', max=length)
    start = timeit.default_timer()

    # direction of atom relative to the previous ones
    directions = ['L', 'R', 'S']

    # """
    # Setting this variable right allows the program to run well. 
    # It makes sure no branches are explored that couldn't possibly get a score lower than this lowest known score.
    # Setting it too high means the program is slow as it will start constructing 3^(n-1) branches.
    # Setting it too low results in finding no solutions.
    # You can find the lowest known score by running one of the other faster but less certain programs first.
    # """

    # variables to keep track of the optimal configuration
    best_score = 0
    best_x = []
    best_y = []
    loop_time = []

    keep_going = True
    # loop for constructing the tree
    for i, p in enumerate(protein):
        loop_start = timeit.default_timer()
        # stops when the remaining atoms can't add to the score
        if keep_going:

            lowest_score_k = 0
            average_scores_k = [0]
            potential_score = 0

            # keep track of how many nodes are added per depth level
            breadth_counter = 0

            # place the first atom
            # intermezzo: the loop relies heavily on the globals()[...] function
            # it creates a transforms a stirng into an actual variable and allows for dynamic naming and referencing of the nodes
            if i == 0:
                name = f'{p}{i}{i}'
                globals()[name] = Node('start')
                bar.next()
                continue

            # second atom is placed in one direction only, others would merely be a rotaion around the axis with no score advantage
            if i == 1:
                name = f'{p}{i}{i}'
                globals()[name] = Node('R', parent=globals()[f'{protein[0]}00'])
                bar.next()
                continue

            # determines how many nodes there are in in the previous tree layer
            parent_counter = len(list(PreOrderIter(globals()[f'{protein[0]}00'], filter_=lambda node: node.is_leaf)))

            # for each node, create a new one in every direction
            for j in range(len(directions) * parent_counter):

                # get the name of the parent node
                parent = f'{protein[i - 1]}{i - 1}{j}'

                # see if that parent exists, if not its score was too high and no further node along that branch has to be made
                try:
                    globals()[parent]
                except:
                    continue

                # if it exist contruct a new node for each direction
                for d in directions:

                    # name of the node to be
                    name = f'{p}{i}{breadth_counter}'

                    # see which nodes are visited already to get to this point in the tree
                    nodes_visited = [node.name for node in globals()[parent].path]
                    # the atoms still to be placed
                    nodes_to_visit = protein[i:]

                    potential_nodes_visited = copy.deepcopy(nodes_visited)
                    potential_nodes_visited.append(d)

                    # score of the protien so far
                    partial_score = partial_score_func(nodes_visited, proti)

                    # lowest score possible given the remaining atoms
                    possible_score = possible_score_func_dee(nodes_to_visit, partial_score, lowest_known_score, proti)

                    if partial_score < lowest_known_score:
                        lowest_known_score = copy.deepcopy(partial_score)

                    # new node is only made if it is possible to get a score lower than the lowest known score
                    if partial_score + possible_score >= lowest_known_score and possible_score != 0:
                        breadth_counter += 1
                        continue

                    potential_score = partial_score_func(potential_nodes_visited, proti)

                        # create new node
                    globals()[name] = Node(d, parent=globals()[parent])
                    breadth_counter += 1

                    average_scores_k.append(potential_score)
                    if potential_score < lowest_score_k:
                        lowest_score_k = copy.deepcopy(potential_score)

                    # this part only runs when the last atom is placed or the remaining atoms can't add to the score
                    if i == length - 1 or possible_score == 0:
                        # determine the path of the string
                        nodes_visited = [node.name for node in globals()[name].path]
                        pos_x, pos_y = nodes_to_xy(nodes_visited)

                        # disregard foldings onto itself
                        if double_xy(pos_x, pos_y):
                            continue

                        # get the structures score
                        current_score = score(pos_x, pos_y, proti)

                        # save if it is an improvemnt
                        if current_score <= best_score:
                            best_x = copy.deepcopy(pos_x)
                            best_y = copy.deepcopy(pos_y)
                            best_score = copy.deepcopy(current_score)

                        # stop running if the remaining atoms can't add to the score
                        if possible_score == 0:
                            keep_going = False

        # log how long each iteration takes
        loop_stop = timeit.default_timer()
        loop_time.append(loop_stop - loop_start)

        # update the progress bar
        bar.next()

    # stop the progres bar and timer
    bar.finish()
    stop = timeit.default_timer()
    runtime = stop - start
    # print some interesting information and plot result

    print(f'Length: {proti.length}')
    print(f'Score: {best_score}')
    print(f'Time: {runtime}')
    print(f'Conformation: \nx:{best_x}\ny:{best_y}')
    if len(nodes_to_visit) > 1:
        print(f'{nodes_to_visit[1:]} are not placed since they would not add to the score')
    if len(best_x) == 0:
        print('No stable solution')
    else:
        dee_plot(best_x, best_y, best_score, loop_time, runtime, proti)




if __name__ == "__main__":
    main()