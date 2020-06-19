""" 
BFBAB3D.py

Minor Programmeren
Team Proti

Attempt to implement a branch and bound protein folding algorithm as described by Mao Chen and Wen-Qi Huang 
in 'Branch and Bound Algorithm for the Protein Folding Problem in the HP Lattice Model'.
Basis for breadth first structure from Bas Terwijn's lecture.
"""
import copy
import timeit
import random
import queue
from helpers import *

def run(proti):
    # start timer
    start = timeit.default_timer()

    # specifications for depth first tree building
    depth = proti.length - 2
    q = queue.Queue()
    q.put('')
    final_configurations = []
    
    # keep track of scores per substring
    lowest_score_k = {}
    all_scores_k = {}
    lowest_score = 0
    
    # (0,1) probabilities of pruning a path, lower is more exact but less fast
    p1 = 1
    p2 = 1

    # set inital values
    for i in range(proti.length + 1):
        lowest_score_k[i] = 0
        all_scores_k[i] = [0]

    # create a breadth first tree
    while not q.empty():

        state = q.get()
        # if all aminos are placed, put the string in a list

        if len(state) == depth and not xyz_double(state):
            final_configurations.append(state)
      
        if len(state) < depth:
            for i in ['L', 'R', 'S', 'U', 'D']:

                # substring
                child = copy.deepcopy(state) 

                # string after potentialy placing the next amino
                child += i 

                # discard the string folding into themselves
                if xyz_double(child):
                    continue
                
                # identify how for into the string it is
                k = len(child) + 1
                
                # P's are always placed, rest have some conditions
                if not proti.listed[k] == 'P':

                    # score if placed 
                    score = xyz_score_func(child, proti)

                    # min score to get from remaining aminos
                    possible_score = possible_score_func(proti.listed[k+1:], score, proti.min_score)

                    if score + possible_score > lowest_score:
                        continue

                    # random number between 0 and 1
                    r = random.random()

                    # avergage of all strings of the same length
                    average_score_k = sum(all_scores_k[k]) / len(all_scores_k[k])

                    # conditions for pruning
                    if score > average_score_k and r < p1:
                        continue
                    elif (average_score_k >= score and score > lowest_score_k[k]) and r < p2:
                        continue
                    
                    # add to tree
                    q.put(child) 
                    all_scores_k[k].append(score)

                    if score < lowest_score_k[k]:
                        lowest_score_k[k] = copy.deepcopy(score)
                    
                    if score < lowest_score:
                        lowest_score = copy.deepcopy(score)
                        
                else:
                    q.put(child) 

    
    lowest_score = 0

    # weed out the best configuration from the remaining strings
    for c in final_configurations:
        if xyz_score_func(c, proti) < lowest_score:
            best_config = copy.deepcopy(c)
            lowest_score = copy.deepcopy(xyz_score_func(c, proti))

    # plot the result
    stop = timeit.default_timer()
    print(f'Strings made: {len(final_configurations)}')
    print(f'Length: {proti.length}')
    print(f'Score: {lowest_score}')
    print(f'Runtime: {stop - start}')
    xyz_plot(best_config, lowest_score, proti)


if __name__ == "__main__":
    main()