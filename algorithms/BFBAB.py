""" 
BFBAB.py

Minor Programmeren
Team Proti

Attempt to implement a branch and bound protein folding algorithm as described
by Mao Chen and Wen-Qi Huang in 'Branch and Bound Algorithm for the Protein 
Folding Problem in the HP Lattice Model'.
Basis for breadth first structure from Bas Terwijn's lecture.
"""

# import modules
from helpers import *
import copy
import timeit
import random
import queue
from progress.bar import Bar

def run(proti):

    # start timing script run time
    start = timeit.default_timer()

    p1 = 0.99
    p2 = 0.98

    # initialize progress bar
    bar = Bar('Progress', max=proti.length)
    bar.next()
    bar.next()
    k = 0

    # specifications for breadth first tree building
    depth = proti.length - 2
    q = queue.Queue()
    q.put(('',0))
    final_configurations = []

    # keep track of scores per substring
    lowest_score_k = {}
    all_scores_k = {}
    lowest_score = 0

    # set inital values
    for i in range(proti.length + 1):
        lowest_score_k[i] = 0
        all_scores_k[i] = [0]

    # create a breadth first tree
    while not q.empty():

        state = q.get()

        # if all aminos are placed, put the string in a list
        if len(state[0]) == depth and not double(state[0]):
            final_configurations.append(state)
      
        if len(state[0]) < depth:
            
            for i in ['L', 'R', 'S']:

                child = copy.deepcopy(state) 

                temp_list = list(child)
                temp_list[0] += i
                child = tuple(temp_list)

                if double(child[0]):
                    continue
                
                if len(child[0]) + 1 > k:
                    bar.next()
          
                k = len(child[0]) + 1

                # P's are always placed, rest have some conditions
                if not proti.listed[k] == 'P':
                    score = child[1] + score_func(child[0], proti)
                    
                    possible_score = possible_score_func(proti.listed[k+1:], \
                                                         score, proti.min_score)

                    if score + possible_score > lowest_score:
                        continue

                    # random number between 0 and 1
                    r = random.random()
                    
                    average_score_k = sum(all_scores_k[k]) / len(all_scores_k[k])

                    # conditions for pruning
                    if score > average_score_k and r < p1:
                        continue
                    elif (average_score_k >= score and score > lowest_score_k[k])\
                         and r < p2:
                        continue
                    
                    temp_list = list(child)
                    temp_list[1] = score
                    child = tuple(temp_list)  

                    # add to tree
                    q.put(child) 

                    all_scores_k[k].append(score)

                    if score < lowest_score_k[k]:
                        lowest_score_k[k] = copy.deepcopy(score)
                    
                    if score < lowest_score:
                        lowest_score = copy.deepcopy(score)
                        
                else:
                    q.put(child) 

    if len(final_configurations) == 0:
        bar.finish()
        print('No conformations found.')
        return

    lowest_score = 0
    
    # weed out the best configuration from the remaining strings
    for c in final_configurations:
        if c[1] < lowest_score:
            best_config = copy.deepcopy(c[0])
            lowest_score = copy.deepcopy(c[1])

    # plot the result
    stop = timeit.default_timer()
    total_time = stop - start
    bar.finish()
    print(f'Length: {proti.length}')
    print(f'Score: {lowest_score}')
    print(f'Total runtime: {total_time}')
    print(f'Configuration: {best_config}')
    # plot(best_config, lowest_score, total_time, proti)

    return total_time, lowest_score, best_config

