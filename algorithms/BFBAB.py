""" 
BFBAB.py

Minor Programmeren
Team Proti

Implementation of a branch and bound protein folding algorithm as described
by Mao Chen and Wen-Qi Huang in 'Branch and Bound Algorithm for the Protein 
Folding Problem in the HP Lattice Model' (2005).
"""

# import modules

from generalhelpers import direction_to_xy, double, score_con, plot
import copy
import timeit
import random
import queue
from progress.bar import Bar
from deehelpers import possible_score_func_dee


def run(proti):

    # start timing script run time
    start = timeit.default_timer()

    # pruning chance, higher values mean a faster program but the results may not be as good
    p1 = 0.99
    p2 = 0.98

    # initialize progress bar
    bar = Bar('Progress', max=proti.length)
    bar.next()
    bar.next()
    k = 0
    amino_time = []

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

    amino_start = timeit.default_timer()
    # create a breadth first tree
    while not q.empty():

        state = q.get()

        # if all aminos are placed, put the string in a list
        state_x, state_y = direction_to_xy(state[0])
        if len(state[0]) == depth and not double(state_x, state_y):
            final_configurations.append(state)
      
        if len(state[0]) < depth:
            
            for i in ['L', 'R', 'S']:

                child = copy.deepcopy(state) 

                temp_list = list(child)
                temp_list[0] += i
                child = tuple(temp_list)

                child_x, child_y = direction_to_xy(child[0])
                if double(child_x, child_y):
                    continue
                
                
                if len(child[0]) + 1 > k:
                    amino_stop = timeit.default_timer()
                    amino_time.append(amino_stop-amino_start)
                    bar.next()
                    amino_start = timeit.default_timer()


                k = len(child[0]) + 1

                # P's are always placed, rest have some conditions
                if not proti.listed[k] == 'P':
                    score = child[1] + score_con(child_x, child_y, proti)
                    
                    possible_score = possible_score_func_dee(proti.listed[k+1:], \
                                                         score, proti.min_score, proti)

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

    best_x, best_y = direction_to_xy(best_config)
    # plot the result
    stop = timeit.default_timer()
    total_time = stop - start
    bar.finish()
    print(f'Length: {proti.length}')
    print(f'Score: {lowest_score}')
    print(f'Time: {total_time}')
    print(f'Conformation: {best_config}')
    plot(proti, lowest_score, best_x, best_y, 'Time per amino placement', 'Amino', 'Time[s]', scores=amino_time)
    return total_time, lowest_score, best_config
