""" 
BFBA.py

Minor Programmeren
Team Proti

Attempt to implement a branch and bound protein folding algorithm as described by Mao Chen and Wen-Qi Huang 
in 'Branch and Bound Algorithm for the Protein Folding Problem in the HP Lattice Model'.
Basis for breadth first structure from Bas Terwijn's lecture.
"""
from helpers import *
import copy
import timeit
import random
import queue
from progress.bar import Bar

def run(protein_list):
    start = timeit.default_timer()
    protein = protein_list
    length = len(protein)


    # set theoretical lower bound on score
    even = protein[::2]
    odd = protein[1::2]


    H_count = protein.count('H')
    C_count = protein.count('C')
    min_score = 2 * max([- even.count('H') - 5 * even.count('C'), - odd.count('H') - 5 * odd.count('C')])


    # list of timer to see why program so slow
    possible_list = []
    score_list = []
    direction_list = []
    double_list = []


    bar = Bar('Progress', max=length)
    bar.next()
    bar.next()
    k = 0

    # specifications for depth first tree building
    depth = length - 2
    q = queue.Queue()
    q.put(('',0))
    final_configurations = []

    # keep track of scores per substring
    lowest_score_k = {}
    all_scores_k = {}
    lowest_score = 0
    
    # (0,1) probabilities of pruning a path, lower is more exact but less fast
    p1 = 0.99
    p2 = 0.9

    # set inital values
    for i in range(length + 1):
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
                if not protein[k] == 'P':
                    score = child[1] + score_func(child[0], protein)
                    
                    possible_score = possible_score_func(protein[k+1:], score, min_score)

                    if score + possible_score > lowest_score:
                        continue

                    # random number between 0 and 1
                    r = random.random()
                    
                    average_score_k = sum(all_scores_k[k]) / len(all_scores_k[k])

                    # conditions for pruning
                    if score > average_score_k and r < p1:
                        continue
                    elif (average_score_k >= score and score > lowest_score_k[k]) and r < p2:
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
    # print(f'Strings made: {len(final_configurations)}')
    total_time = stop - start
    bar.finish()
    print(f'Length: {length}')
    print(f'Score: {lowest_score}')
    print(f'Total runtime: {total_time}')
    print(f'Possible score function time: {round(sum(possible_list) / total_time * 100,1)}%')
    print(f'Score function time: {round(sum(score_list) / total_time * 100,1)}%')
    print(f'Direction function time: {round(sum(direction_list) / total_time * 100,1)}%')
    print(f'Double function time: {round(sum(double_list) / total_time * 100,1)}%')
    print(f'Conformation: {best_config}')
    plot(best_config, lowest_score, total_time, possible_list, score_list, direction_list, double_list, protein)

