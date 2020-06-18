"""
br_babko.py

Minor Programmeren
Team Proti

Branch and bound protein folding algorithm as described by Mao Chen and Wen-Qi Huang
in 'Branch and Bound Algorithm for the Protein Folding Problem in the HP Lattice Model'.
Breadth first with command line input interface.
"""
from helpers import *
import re
import copy
import timeit
import random
import queue
import sys

def main():
    # Input for protein variable
    algos = 'BFBAB'
    print(f'algorithms available: {algos}')
    usage = input("please input algorithm name followed by an underscore and a single protein string")
    usage = re.split('_+', usage)
    algo = usage[0]
    protein = list(usage[1])
    length = len(protein)

    allow_chars = ['P','C','H']
    for i in protein:
        if i not in allow_chars:
            sys.exit("Usage parameters not met. Please try again")

    if algo not in algos:
        sys.exit("Usage parameters not met. Please try again")

    if algo == 'BFBAB':
        start = timeit.default_timer()
        # specifications for breadth first tree building
        depth = length - 2
        q = queue.Queue()
        q.put('')
        final_configurations = []

        # keep track of scores per substring
        lowest_score_k = {}
        all_scores_k = {}
        lowest_score = 0

        # set theoretical lower bound on score
        even = protein[::2]
        odd = protein[1::2]
        min_score = 2 * max([- even.count('H') - 5 * even.count('C'), - odd.count('H') - 5 * odd.count('C')])

        # (0,1) probabilities of pruning a path, lower is more exact but less fast
        p1 = 1
        p2 = 1

        # set initial values
        for i in range(length + 1):
            lowest_score_k[i] = 0
            all_scores_k[i] = [0]

        # create a breadth first tree
        while not q.empty():

            state = q.get()

            # if all aminos are placed, put the string in a list
            if len(state) == depth and not double(state):
                final_configurations.append(state)

            if len(state) < depth:
                for i in ['L', 'R', 'S']:

                    # substring
                    child = copy.deepcopy(state)

                    # string after potentially placing the next amino
                    child += i

                    # discard the string folding into themselves
                    if double(child):
                        continue

                    # identify how for into the string it is
                    k = len(child) + 1

                    # P's are always placed, rest have some conditions
                    if not protein[k] == 'P':
                        # score if placed.
                        score = score_func(child, protein)

                        #print(f'{child}')

                        # min score to get from remaining aminos
                        possible_score = possible_score_func(protein[k + 1:], score, min_score)

                        if score + possible_score > lowest_score:
                            continue

                        # random number between 0 and 1
                        r = random.random()

                        # average of all strings of the same length
                        average_score_k = sum(all_scores_k[k]) / len(all_scores_k[k])

                        # conditions for pruning
                        if score > average_score_k and r < p1:
                            continue
                        elif (average_score_k >= score > lowest_score_k[k]) and r < p2:
                            continue
                        # elif when stringlength == half of total and child is straight: continue

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
            if score_func(c, protein) < lowest_score:
                best_config = copy.deepcopy(c)
                lowest_score = copy.deepcopy(score_func(c, protein))

        # plot the result
        stop = timeit.default_timer()
        # print(f'Strings made: {len(final_configurations)}')
        print(f'Length: {length}')
        print(f'Score: {lowest_score}')
        print(f'Runtime: {stop - start}')
        plot(best_config, lowest_score, protein)




if __name__ == "__main__":
    main()