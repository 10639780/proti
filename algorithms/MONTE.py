"""
MONTE.py

Minor Programmeren
Team Proti

Folds protein into the (probably) most stable state using a Monte Carlo algorithm. 
Inspired by Ramji T. Venkatasubramanian's Computational Nanomechanics assignment.
"""

# import modules
from generalhelpers import score_it, random_rotation_xy, double, output, plot
import random
import math
import copy
import timeit
from progress.bar import Bar


def run(proti):

    # start timing code
    start = timeit.default_timer()

    # create the initial straight string
    pos_x = [i for i in range(proti.length)]
    pos_y = [0] * proti.length

    rotation_counter = 0
    N = 5000

    # initialize progress bar
    bar = Bar('Progress', max=N/1000)

    # lists to keep track of the scores of each rotation and remember the one
    # with the best score
    lowest_score = 0
    best_x = []
    best_y = []
    scores = []

    # probability functions depens on temperature and the boltzmann constant, 
    # can be set to their actual values if you want to be physically responsible
    boltzmann = 1
    prob = []
    
    # loop that keeps folding the protein
    while rotation_counter < N:

        _max = 10
        scale = N/20
        center = 9 * N / 40
        temperature = _max / (1 + math.exp((rotation_counter - center) / scale)) + 0.5

        # a copy is made in case the fold is invalid or unfavourable
        log_pos_x = copy.deepcopy(pos_x)
        log_pos_y = copy.deepcopy(pos_y)

        # protein is folded randomly
        rotating_amino = random.randint(0, proti.length - 1)
        random_rotation_xy(list_x=pos_x, list_y=pos_y, n=rotating_amino, proti=proti)
        
        # check whether the protein has not folded onto itself
        if double(pos_x, pos_y):
            # if it is folded wrongly, restore to the previous configuration
            pos_x = log_pos_x
            pos_y = log_pos_y
            continue

        # calculate the scores of the old and new structure
        old_score = score_it(proti=proti, list_x=log_pos_x, list_y=log_pos_y)
        new_score = score_it(proti=proti, list_x=pos_x, list_y=pos_y)
        
        # keep track of each score
        scores.append(old_score)

        # if a score beats the old one, remember that structure 
        if new_score < lowest_score:
            best_x = copy.deepcopy(pos_x)
            best_y = copy.deepcopy(pos_y)
            lowest_score = copy.deepcopy(new_score)

        # probability function to determine whether a fold will be 'accepted'
        p = math.exp(-(new_score - old_score)/(temperature * boltzmann))

        # the treshhold for acceptance varies and is randomly determined
        treshhold = random.random()
        if p < treshhold:
            pos_x = log_pos_x
            pos_y = log_pos_y

        rotation_counter += 1
        
        # print statement for time indication in long calculations
        if rotation_counter % 1000 == 0:
            bar.next()
            prob.append(p)

    bar.finish()
    stop = timeit.default_timer()
    print('Runtime:', stop - start, 'seconds')

    # the best structure is copied to a csv file and shown in a graph
    output(proti=proti, score=lowest_score, list_x=best_x, list_y=best_y)
    plot(proti, lowest_score, best_x, best_y, 'Score after rotation', 'Rotation',
        'Score', scores=scores)
