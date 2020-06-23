""" 
genetic.py

Minor Programmeren
Team Proti

Genetic protein folding algorithm based on 'An Efficient Genetic Algorithm for Predicting
Protein Tertiary Structures in the 2D HP Model' by Thang N. Bui and Gnanasekaran Sundarraj 
and on 'Genetic Algorithm for Predicting Protein Folding in the 2D HP Model' by Eyal Halm.
"""

import random
import operator
import copy
import matplotlib.pyplot as plt
from progress.bar import Bar
import timeit
import math
from helpers import *

def run(proti):
    # some constants
    population_size = 300
    N = 6000
    sample_size = round(population_size * 0.1)
    _max = 0.4
    _min = 0.2
    scale = N/20
    center = 9 * N / 40
    bar = Bar('Progress', max=N/1000)
    start = timeit.default_timer()
    best_yet = []

    # initialize by creating a population of random conformations
    conformations_list = initial_population(population_size, proti)
    initial_time = timeit.default_timer()

    for i in range(N): 
        # determine the chance that a mutation occurs
        mutation_chance = _max / (1 + math.exp((i - center )/ scale)) + _min

        # choose 2 parents from the population, make 2 children
        parent1, parent2, child1, child2  = make_child(conformations_list, sample_size)   

        # mutate the children 
        mutation1, mutation2 = mutate(child1, child2, mutation_chance, proti)

        # replace the weakest parent/member of population with the best child
        replace(parent1, parent2, mutation1, mutation2, conformations_list)

        # keep track of lowest score
        best_yet.append(min([parent1, parent2, mutation1, mutation2], key=operator.itemgetter(1))[1])

        # update the progress bar
        if i % 1000 == 0:
            bar.next()

    # output
    best_conformation =  min(conformations_list, key=operator.itemgetter(1))
    bar.finish()
    stop = timeit.default_timer()
    initialization_time = initial_time - start
    runtime = stop - initial_time
    total_time = stop - start
    print(f'Length: {proti.length}')
    print(f'Score: {best_conformation[1]}')
    print(f'Initialization: {initialization_time}')
    print(f'Runtime: {runtime}')
    print(f'Total: {total_time}')
    print(f'Conformation: {best_conformation[0]}')
    genetic_plot(best_conformation[0], best_conformation[1], best_yet, proti)

    return total_time, best_conformation[1], best_conformation[0]