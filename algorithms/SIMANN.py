"""
SIMANN.py

Minor Programmeren 
Team Proti 

Tries to find the most stable configuration of the protein in 2D
with a simulated annealing algorithm in 2D
"""

from helpers import *
import random
import copy
import timeit
from progress.bar import Bar

def run(proti):
     # start timing the run of the code
    start = timeit.default_timer()

    # initialize the protein in a straight configuration
    x = [i for i in range(proti.length)]
    y = [0] * proti.length
    
    # high number of iterations for optimising result
    iterations = 1000
    rotations = 0

    # initialize progress bar
    bar = Bar('Progress', max = iterations / 1000)

    # list to keep track of best configuration and scores
    lowest_score = 0
    best_x = []
    best_y = []
    scores = []

    # set start temperature for annealing
    start_temp = 10

    while rotations < iterations:

        # remember the previous configuration
        backup_x = copy.deepcopy(x)
        backup_y = copy.deepcopy(y)

        # remember the previous score
        old_score = score(backup_x, backup_y, proti)
        scores.append(old_score)

        # fold protein at a random amino
        rotating_amino = random.randint(0, proti.length - 1)
        random_rotation(x, y, rotating_amino, proti)

        # if protein folded into itself restore and go back
        if double_xy(x, y):
            x = backup_x
            y = backup_y
            continue

        # get the score of the current configuration
        new_score = score(x, y, proti)

        # drop temperature gradually
        temperature = start_temp * (0.997 ** rotations)
        acceptance_chance = 2 ** (-(new_score - old_score) / temperature)
        treshhold = random.random()

        # if new score is worse and it can't be accepted restore backup
        if new_score > old_score and acceptance_chance > treshhold:
            x = backup_x
            y = backup_y
            new_score = old_score

        # check if a lower score is found and remember configuration
        if new_score < lowest_score:
            best_x = copy.deepcopy(x)
            best_y = copy.deepcopy(y)
            lowest_score = copy.deepcopy(new_score)

        rotations += 1

        # continue the progress bar every thousand configurations
        if rotations % 1000 == 0:
            bar.next()

    # finish the progress bar and th timer and show information
    bar.finish()
    stop = timeit.default_timer()
    print('Runtime', stop - start, 'seconds')

    # render the output and plot the figure
    output(best_x, best_y, lowest_score, proti)
    plot_m(best_x, best_y, lowest_score, scores, proti)

if __name__ == "__main__":
    main()
