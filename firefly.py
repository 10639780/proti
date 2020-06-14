"""
firefly.py

Minor Programmeren
Team Proti

Attempt to implement the firefly algorithm as described by Yudong Zhang, LenanWu, and Shuihua Wang
in 'Solving Two-Dimensional HP Model by Firefly Algorithm and Simplified Energy Function',
and Neal Lesh, Michael Mitzenmacher and Sue Whitesides 
in 'A Complete and Effective Move Set for Simplified Protein Folding'.
"""
import pandas as pd 
from scipy.spatial import distance_matrix
from scipy.spatial.distance import squareform
import timeit
import numpy as np 
import math
import random
import copy
import matplotlib.pyplot as plt 
from progress.bar import Bar

protein = ['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H'] # official 8, mc -3
# protein = ['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'P', 'H'] # official 14, mc -6
protein = ['H', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'H', 'P', 'P', 'H', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'H'] # official 20, mc -9
protein = ['P', 'P', 'H', 'P', 'P', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H'] # opt -21
protein = ['H', 'P', 'P', 'H', 'P', 'P', 'H', 'P', 'H', 'P', 'P', 'H', 'P']

length = len(protein)


def main():

    N = 100
    iterations = 1000
    strings_created = 0
    routes = []
    lowest_score = 0
    best_route = []
    bar = Bar('Progress', max=iterations)

    # create a swarm of N protein strings that dont fold into themselves
    while strings_created < N:
        route = []

        for i in range(length - 2):
            route.append(random.randint(0,2))
        
        pos_x, pos_y = direction_to_xy(route)
        if not double(pos_x, pos_y):
            routes.append(route)
            strings_created += 1
       
    for k in range(iterations):
        # determine which route has the lowest score 
        for r in routes:

            pos_x, pos_y = direction_to_xy(r)

            while double(pos_x, pos_y):
                route = []

                for i in range(length - 2):
                    route.append(random.randint(0,2))
                
                pos_x, pos_y = direction_to_xy(route)
                if not double(pos_x, pos_y):
                    r = route

            score = score_func(pos_x, pos_y)

            if score <= lowest_score:
                lowest_score = copy.deepcopy(score)
                best_route = copy.deepcopy(r)

        # bend the rest so as to be more like the best score
        for r in routes:

            log_r = copy.deepcopy(r)
            r_x, r_y = direction_to_xy(r)
            b_r_x, b_r_y = direction_to_xy(best_route)
            before = similarities(r_x, r_y, b_r_x, b_r_y)
            invalid = True
            for i in range(10):
                while invalid:
                    r[random.randint(0, len(r) - 1)] = random.randint(0, 2)
                    r_x, r_y = direction_to_xy(r)
                    if not double(r_x, r_y):
                        invalid = False

            r_x, r_y = direction_to_xy(r)
            after = similarities(r_x, r_y, b_r_x, b_r_y)
        
            if after < before:
                r = log_r
        
        bar.next()
        # print(lowest_score)

     
    bar.finish()
    stop = timeit.default_timer()
    print(f'Runtime: {round(stop - start, 20)} seconds')



    pos_x, pos_y = direction_to_xy(best_route)
    score = score_func(pos_x, pos_y)
    plot(pos_x, pos_y, score)



def score_func(list_x, list_y):
    """Given the coordinates of a protein string, calculate the score of the shape."""

    # list to place the 'already scored' atoms into
    coordinates = []
    directions = [[-1,0],[0,1],[1,0],[0,-1]]
    length = len(list_x)
    score = 0

    for i in range(length):

        # P's dont interact so can skip those cases
        if not protein[i] == 'P':

            # for every atom look around in all 4 directions
            for d in directions:

                # check whether one of the previously placed atoms is in the vicinity and determine the score of the interaction with it and the current atom
                for j in range(len(coordinates)):

                    if [list_x[i] + d[0], list_y[i] + d[1]] == coordinates[j] and not(list_x[i] + d[0] == list_x[i-1] and list_y[i] + d[1] == list_y[i-1]):
                        
                        if protein[i] == 'H':
                            if protein[j] == 'H' or protein[j] == 'C':
                                score += -1

                        if protein[i] == 'C':
                            if protein[j] == 'C':
                                score += -5
                            if protein[j] == 'H':
                                score += -1 

        # place in the list with coordinates
        coordinates.append([list_x[i], list_y[i]])
    
    return score


def plot(list_x, list_y, score):
    """Makes a graph of two lists list_x, list_y."""
    
    # differentiate between types of atom
    red_dots_x = []
    red_dots_y = []
    blue_dots_x = []
    blue_dots_y = []
    yellow_dots_x = []
    yellow_dots_y = []

    # search through protein and place each atom in the appropiate list
    for x, y, p in zip(list_x, list_y, protein):

        if p == 'H':
            red_dots_x.append(x)
            red_dots_y.append(y)
        if p == 'P':
            blue_dots_x.append(x)
            blue_dots_y.append(y)       
        if p == 'C':
            yellow_dots_x.append(x)
            yellow_dots_y.append(y)

   # create graphs with colors
    fig, ax1 = plt.subplots(1, 1, figsize=(7, 7))

    ax1.plot(list_x, list_y, '--', color='darkgrey')
    ax1.plot(red_dots_x, red_dots_y, 'or', markersize=17)
    ax1.plot(blue_dots_x, blue_dots_y, 'ob', markersize=17)
    ax1.plot(yellow_dots_x, yellow_dots_y, 'oy', markersize=17)
    ax1.set_title(f'Folded protein of length {length}, score: {score}')

    plt.show()



def direction_to_xy(nodes_visited):
    """Converts a series of string with directions like ['left', 'right'] to lists with xy positions."""
   
    pos_x = [0,1]
    pos_y = [0,0]

    # go over every node
    for i, n in enumerate(nodes_visited):
        
        # previous direction is determined
        delta_x = pos_x[-1] - pos_x[-2]
        delta_y = pos_y[-1] - pos_y[-2]

        # rotation matrices used to turn into the desired direction
        if n == 0: #straight
            pos_x.append(pos_x[-1] + delta_x)
            pos_y.append(pos_y[-1] + delta_y)
        elif n == 1: # left
            pos_x.append(pos_x[-1] - delta_y)
            pos_y.append(pos_y[-1] + delta_x)
        elif n == 2: # right
            pos_x.append(pos_x[-1] + delta_y )
            pos_y.append(pos_y[-1] - delta_x)

    return pos_x, pos_y

  
def double(list_x, list_y):
    """Checks whether two atoms occupy the same point."""

    coordinates = []

    # see if a coordinate is already in the list, then add that coordinate to the list
    for x, y in zip(list_x, list_y):
        if [x,y] in coordinates:
            return True
        coordinates.append([x,y])
    
    return False


def similarities(list_x1, list_y1, list_x2, list_y2):
    d = 0

    for x1, y1, x2, y2 in zip(list_x1, list_y1, list_x2, list_y2):
        d += math.sqrt((x2-x1)**2 +(y2-y1)**2)

    return d

def simplified_score(list_x, list_y):

    coords = []

    for x, y, p in zip(list_x, list_y, protein):
        if p != 'P':
            coords.append([x,y])

    df = pd.DataFrame(coords)
    dm = pd.DataFrame(distance_matrix(df.values, df.values))
    upper_sum = squareform(dm, checks=False).sum()

    return -upper_sum


def random_rotation(list_x, list_y):
    """Rotates the string 90 degrees to the left or right from the nth atom onwards."""

    n = random.randint(0, length - 1)

    rotation_point_x = list_x[n]
    rotation_point_y = list_y[n]

    # left or right rotation are randomly chosen
    p = random.random()

    # calculates the new positions for the remainder of the string using the equations from a 2D rotation matrix
    for i in range(n + 1, length):
       
        relative_x = list_x[i] - rotation_point_x
        relative_y = list_y[i] - rotation_point_y

        if p > 0.5:
            # rotate left
            list_x[i] = rotation_point_x - relative_y
            list_y[i] = rotation_point_y + relative_x
        else:
            # rotate right
            list_x[i] = rotation_point_x + relative_y
            list_y[i] = rotation_point_y - relative_x

if __name__ == "__main__":
    start = timeit.default_timer()
    main()
