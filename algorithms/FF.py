"""
FF.py

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
from helpers import *

def run(proti):
    length = proti.length
    strings_created = 0
    lowest_score = 0
    population = 100
    N = 20

    routes = []    
    best_route = []

    start = timeit.default_timer()
    bar = Bar('Progress', max=N/100)

    # create a swarm of protein strings that dont fold into themselves
    while strings_created < population:
        route = []

        for i in range(length - 2):
            route.append(['S', 'L', 'R'][random.randint(0,2)])
        
        if not double(route):
            routes.append(route)
            strings_created += 1
       
    for k in range(N):
        # determine which route has the lowest score 
        for r in routes:

            # pos_x, pos_y = direction_to_xy(r)

            while double(r):
                route = []

                for i in range(length - 2):
                    route.append(['S', 'L', 'R'][random.randint(0,2)])
                
                # pos_x, pos_y = direction_to_xy(route)
                if not double(pos_x, pos_y):
                    r = route

            score = score_whole_string(r, proti)

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
                    r[random.randint(0, len(r) - 1)] = ['S', 'L', 'R'][random.randint(0,2)]
                    r_x, r_y = direction_to_xy(r)
                    if not double_xy(r_x, r_y):
                        invalid = False

            r_x, r_y = direction_to_xy(r)
            after = similarities(r_x, r_y, b_r_x, b_r_y)
        
            if after < before:
                r = log_r

        if k % 100 == 0:
            bar.next()

    pos_x, pos_y = direction_to_xy(best_route)
    score = score_whole_string(best_route, proti)
    bar.finish()
    stop = timeit.default_timer()
    runtime = stop - start
    route_string = ''.join(best_route)
    print(f'Length: {proti.length}')
    print(f'Score: {score}')
    print(f'Time: {runtime}')
    print(f'Conformation: {route_string}')
    ff_plot(pos_x, pos_y, score, proti)
