"""
Minor programmeren

Team Proti

Case Protein Pow(d)er
"""

import matplotlib.pyplot as plt
from random import randint

DIMENSION = 2
PROTEIN = 'HHPHHHPH' # later function to automatically get this out of file
LENGTH = len(PROTEIN)

def direction():
    
    direction = 0

    while direction == 0:
        direction = randint(-DIMENSION, DIMENSION)

    return direction

def step(x, y):

    direction = direction()
    if direction in [-1, 1]:
        step = x[-1] + direction        

    else:
        step = y[-1] + direction

    return direction, step
        
def fold_protein():
    
    x_coords = [LENGTH] 
    y_coords = [LENGTH]

    first_step = True

    for i in range(len(LENGTH) - 1):
        
        step = step()
        if first_step:


    






