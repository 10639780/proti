"""
monteko.py

Minor Programmeren
Team Proti

Outputs plot based on simplified BT string.
"""
import numpy as np 
import matplotlib.pyplot as plt 
import random
import math
import copy
import csv
import sys
import re

protein = []
length = len(protein)

# Declare specific folding config with BT-string

BT_string = input("please input a single BT-string, output: images/BT-string.png") 

def main():
    # lists that keep track of coordinates.
    # Starting point can be altered so that a config that crosses either x=0 or y=0, is pruned. 
    
    best_x = [5]
    best_y = [10]
        
    for char in BT_string:
        if char.isalpha() == True: 
            protein.append(char)
    
    BT_list = re.findall('[A-Z][^A-Z]*', BT_string)
    BTS_list = []
    for i in BT_list:
        i = re.sub('[HPC]', '', i)
        BTS_list.append(int(i))
    
    for i in BTS_list: 
        if i > 0:
            if i == 1:
                best_x.append(best_x[-1] + 1)
                best_y.append(best_y[-1])
            if i == 2:
                best_y.append(best_y[-1] + 1)
                best_x.append(best_x[-1])
        if i < 0:
            if i == -1:
                best_x.append(best_x[-1] - 1)
                best_y.append(best_y[-1])
            if i == -2:
                best_y.append(best_y[-1] - 1)
                best_x.append(best_x[-1])
        if i == 0:
            pass
    
    
    

    # the best structure is shown in a graph
    plot(best_x, best_y)

def plot(list_x, list_y):
    """Makes a graph of two lists list_x, list_y."""
    
    # differentiate between types of elements
    red_dots_x = []
    red_dots_y = []
    blue_dots_x = []
    blue_dots_y = []
    yellow_dots_x = []
    yellow_dots_y = []

    # search through protein and place each atom in the appropriate list
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
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3, 15))
    groups = ['H2H1P-2', 'H2P-1H-1P-2P1H-2']
    colorcounter = 0
    for i in groups:
        if i in BT_string:
            colorcounter += 1
            
    if colorcounter == 1:
        ax1.plot(list_x, list_y, '--', color='orange')
    
    if colorcounter == 2:
        ax1.plot(list_x, list_y, '--', color='purple')
        
    elif colorcounter == 0:
        ax1.plot(list_x, list_y, '--', color='darkgrey')
    
    ax1.plot(red_dots_x, red_dots_y, 'or')
    ax1.plot(blue_dots_x, blue_dots_y, 'ob')
    ax1.plot(yellow_dots_x, yellow_dots_y, 'oy')
    ax1.set_title(f'Folded protein')
    fig.delaxes(ax2)
    fig.savefig('images/' + '{}.png'.format(BT_string) , bbox_inches='tight')
    #plt.show()
    print(f'see images folder')
    
if __name__ == "__main__":
    main()