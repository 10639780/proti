"""
monte.py

Minor Programmeren
Team Proti

Folds protein into the (probably) most stable state using a Monte Carlo algorithm. 
"""
import numpy as np 
import matplotlib.pyplot as plt 
import random
import math
import copy
from mpl_toolkits import mplot3d

# string and length are used by most functions so declare as global variable
# protein = ['H','H','P','H','H','H','P','H']
<<<<<<< HEAD
protein = ['H', 'C', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'C', 'H', 'C', 'H', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'H', 'H', 'H', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'H']
# protein = ['H', 'H', 'H', 'P', 'C', 'C', 'H', 'P', 'C', 'C', 'P', 'H']
=======
# protein = ['H', 'C', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'C', 'H', 'C', 'H', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'H', 'H', 'H', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'H']
protein = ['C', 'P', 'P', 'C', 'H', 'P', 'P', 'C', 'H', 'P', 'P', 'C', 'P', 'P', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'P', 'C', 'H', 'P', 'P', 'C', 'P', 'C', 'H', 'P', 'P', 'H', 'P', 'C']
>>>>>>> 9698fda472b99522490cb347b94db3009cf5b660
length = len(protein)

def main():

    # create the initial straight string
    pos_x = [i for i in range(length)]
    pos_y = [0] * length
    pos_z = [0] * length

    # number of iterations, higher N is better result
    N = 100000
    rotation_counter = 0

    # lists to keep track of the scores of each rotation and remember the one with the best score
    lowest_score = 0
    best_x = []
    best_y = []
    best_z = []
    scores = []

    # probability functions depens on temperature and the boltzmann constant, can be set to their actual values if you want to be physically responsible
    temperature = 1
    boltzmann = 1
    
    # loop that keeps folding the protein
    while rotation_counter < N:

        # a copy is made in case the fold is invalid or unfavourable
        log_pos_x = copy.deepcopy(pos_x)
        log_pos_y = copy.deepcopy(pos_y)
        log_pos_z = copy.deepcopy(pos_z)

        # protein is folded randomly
        random_rotation(pos_x, pos_y, pos_z, random.randint(0, length - 1))

        # check whether the protein has not folded onto itself
        if double(pos_x, pos_y, pos_z):
            # if it is folded wrongly, restore to the previous configuration
            pos_x = log_pos_x
            pos_y = log_pos_y
            pos_z = log_pos_z
            continue
        
        # calculate the scores of the old and new structure
        old_score = score(log_pos_x, log_pos_y, log_pos_z)
        new_score = score(pos_x, pos_y, pos_z)

        # keep track of each score
        scores.append(old_score)

        # if a score beats the old one, remember that structure 
        if new_score < lowest_score:
            best_x = copy.deepcopy(pos_x)
            best_y = copy.deepcopy(pos_y)
            best_z = copy.deepcopy(pos_z)
            lowest_score = copy.deepcopy(new_score)

        # probability function to determine whether a fold will be 'accepted' or not, a lower score relative to the previous configuration increases the changes of adoption
        p = math.exp(-(new_score - old_score)/(temperature * boltzmann))

        # the treshhold for acceptance varies and is randomly determined
        treshhold = random.random()
        if p < treshhold:
            pos_x = log_pos_x
            pos_y = log_pos_y
            pos_z = log_pos_z

        rotation_counter += 1

        # print statement for time indication in long calculations
        if rotation_counter % 1000 == 0:
            print(f'{rotation_counter / N * 100}%')

    # the best structure is copied to a csv file and shown in a graph
    output(best_x, best_y, best_z, lowest_score)
    plot(best_x, best_y, best_z, lowest_score, scores)

def output(list_x, list_y, list_z, score):
    """Prints the folded string to a csv file in the Bas Terwijn style."""

    numbers = []
   
    for i in range(len(protein)-1):
        
        # new position is compared to the old
        delta_x = list_x[i+1] - list_x[i]
        delta_y = list_y[i+1] - list_y[i]
        delta_z = list_z[i+1] - list_z[i]

        # conversion between coordinates  and bas terwijn numbers
        if delta_x == 1:
            number = -2
        elif delta_x == -1:
            number = 2
        elif delta_y == 1:
            number = 1
        elif delta_y == -1:
            number = -1
        elif delta_z == 1:
            number = 3
        elif delta_z == -1:
            number = -3
        numbers.append(number)

    # add 0 to signal end of protein
    numbers.append(0)

    # write the list to a file
    f = open('output.csv', 'w')
    f.write('amino,fold\n')
    for p, n in zip(protein, numbers):
        f.write(f'{p}, {n}\n')
    f.write(f'score,{score}') 
    f.close()
  

def score(list_x, list_y, list_z):
    """Given the coordinates of a protein string, calculate the score of the shape."""

    # list to place the 'already scored' atoms into
    coordinates = []
    directions = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]
    score = 0

    for i in range(length):

        # P's dont interact so can skip those cases
        if not protein[i] == 'P':

            # for every atom look around in all 6 directions
            for d in directions:

                # check whether one of the previously placed atoms is in the vicinity and determine the score of the interaction with it and the current atom
                for j in range(len(coordinates)):

                    if [list_x[i] + d[0], list_y[i] + d[1], list_z[i] + d[2]] == coordinates[j] and not(list_x[i] + d[0] == list_x[i-1] and list_y[i] + d[1] == list_y[i-1] and list_z[i] + d[2] == list_z[i-1]):
                        
                        if protein[i] == 'H':
                            if protein[j] == 'H' or protein[j] == 'C':
                                score += -1

                        if protein[i] == 'C':
                            if protein[j] == 'C':
                                score += -5
                            if protein[j] == 'H':
                                score += -1 

        # place in the list with coordinates
        coordinates.append([list_x[i], list_y[i], list_z[i]])
    
    return score


def double(list_x, list_y, list_z):
    """Checks whether two atoms occupy the same point."""

    coordinates = []

    # see if a coordinate is already in the list, then add that coordinate to the list
    for x, y, z in zip(list_x, list_y, list_z):
        if [x,y,z] in coordinates:
            return True
        coordinates.append([x,y,z])
    
    return False


def plot(list_x, list_y, list_z, score, scores):
    """Makes a graph of two lists list_x, list_y."""
    
    # differentiate between types of atom
    red_dots_x = []
    red_dots_y = []
    red_dots_z = []
    blue_dots_x = []
    blue_dots_y = []
    blue_dots_z = []
    yellow_dots_x = []
    yellow_dots_y = []
    yellow_dots_z = []

    # search through protein and place each atom in the appropiate list
    for x, y, z, p in zip(list_x, list_y, list_z, protein):

        if p == 'H':
            red_dots_x.append(x)
            red_dots_y.append(y)
            red_dots_z.append(z)
        if p == 'P':
            blue_dots_x.append(x)
            blue_dots_y.append(y) 
            blue_dots_z.append(z)       
        if p == 'C':
            yellow_dots_x.append(x)
            yellow_dots_y.append(y)
            yellow_dots_z.append(z)


    # colored graph of protein
    fig = plt.figure()
    ax = plt.axes(projection='3d')
  
    ax.plot(list_x, list_y, list_z, '--', color='darkgrey')
    ax.plot(red_dots_x, red_dots_y, red_dots_z, 'or')
    ax.plot(blue_dots_x, blue_dots_y, blue_dots_z, 'ob')
    ax.plot(yellow_dots_x, yellow_dots_y, yellow_dots_z, 'oy')
    ax.set_title(f'Folded protein, score: {score}')
    plt.savefig("monte3dfig.png")
    plt.show()

    # graph of scores
    plt.plot(scores)
    plt.title('Scores of the configurations after each rotatation')
    plt.xlabel('Rotation')
    plt.ylabel('Score')
    plt.savefig("monte3dstats.png")
    plt.show()


def random_rotation(list_x, list_y, list_z, n):
    """Rotates the string 90 degrees to the left, right, up or down from the nth atom onwards."""

    rotation_point_x = list_x[n]
    rotation_point_y = list_y[n]
    rotation_point_z = list_z[n]

    # rotation direction is chosen at random
    p = random.random()

    # calculates the new positions for the remainder of the string using the equations from a 3D rotation matrix
    for i in range(n + 1, length):
       
        relative_x = list_x[i] - rotation_point_x
        relative_y = list_y[i] - rotation_point_y
        relative_z = list_z[i] - rotation_point_z

        if p < 0.25:
            # rotate left
            list_x[i] = rotation_point_x - relative_y
            list_y[i] = rotation_point_y + relative_x
            list_z[i] = rotation_point_z + relative_z
        elif p > 0.25 and p < 0.5:
            # rotate right
            list_x[i] = rotation_point_x + relative_y
            list_y[i] = rotation_point_y - relative_x
            list_z[i] = rotation_point_z + relative_z
        elif p > 0.5 and p < 0.75:
            # rotate up
            list_x[i] = rotation_point_x - relative_z
            list_y[i] = rotation_point_y + relative_y
            list_z[i] = rotation_point_z + relative_x        
        elif p > 0.75:
            #rotate down
            list_x[i] = rotation_point_x + relative_z
            list_y[i] = rotation_point_y + relative_y
            list_z[i] = rotation_point_z - relative_x   

if __name__ == "__main__":
    main()