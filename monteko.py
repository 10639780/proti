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
import csv

# string and length are used by most functions so declare as global variable
# protein = ['H','H','P','H','H','H','P','H']
# protein = ['H', 'C', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'C', 'H', 'C', 'H', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'H', 'H', 'H', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'H']
protein = ['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'P', 'H']
length = len(protein)

def main():
    for i in range(500):
        # create the initial straight string
        pos_x = [i for i in range(length)]
        pos_y = [0] * length

        # number of iterations, higher N is better result
        N = 10000
        rotation_counter = 0

        # lists to keep track of the scores of each rotation and remember the one with the best score
        lowest_score = 0
        best_x = []
        best_y = []
        scores = []

        # probability functions depens on temperature and the boltzmann constant, can be set to their actual values if you want to be physically responsible
        temperature = 1
        boltzmann = 1
    
        # loop that keeps folding the protein
        while rotation_counter < N:

            # a copy is made in case the fold is invalid or unfavourable
            log_pos_x = copy.deepcopy(pos_x)
            log_pos_y = copy.deepcopy(pos_y)

            # protein is folded randomly
            random_rotation(pos_x, pos_y, random.randint(0, length - 1))

            # check whether the protein has not folded onto itself
            if double(pos_x, pos_y):
                # if it is folded wrongly, restore to the previous configuration
                pos_x = log_pos_x
                pos_y = log_pos_y
                continue
        
            # calculate the scores of the old and new structure
            old_score = score(log_pos_x, log_pos_y)
            new_score = score(pos_x, pos_y)

            # keep track of each score
            scores.append(old_score)

            # if a score beats the old one, remember that structure 
            if new_score < lowest_score:
                best_x = copy.deepcopy(pos_x)
                best_y = copy.deepcopy(pos_y)
                lowest_score = copy.deepcopy(new_score)

            # probability function to determine whether a fold will be 'accepted' or not, a lower score relative to the previous configuration increases the changes of adoption
            p = math.exp(-(new_score - old_score)/(temperature * boltzmann))

            # the treshhold for acceptance varies and is randomly determined
            treshhold = random.random()
            if p < treshhold:
                pos_x = log_pos_x
                pos_y = log_pos_y

            rotation_counter += 1

            # print statement for time indication in long calculations
            if rotation_counter % 1000 == 0:
                print(f'{rotation_counter / N * 100}%')

        # the best structure is copied to a csv file and shown in a graph
        output(best_x, best_y, lowest_score)
        plot(best_x, best_y, lowest_score, scores)
        i += 1

def output(list_x, list_y, score):
    """Prints the folded string to a csv file in the Bas Terwijn style."""

    numbers = []
   
    for i in range(len(protein)-1):
        
        # new position is compared to the old
        delta_x = list_x[i+1] - list_x[i]
        delta_y = list_y[i+1] - list_y[i]

        # conversion between coordinates  and bas terwijn numbers
        if delta_x == 1:
            number = -2
        elif delta_x == -1:
            number = 2
        elif delta_y == 1:
            number = 1
        elif delta_y == -1:
            number = -1
        numbers.append(number)

    # add 0 to signal end of protein
    numbers.append(0)

    # write the list to a file
    f = open('montetest.csv', 'w')
    # f.write('amino,fold\n')
    for p, n in zip(protein, numbers):
        if n == 0:
            f.write(f'{p}, {score}\n')
        else:
            f.write(f'{p}, {n}\n')
    f.close()
    
    # the following dict logs pre-pruning unique outcomes
    resultstr = ""
    resultdict = {}
    
    with open('montedict.csv', 'r') as resultdict_csv:
        csv_reader = csv.reader(resultdict_csv)

        # Print to terminal, for testing
        for key, score in csv_reader:
            # Insert data into dict
            resultdict[key] = score 
    
    for p, n in zip(protein, numbers):
        resultstr = resultstr + str(p) + str(n)

    if resultstr not in resultdict:
        resultdict[resultstr] = score 
        print(f'unique') 
        f = open('montedict.csv', 'w')
        
        for key in resultdict:
            value = resultdict.get(key)
            f.write(f'{key},{value}\n')
        
        f.close()
    else:
        print(f'not unique') 
    
        
def score(list_x, list_y):
    """Given the coordinates of a protein string, calculate the score of the shape."""

    # list to place the 'already scored' atoms into
    coordinates = []
    directions = [[-1,0],[0,1],[1,0],[0,-1]]
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


def double(list_x, list_y):
    """Checks whether two atoms occupy the same point."""

    coordinates = []

    # see if a coordinate is already in the list, then add that coordinate to the list
    for x, y in zip(list_x, list_y):
        if [x,y] in coordinates:
            return True
        coordinates.append([x,y])
    
    return False


def plot(list_x, list_y, score, scores):
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
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 9))

    ax1.plot(list_x, list_y, '--', color='darkgrey')
    ax1.plot(red_dots_x, red_dots_y, 'or')
    ax1.plot(blue_dots_x, blue_dots_y, 'ob')
    ax1.plot(yellow_dots_x, yellow_dots_y, 'oy')
    ax1.set_title(f'Folded protein, score: {score}')

    ax2.plot(scores)
    ax2.set_title('Scores of the configurations after each rotation')
    ax2.set(xlabel='Rotation', ylabel='Score')
    
    # comment close function to show figure 
    plt.close()
    #plt.show()


def random_rotation(list_x, list_y, n):
    """Rotates the string 90 degrees to the left or right from the nth atom onwards."""

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
    main()