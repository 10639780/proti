import numpy as np 
import matplotlib.pyplot as plt 
import random

protein = ['H','H','P','H','H','H','P','H']
length = len(protein)

def main():

    pos_x = [i for i in range(length)]
    pos_y = [0] * length
    rotation_counter = 0
    double_counter = 0

    while rotation_counter < 10:

        log_pos_x = pos_x
        log_pos_y = pos_y

        random_rotation(pos_x, pos_y, random.randint(0, length - 1))

        if double(pos_x, pos_y):
            double_counter += 1
            pos_x = log_pos_x
            pos_y = log_pos_y
        else:
            rotation_counter += 1

    # plot(pos_x, pos_y)

    output(pos_x, pos_y)


def output(log_x, log_y):
    """Prints the folded string to the terminal in the Bas Terwijn style."""

    numbers = []
   
    for i in range(len(protein)-1):
        
        # new position is compared to the old
        delta_x = log_x[i+1] - log_x[i]
        delta_y = log_y[i+1] - log_y[i]

        # conversion between matrix indeces and bas terwijn numbers
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


    f = open('output.csv', 'w')
    f.write('amino,fold\n')
    for p, n in zip(protein, numbers):
        f.write(f'{p}, {n}\n')
    total_score = score(log_x, log_y)
    f.write(f'score,{total_score}') 
    f.close()
  

def score(list_x, list_y):

    coordinates = []
    directions = [[-1,0],[0,1],[1,0],[0,-1]]
    score = 0

    for i in range(length):

        xi = list_x[i]
        yi = list_y[i]
        pi = protein[i]

        if not pi == 'P':

            for d in directions:

                for j in range(len(coordinates)):

                    if [xi + d[0], yi + d[1]] == coordinates[j] and not(xi + d[0] == list_x[i-1] and yi + d[1] == list_y[i-1]):
                        
                        if pi == 'H':
                            if protein[j] == 'H' or protein[j] == 'C':
                                score += -1

                        if pi == 'C':
                            if protein[j] == 'C':
                                score += -5
                            if protein[j] == 'H':
                                scpre += -1 
        
        coordinates.append([xi, yi])
    
    return score


def double(list_x, list_y):
    """Checks whether two atoms occupy the same point."""

    coordinates = []

    for x, y in zip(list_x, list_y):
        if [x,y] in coordinates:
            return True
        coordinates.append([x,y])
    
    return False


def plot(list_x, list_y):
    """Makes a graph of two lists list_x, list_y."""
    plt.plot(list_x, list_y, '--')
    plt.plot(list_x, list_y, 'o')
    plt.show()


def random_rotation(x, y, n):
    """Rotates the string 90 degrees to the left from the nth atom onwards."""

    rotation_point_x = x[n]
    rotation_point_y = y[n]

    # left or right rotation are randomly chosen
    p = random.random()

    # calculates the new positions for the remainder of the string using the equations from a 2D rotation matrix
    for i in range(n + 1, length):
       
        relative_x = x[i] - rotation_point_x
        relative_y = y[i] - rotation_point_y

        if p > 0.5:
            # rotate left
            x[i] = rotation_point_x - relative_y
            y[i] = rotation_point_y + relative_x
        else:
            # rotate right
            x[i] = rotation_point_x + relative_y
            y[i] = rotation_point_y - relative_x
 

if __name__ == "__main__":
    main()