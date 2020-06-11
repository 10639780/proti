"""
prot.py

Minor Programmeren
Team Proti

Folds a protein string randomly and determines the resulting stability.
I don't know how to make that if __main__ stuff work without losing access the all variables so all functions are now above and the main code is at the bottom.
I used x,y for many of the coordinate placeholders, but since all values are stored in a matrix i,j would be more appropiate, but since they are already used in the loops and I'm to lazy rewrite everything in termns of m,n I chose to stick with x,y which hopefully is not too confusing.
"""
import numpy as np 
from random import randint
import sys
import matplotlib.pyplot as plt 
import csv

def possible_moves(x,y,_type):
    """Determines which moves are possible from a given position. 
    If possible function will also determine the points gained by making that move."""

    score = 0
    moves = []
    directions = [[x-1,y], [x,y+1], [x+1,y], [x,y-1]]

    # check if the target tile is still empty, and whether that tile doesn't have 4 filled neighbours already, which would result in the string getting wrappped into itself
    for direction in directions:
        if grid[direction[0]][direction[1]] == 0 and neighbour_count(direction[0],direction[1]) < 4:
            # calculate points gained by making this move
            score = points(direction[0],direction[1],x,y, _type)
            # place move and associated score in the moves list
            moves.append([direction[0],direction[1], score])

    return moves


def neighbour_count(x,y):
    """Counts how many non-empty neighbours a given position has."""

    counter = 0
    directions = [[x-1,y], [x,y+1], [x+1,y], [x,y-1]]

    for direction in directions:
        # check whether the neighbouring tile is occupied, if so add one to counter
        if grid[direction[0]][direction[1]] != 0:
            counter += 1

    return counter


def display(score):
    """Gives a graphical display of the protein string after folding."""

    # make two coordinate list for both tpyes of protein
    red_dots_x = []
    red_dots_y = []
    blue_dots_x = []
    blue_dots_y = []
    yellow_dots_x = []
    yellow_dots_y = []

    # iterate through the grid and add coordinates of the proteins found
    for i in range(size):
        for j in range(size):

            if grid[i][j] == 'H':
                red_dots_x.append(i)
                red_dots_y.append(j)
            
            if grid[i][j] == 'P':
                blue_dots_x.append(i)
                blue_dots_y.append(j)
            
            if grid[i][j] == 'C':
                yellow_dots_x.append(i)
                yellow_dots_y.append(j)

    # make a colorcoded graph
    plt.plot(log_x,log_y,'--')
    plt.plot(yellow_dots_x, yellow_dots_y, 'oy', label='Cysteine')
    plt.plot(red_dots_x, red_dots_y, 'or', label='Hydrophobe')
    plt.plot(blue_dots_x,blue_dots_y,'ob', label='Polar')
   
    plt.title(f'Folded protein string, score: {score}')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.text('left', 'bottom', f'Score: {score}')

    plt.show()


def points(x, y, old_x, old_y, _type):
    """Determines how many points are gained by making a move to x, y, coming from old_x, old_y."""
    
    counter = 0
    HH_HC_score = -1
    CC_score = -5
    directions = [[x-1,y], [x,y+1], [x+1,y], [x,y-1]]

    #  add to the score if the target position is occupied by an H or C atom
    if _type == 'H':

        for direction in directions:
            if (grid[direction[0]][direction[1]] == 'H' or grid[direction[0]][direction[1]] == 'C') and not (direction[0] == old_x and direction[1] == old_y):
                counter += HH_HC_score            

    # for C atom differentiate between the cases of neighbouring with a C or H atom
    if _type == 'C':

        for direction in directions:
            if grid[direction[0]][direction[1]] == 'H' and not (direction[0] == old_x and direction[1] == old_y):
                counter += HH_HC_score   

        for direction in directions:
            if grid[direction[0]][direction[1]] == 'C' and not (direction[0] == old_x and direction[1] == old_y):
                counter += CC_score  

    return counter


def make_grid(size):
    """Creates a grid of specified size."""

    grid = []
    for i in range(size):
        row = []
        for j in range(size):
            row.append(0)
        grid.append(row)

    return grid

def output(log_x, log_y, protein):
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
    
    # write output to a csv file
    for p, n in zip(protein, numbers):
        print(f'{p}, {n}')
        f.write(f'{p}, {n}\n')

def lowest_score(moves):
    """Determines which move has te highest score."""

    best_moves = []
    best_score = 0

    for move in moves:
        
        # take the score assiciated with each move
        score = move[2]
    
        # if score is equal to the best add it to the list of best moves
        if score == best_score:
            best_moves.append(move)
        # when a lower score is found, empty the list and place the new score in there
        elif score < best_score:
            best_moves = [move]
            best_score = score 
        
    return best_moves


# enter protein string
# protein = ['H','H', 'P','H','H','H','P','H']
# protein = ['H', 'C', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'C', 'H', 'C', 'H', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'H', 'H', 'H', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'H']
protein = ['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'P', 'H']
# protein = ['H', 'H', 'H', 'P', 'C', 'C', 'H', 'P', 'C', 'C', 'P', 'H']
protein = ['H', 'P', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'P']
protein = ['P', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'P', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'P'] # official 36
protein = ['P', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'P']

length = len(protein)
size = 2 * length + 1

# some initialization variables
total_score = 0
first_loop = True

# create grid and set begin position of protein string
grid = make_grid(size)
current_position = [length , length, 0]

# log of the moves made, needed for correct output
log_x = []
log_y = []

# loop over the entire protein string
for atom in protein:

    # first atom is just placed in the middle and nothing more
    if first_loop:
        grid[current_position[0]][current_position[1]] = atom
        log_x.append(current_position[0])
        log_y.append(current_position[1])
        first_loop = False
        continue

    # determine all valid directions
    moves = possible_moves(current_position[0],current_position[1],atom)

    # make sure the next atom has somewhere to go
    if len(moves) == 0:
        print ('Protein is stuck.') 
        sys.exit()

    # # choose one direction at random (could just as well choose the one with the highest score, which would be the third element in the list's elements)
    # current_position = moves[randint(0,len(moves) - 1)]

    # out of the valid possible moves, choose the one(s) with the highest score
    best_moves = lowest_score(moves)

    # out of the highest scored moves, choose one at random and update the position
    current_position = best_moves[randint(0, len(best_moves) - 1)]
    
    # keep track of the score associated with the move being made
    total_score += current_position[2]

    # add atom to the grid
    grid[current_position[0]][current_position[1]] = atom
    log_x.append(current_position[0])
    log_y.append(current_position[1])

# optional displays of the result (comment out the ones you dont need)
# print(np.array(grid)) # grid in terminal
print(f'Score: {total_score}') # final score
display(total_score) # graphical display

# correct output of results (don't comment out)
f = open('output.csv', 'w')
f.write('amino,fold\n')
output(log_x, log_y, protein) 
f.write(f'score,{total_score}') 
f.close()
