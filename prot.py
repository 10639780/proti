"""
prot.py

Minor Programmeren
Team Proti

Folds a protein string randomly and determines the resulting stability.
I don't know how to make that if __main__ stuff  work without losing access the all variables so all functions are now above and the main code is at the bottom.
I used x,y for many of the coordinate placeholders, but since all values are stored in a matrix i,j would be more appropiate, but since they are already used in the loops and I'm to lazy rewrite everything in termns of m,n I chose to stick with x,y which hopefully is not too confusing.
"""
import numpy as np 
from random import randint
import sys
import matplotlib.pyplot as plt 


def possible_moves(x,y,_type):
    """Determines which moves are possible from a given position. 
    If possible function will also determine the points gained by making that move."""

    # list to hold all valid moves
    moves = []
   
    # check if the target tile is still empty, and whether that tile doesn't have 4 filled neighbours already, which would result in the string getting wrappped into itself
    if grid[x-1][y] == 0 and neighbour_count(x-1,y) < 4:
        score = 0
        # if the protein to be placed is hydrophobe, determine what the score would be if placed in this position
        if _type == 'H':
            score = points(x-1,y,x,y)
        # add to list of valid moves
        moves.append([x-1, y, score])
    
    if grid[x][y+1] == 0 and neighbour_count(x,y+1) < 4:
        score = 0
        if _type == 'H':
            score = points(x,y+1,x,y)
        moves.append([x, y+1, score])

    if grid[x+1][y] == 0 and neighbour_count(x+1,y) < 4:
        score = 0
        if _type == 'H':
            score = points(x+1,y,x,y)
        moves.append([x+1, y, score])

    if grid[x][y-1] == 0 and neighbour_count(x,y-1) < 4:
        score = 0
        if _type == 'H':
            score = points(x,y-1,x,y)
        moves.append([x, y-1, score])

    return moves


def neighbour_count(x,y):
    """Counts how many non-empty neighbours a given position has."""
    counter = 0

    # goes to all 4 neighbouring tiles, if not empty, add 1 to counter
    if grid[x-1][y] != 0:
        counter += 1
    
    if grid[x][y+1] != 0:
        counter += 1

    if grid[x+1][y] != 0:
        counter += 1

    if grid[x][y-1] != 0:
        counter += 1

    return counter


def display(score):
    """Gives a graphical display of the protein string after folding."""

    # make two coordinate list for both tpyes of protein
    red_dots_x = []
    red_dots_y = []
    blue_dots_x = []
    blue_dots_y = []

    # iterate through the grid and add coordinates of the proteins found
    for i in range(2 * length + 1):
        for j in range(2 * length + 1):

            if grid[i][j] == 'H':
                red_dots_x.append(i)
                red_dots_y.append(j)
            
            if grid[i][j] == 'P':
                blue_dots_x.append(i)
                blue_dots_y.append(j)

    # make a colorcoded graph
    plt.plot(red_dots_x,red_dots_y,'or',blue_dots_x,blue_dots_y,'ob')
    plt.title('Folded protein string')
    plt.text(1, 1, f'Score: {score}')
    plt.xlim(0, 2 * length + 1)
    plt.ylim(0, 2 * length + 1)
    plt.show()


def points(x, y, old_x, old_y):
    """Determines how many points are gained by making a move to x, y, coming from old_x, old_y."""
    
    counter = 0

    # add 1 point if the target position is a hydrophobe and is not the previous protein in the string
    if grid[x-1][y] == 'H' and not (x-1 == old_x and y == old_y):
        counter += 1
    
    if grid[x][y+1] == 'H' and not (x == old_x and y+1 == old_y):
        counter += 1

    if grid[x+1][y] == 'H' and not (x+1 == old_x and y == old_y):
        counter += 1

    if grid[x][y-1] == 'H' and not (x == old_x and y-1 == old_y):
        counter += 1

    return counter


def make_grid(x_width, y_width):
    """Creates a grid of specified size."""

    grid = []
    for i in range(x_width):
        row = []
        for j in range(y_width):
            row.append(0)
        grid.append(row)

    return grid

def output(log_x, log_y, protein_string):
    """Prints the folded string to the terminal in the Bas Terwijn style."""

    numbers = []
   
    for i in range(len(protein_string)-1):
        
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

    numbers.append(0)

    
    for p, n in zip(protein_string, numbers):
        print(f'{p}, {n}')
    
  
# enter protein string
protein_string = ['H','H','P','H','H','H','P','H']
length = len(protein_string)

# some initialization variables
total_score = 0
first_loop = True

# create grid and set begin position of protein string
grid = make_grid(2*length+1,2*length+1)
current_position = [length , length, 0]

# log of the moves made, needed for correct output
log_x = []
log_y = []

# loop over the entire protein string
for protein in protein_string:

    # first protein is just placed in the middle and nothing more
    if first_loop:
        grid[current_position[0]][current_position[1]] = protein
        log_x.append(current_position[0])
        log_y.append(current_position[1])
        first_loop = False
        
        continue

    # determine all valid directions
    moves = possible_moves(current_position[0],current_position[1],protein)

    # make sure the next protein has somewhere to go
    if len(moves) == 0:
        print ('Protein is stuck.') 
        sys.exit()

    # choose one direction at random (could just as well choose the one with the highest score, which would be the third element in the list's elements)
    current_position = moves[randint(0,len(moves) - 1)]

    # keep track of the score associated with the move made
    total_score += current_position[2]

    # add protein to the grid
    grid[current_position[0]][current_position[1]] = protein
    log_x.append(current_position[0])
    log_y.append(current_position[1])


# display the result (comment out the ones you dont need)
print(np.array(grid)) # grid in terminal
print(f'Score: {total_score}') # final score
output(log_x, log_y, protein_string) # Bas Terwijn output
display(total_score) # graphical display



