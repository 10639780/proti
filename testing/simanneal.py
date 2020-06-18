import numpy as np
from random import randint

DIMENSION = 2
PROTEIN = 'HHPH'
LENGTH = len(PROTEIN)

def grid_structure():
    grid_size = 2 * LENGTH + 1
    grid = [[0 for i in range(grid_size)] for j in range(grid_size)]

    return grid

def atom_positions():
    grid = grid_structure()
    start_amino = PROTEIN[0]
    grid[LENGTH][LENGTH] = start_amino
    
    x_coord = [LENGTH]
    y_coord = [LENGTH]
    protein_build = [start_amino]

    first_step = True

    for i in range(LENGTH - 1):
        next_amino = PROTEIN[i + 1]


        if first_step:

            step = valid_step(x_coord[-1], y_coord[-1], grid)
            grid[step[0]][step[1]] = next_amino
            first_step = False
            continue

        step = valid_step(x_coord[-1], y_coord[-1], grid)
        grid[step[0]][step[1]] = next_amino

    return grid


def step_direction(current_x, current_y):
    directions = [-2, -1, 1, 2]
    random_move = randint(0, len(directions) - 1)

    step = directions[random_move]
    
    if step in [-1, 1]:
        new_x = current_x + step
        new_y = current_y

    else:
        new_x = current_x
        new_y = current_y + step

    return new_x, new_y

def valid_step(current_x, current_y, grid):
    new_coords = step_direction(current_x, current_y)

    while grid[new_coords[0]][new_coords[1]] == 0:
        new_coords = step_direction(current_x, current_y)

    return new_coords
    
    
def display_grid(grid):

    # grid_print = '\n'.join([''.join(['{:5}'.format(item) for item in row]) for row in grid])

    grid_print = np.matrix(grid)
    return grid_print


if __name__ == "__main__":
    g = atom_positions()
    print(display_grid(g))
