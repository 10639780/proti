"""
genetichelpers.py

Minor Programmeren 
Team Proti

Algorithms are implemented by calling the functions 
in this script.
"""

# import modules
import matplotlib.pyplot as plt
import random
import operator
import copy
from generalhelpers import double, directions, score_it



def initial_population(n, proti):
    """
    Initialize population for tree.
    """
    
    conformations_list = []
    return_list = []

    for i in range(2 * n):
        conformation_string = create_individual(proti)
        x, y, z = directions(conformation_string)
        score = score_it(proti=proti, list_x=x, list_y=y, list_z=z)
        conformation = (conformation_string, score)
        conformations_list.append(conformation)

    conformations_list.sort(key=operator.itemgetter(1))

    for i in range(n):
        return_list.append(conformations_list[i])

    return return_list


def make_child(conformations_list, sample_size):
    """
    Produces offspring using two parents.
    """
    child_not_made = True

    while child_not_made:

        # choose two parent from a pool of potential parents
        parent1, parent2 = choose_parents(conformations_list, sample_size)

        # use crossover between the parents to make children
        child1 = crossover(parent1, parent2)
        child2 = crossover(parent1, parent2)

        # make sure that the children have valid conformations
        if not (child1[0] == 'impossible' or child2[0] == 'impossible'):
            child_not_made = False
    
    return parent1, parent2, child1, child2


def mutate(child1, child2, p, proti):
    """
    Causes random mutations of the directions in the conformations.
    """

    # convert to lists to make the tuple mutable
    child1_list = list(child1)
    child2_list = list(child2)

    for i in range(len(child1_list)):

        # random number (0,1)
        r1 = random.random()
        r2 = random.random()

        possible = False
        if r1 < p:

            # save copy
            backup = copy.deepcopy(child1_list[i])

            # check if mutation would be possible without making the conformation invalid
            for d in random.sample(['L', 'R', 'S'], 3):
                child1_list[i] = d
                x, y, z = directions(''.join(child1_list))
                if not double(x, y, z):
                    possible = True 
                    break

            # undo the mutation if it would invalidate the conformation
            if not possible:
                child1_list[i] = backup

        possible = False
        if r2 < p:
            backup = copy.deepcopy(child2_list[i])
            for d in random.sample(['L', 'R', 'S'], 3):
                child2_list[i] = d
                x, y, z = directions(''.join(child2_list))
                if not double(x, y, z):
                    possible = True 
                    break
            if not possible:
                child2_list[i] = backup

    child1_list = ''.join(child1_list)
    child2_list = ''.join(child2_list)

    # return the mutated conformations along with their scores
    x1, y1, z1 = directions(child1_list)
    x2, y2, z2 = directions(child2_list)
    return (child1_list, score_it(proti, x1, y1, z1)), \
        (child2_list, score_it(proti, x2, y2, z2))

def create_individual(proti):
    
    conformation = ''
    direction = ['L', 'R', 'S']

    i = 1
    while i < proti.length - 1:
        available_directions = []

        for d in direction: 
            x_pos, y_pos, z_pos = directions(conformation + d)     
            if not double(x_pos, y_pos, z_pos):
                available_directions.append(d)
        
        try:
            random_direction = available_directions[random.randint(0,\
                                                    len(available_directions) - 1)]
            conformation += random_direction
        except:
            i = 0
            conformation = ''
    
        i += 1
  
    return conformation



def replace(parent1, parent2, mutation1, mutation2, conformations_list):
    """
    Replaces the fittest offspring with either the least fit parent or member 
    of the population.
    """

    # determine best child, worst parent
    best_offspring = min([mutation1, mutation2], key=operator.itemgetter(1))
    worst_parent = max([parent1, parent2], key=operator.itemgetter(1))

    # replace if child is better than parent
    if best_offspring[1] < worst_parent[1]:
        index = conformations_list.index(worst_parent)
        conformations_list[index] = best_offspring

    # otherwise look for replacement option in the population
    else:
        worst_population = max(conformations_list, key=operator.itemgetter(1))
        if best_offspring < worst_population:
            index = conformations_list.index(worst_population)
            conformations_list[index] = best_offspring

def choose_parents(conformations_list, sample_size):
    """
    Takes two random samples from the population, chooses the best of each 
    sample to be a parent.
    """

    # make pools
    pool1 = random.sample(conformations_list, sample_size)
    pool2 = random.sample(conformations_list, sample_size)

    # choose the best
    parent1 = min(pool1, key=operator.itemgetter(1))
    parent2 = min(pool2, key=operator.itemgetter(1))

    return parent1, parent2

def crossover(parent1, parent2):
    """
    Uses crossover to producre a child from two parents.
    """

    # select point of adjoinment at random
    n = random.randint(0, len(parent1[0]) - 1)

    # get the part of one parent up to the nth element
    partial_parent1 = parent1[0][:n]

    # get the part of the other parent from the n+1th elment onwards
    partial_parent2 = parent2[0][n+1:]

    # try different adjoinment directions
    for d in random.sample(['L', 'R', 'S'],3):
        
        child_string = partial_parent1 + d + partial_parent2
        x, y, z = directions(child_string)
        # return is the conformation is valid
        if not double(x, y, z):
            return child_string

    # indicate that the conformation is invalid
    return ('impossible', 1)


def genetic_plot(string, score, best_yet, proti):
    """
    Makes a graph of two lists list_x, list_y.
    """
    
    list_x, list_y, list_z = directions(string)

    # differentiate between types of atom
    red_dots_x = []
    red_dots_y = []
    blue_dots_x = []
    blue_dots_y = []
    yellow_dots_x = []
    yellow_dots_y = []

    # search through protein and place each atom in the appropiate list
    for x, y, p in zip(list_x, list_y, proti.listed):

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
    fig, (ax1, ax2) = plt.subplots(2, 1)
  
    ax1.plot(list_x, list_y, '--', color='darkgrey')
    ax1.plot(red_dots_x, red_dots_y, 'or', markersize=17)
    ax1.plot(blue_dots_x, blue_dots_y, 'ob', markersize=17)
    ax1.plot(yellow_dots_x, yellow_dots_y, 'oy', markersize=17)
    ax1.axis('equal')
    ax1.set_title(f'Folded protein of length {proti.length}, score: {score}')

    ax2.plot(best_yet)
    ax2.set_title('Lowest score per iteration')

    plt.show()
