"""
genetic.py

Minor Programmeren
Team Proti

Genetic protein folding algorithm based on 'An Efficient Genetic Algorithm for Predicting
Protein Tertiary Structures in the 2D HP Model' by Thang N. Bui and Gnanasekaran Sundarraj 
and on 'Genetic Algorithm for Predicting Protein Folding in the 2D HP Model' by Eyal Halm.
"""

from itertools import groupby
import random
import operator
import copy
import matplotlib.pyplot as plt
from progress.bar import Bar
import timeit
import math

# protein = 'HHPHHHPHPHHHPH' # 14, opt -6

# protein = 'HPHPPHHPHPPHPHHPPHPH' # 20, otp -9

# protein = 'PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP' # 36, opt -14

# protein = 'HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH' # 50, opt -21

# protein = 'PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP' #36

# protein = 'CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC' # 36

# protein = 'HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH' # 50

# protein = 'HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH' # 50

length = len(protein)

def main():
    # some constants
    population_size = 500
    sample_size = round(population_size * 0.1)
    N = 100000
    _max = 0.4
    _min = 0.2
    scale = N/20
    center = 9 * N / 40
    bar = Bar('Progress', max=N/1000)
    best_yet = []

    # initialization
    conformations_list = initial_population(population_size)
    initial_time = timeit.default_timer()
    # iterative part
    for i in range(N): 
        mutation_chance = _max / (1 + math.exp((i - center )/ scale)) + _min
        parent1, parent2, child1, child2  = make_child(conformations_list, sample_size)    
        mutation1, mutation2 = mutate(child1, child2, mutation_chance)
        replace(parent1, parent2, mutation1, mutation2, conformations_list)
        best_yet.append(min([parent1, parent2, mutation1, mutation2], key=operator.itemgetter(1))[1])
        if i % 1000 == 0:
            bar.next()

    # output
    best_conformation =  min(conformations_list, key=operator.itemgetter(1))
    bar.finish()
    stop = timeit.default_timer()
    initialization_time = initial_time - start
    runtime = stop - initial_time
    total_time = stop - start
    print(f'Length: {length}')
    print(f'Score: {best_conformation[1]}')
    print(f'Initialization: {initialization_time}')
    print(f'Runtime: {runtime}')
    print(f'Total: {total_time}')
    print(f'Conformation: {best_conformation[0]}')
    plot(best_conformation[0], best_conformation[1], best_yet)


def replace(parent1, parent2, mutation1, mutation2, conformations_list):

    best_offspring = min([mutation1, mutation2], key=operator.itemgetter(1))
    worst_parent = max([parent1, parent2], key=operator.itemgetter(1))

    if best_offspring[1] < worst_parent[1]:
        index = conformations_list.index(worst_parent)
        conformations_list[index] = best_offspring
    else:
        worst_population = max(conformations_list, key=operator.itemgetter(1))
        if best_offspring < worst_population:
            index = conformations_list.index(worst_population)
            conformations_list[index] = best_offspring

def mutate(child1, child2, p):
    
    child1_list = list(child1)
    child2_list = list(child2)

    for i in range(len(child1_list)):
        possible = False

        r1 = random.random()
        r2 = random.random()

        possible = False
        if r1 < p:
            backup = copy.deepcopy(child1_list[i])

            for d in random.sample(['L', 'R', 'S'],3):
                
                child1_list[i] = d
                if not double(''.join(child1_list)):
                    possible = True 
                    break
            
            if not possible:
                child1_list[i] = backup

        possible = False
        if r2 < p:
            backup = copy.deepcopy(child2_list[i])

            for d in random.sample(['L', 'R', 'S'],3):
                
                child2_list[i] = d
                if not double(''.join(child2_list)):
                    possible = True 
                    break
            
            if not possible:
                child2_list[i] = backup

    child1_list = ''.join(child1_list)
    child2_list = ''.join(child2_list)

    return (child1_list,score_func(child1_list)), (child2_list, score_func(child2_list))


def make_child(conformations_list, sample_size):

    child_not_made = True
    child2_not_made = True
    while child_not_made:
        parent1, parent2 = choose_parents(conformations_list, sample_size)
        child1 = crossover(parent1, parent2)
        child2 = crossover(parent1, parent2)
        if not (child1[0] == 'impossible' or child2[0] == 'impossible'):
            child_not_made = False
    
    return parent1, parent2, child1, child2


def crossover(parent1, parent2):
    n = random.randint(0, len(parent1[0]) - 1)

    partial_parent1 = parent1[0][:n]
    partial_parent2 = parent2[0][n+1:]

    for d in random.sample(['L', 'R', 'S'],3):
    
        child_string = partial_parent1 + d + partial_parent2
        if not double(child_string):
            return child_string

    return ('impossible', 1)


def choose_parents(conformations_list, sample_size):

    pool1 = random.sample(conformations_list, sample_size)
    pool2 = random.sample(conformations_list, sample_size)

    parent1 = min(pool1, key=operator.itemgetter(1))
    parent2 = min(pool2, key=operator.itemgetter(1))

    return parent1, parent2


def initial_population(n):
    
    conformations_list = []
    return_list = []

    for i in range(2 * n):
        conformation_string = create_individual()
        score = score_func(conformation_string)
        conformation = (conformation_string, score)
        conformations_list.append(conformation)

    conformations_list.sort(key=operator.itemgetter(1))

    for i in range(n):
        return_list.append(conformations_list[i])

    return return_list

def create_individual():
    
    conformation = ''
    directions = ['L', 'R', 'S']

    i = 1
    while i < length - 1:
        available_directions = []

        for d in directions:      
            if not double(conformation + d):
                available_directions.append(d)
        
        try:
            random_direction = available_directions[random.randint(0, len(available_directions) - 1)]
            conformation += random_direction
        except:
            i = 0
            conformation = ''
    
        i += 1
  
    return conformation

def direction_to_xy(string):
    """Converts a series of string with directions like ['L', 'R'] to lists with xy positions."""
    # first two aminos already placed
    pos_x = [0,1]
    pos_y = [0,0]

    # go over every node
    for s in string:
    
        # previous direction is determined
        delta_x = pos_x[-1] - pos_x[-2]
        delta_y = pos_y[-1] - pos_y[-2]

        # rotation matrices used to turn into the desired direction
        if s == 'S':
            pos_x.append(pos_x[-1] + delta_x)
            pos_y.append(pos_y[-1] + delta_y)
       
        elif s == 'L':
            pos_x.append(pos_x[-1] - delta_y)
            pos_y.append(pos_y[-1] + delta_x)
         
        elif s == 'R':
            pos_x.append(pos_x[-1] + delta_y )
            pos_y.append(pos_y[-1] - delta_x)
        
    return pos_x, pos_y


def double(string):

    """Checks whether two atoms occupy the same point."""
    list_x, list_y = direction_to_xy(string)
    coordinates = []
    for x, y in zip(list_x, list_y):
        coordinates.append((x,y))

    if len(coordinates) == len(set(coordinates)):
        return False
    
    return True

def score_func(string):
    """Given the coordinates of a protein string, calculate the score of the shape."""
    list_x, list_y = direction_to_xy(string)
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


def plot(string, score, best_yet):
    """Makes a graph of two lists list_x, list_y."""
    
    list_x, list_y = direction_to_xy(string)

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
    fig, (ax1, ax2) = plt.subplots(2, 1)
  
    ax1.plot(list_x, list_y, '--', color='darkgrey')
    ax1.plot(red_dots_x, red_dots_y, 'or', markersize=17)
    ax1.plot(blue_dots_x, blue_dots_y, 'ob', markersize=17)
    ax1.plot(yellow_dots_x, yellow_dots_y, 'oy', markersize=17)
    ax1.axis('equal')
    ax1.set_title(f'Folded protein of length {length}, score: {score}')

    ax2.plot(best_yet)
    ax2.set_title('Lowest score per iteration')

    plt.show()


if __name__ == "__main__":
    start = timeit.default_timer()
    main()
    