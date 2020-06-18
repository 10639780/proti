""" 
bab.py

Minor Programmeren
Team Proti

Attempt to implement a branch and bound protein folding algorithm as described by Mao Chen and Wen-Qi Huang 
in 'Branch and Bound Algorithm for the Protein Folding Problem in the HP Lattice Model'.
Basis for depth first structure from Bas Terwijn's lecture.
"""
import copy
import timeit
import matplotlib.pyplot as plt 
import random
import queue
from mpl_toolkits import mplot3d

# bunch of test strings

# protein = ['H', 'H', 'P', 'H']
# protein = ['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H'] # official 8, mc -3
# protein = ['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'P', 'H'] # official 14, mc -6
protein = ['H', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'H', 'P', 'P', 'H', 'P', 'H', 'H', 'P', 
  'P', 'H', 'P', 'H'] # 20, opt -9, bench 0.21, -8
# protein = ['H', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'P', 'P', 
#   'H', 'P', 'P', 'H', 'P', 'P', 'H','H'] # 24,  opt -9
# protein = ['P', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 
#   'P', 'H', 'H', 'P', 'P', 'P', 'P', 'H', 'H'] # 25, opt -8
# protein = ['P', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'P', 'H', 'H', 
#   'H', 'H', 'H', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'P'] # 36,  opt -14, bench 4.6, -12
# protein = ['P', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'P', 
#   'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'P', 'P', 'P', 'P', 'P', 'H', 'H', 
#   'P', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'H', 'H', 'H', 'H'] # 48, opt -23
# protein = ['P', 'P', 'H', 'P', 'P', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 
#   'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 
#   'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H']  # 50, opt -21 low -16
# protein = ['H', 'H', 'P', 'H', 'C', 'H', 'P', 'C', 'P', 'C', 'H'] #  mc -3
# protein = ['C', 'P', 'P', 'C', 'H', 'P', 'P', 'C', 'H', 'P', 'P', 'C', 'P', 'P', 'H', 'H', 'H', 
    # 'H', 'H', 'H', 'C', 'C', 'P', 'C', 'H', 'P', 'P', 'C', 'P', 'C', 'H', 'P', 'P', 'H', 'P', 'C'] # 36, mc -35
# protein = ['H', 'C', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'C', 'H', 'C', 'H', 'P', 'H', 'P', 'P', 'P', 
#     'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 
#     'H', 'H', 'H', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'H']  # 50 low -28
# protein = ['H', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'P', 'P',
#      'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 
#         'H', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H'] # 50, low -19
# protein = ['H', 'C', 'P', 'H', 'P', 'H', 'P', 'H', 'C', 'H', 'H', 'H', 'H', 'P', 'C', 'C', 'P', 'P',
#      'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'C', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P',
#      'H', 'H', 'H', 'H', 'C', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H']  # 50 low -27
# protein = ['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'P', 'P',
#      'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 
#      'H', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 
#      'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'] # 64, opti -42
# protein = ['P', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'H', 'H', 'H', 'P', 'P', 'H', 'H', 'H', 'P', 'H', 
#     'H', 'P', 'H', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'H', 'H', 'H',
#      'H', 'H', 'H', 'P', 'P', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P', 'P',
#       'H', 'P', 'H', 'H', 'P', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'P', 'H',
#        'H', 'H', 'P', 'H', 'H', 'P', 'H', 'P', 'P', 'H', 'P', 'H', 'H', 'H', 'P', 'P', 'P', 'P', 'P',
#         'P', 'H', 'H', 'H'] # 100, opt -50
# protein = ['H', 'H', 'H', 'H', 'P', 'P', 'P', 'P', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
#      'H', 'P', 'P', 'P', 'P', 'P', 'P', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'P',
#       'P', 'P', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'P', 'P', 'H', 'H', 'H',
#        'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'H',
#         'H', 'P', 'P', 'H', 'P', 'H'] # 85, opt -53
length = len(protein)

# set theoretical lower bound on score
even = protein[::2]
odd = protein[1::2]
min_score = 4 * max([- even.count('H') - 5 * even.count('C'), - odd.count('H') - 5 * odd.count('C')])

def main():

    # specifications for depth first tree building
    depth = length - 2
    q = queue.Queue()
    q.put('')
    final_configurations = []
    
    # keep track of scores per substring
    lowest_score_k = {}
    all_scores_k = {}
    lowest_score = 0
    
    # (0,1) probabilities of pruning a path, lower is more exact but less fast
    p1 = 1
    p2 = 1

    # set inital values
    for i in range(length + 1):
        lowest_score_k[i] = 0
        all_scores_k[i] = [0]

    # create a breadth first tree
    while not q.empty():

        state = q.get()

        # if all aminos are placed, put the string in a list
        if len(state) == depth and not double(state):
            final_configurations.append(state)
      
        if len(state) < depth:
            for i in ['L', 'R', 'S', 'U', 'D']:

                # substring
                child = copy.deepcopy(state) 

                # string after potentialy placing the next amino
                child += i 

                # discard the string folding into themselves
                if double(child):
                    continue
                
                # identify how for into the string it is
                k = len(child) + 1
                
                # P's are always placed, rest have some conditions
                if not protein[k] == 'P':
                    
                    # score if placed 
                    score = score_func(child)

                    # min score to get from remaining aminos
                    possible_score = possible_score_func(protein[k+1:], score)

                    if score + possible_score > lowest_score:
                        continue

                    # random number between 0 and 1
                    r = random.random()

                    # avergage of all strings of the same length
                    average_score_k = sum(all_scores_k[k]) / len(all_scores_k[k])

                    # conditions for pruning
                    if score > average_score_k and r < p1:
                        continue
                    elif (average_score_k >= score and score > lowest_score_k[k]) and r < p2:
                        continue
                    
                    # add to tree
                    q.put(child) 
                    all_scores_k[k].append(score)

                    if score < lowest_score_k[k]:
                        lowest_score_k[k] = copy.deepcopy(score)
                    
                    if score < lowest_score:
                        lowest_score = copy.deepcopy(score)
                        
                else:
                    q.put(child) 

    
    lowest_score = 0

    # weed out the best configuration from the remaining strings
    for c in final_configurations:
        if score_func(c) < lowest_score:
            best_config = copy.deepcopy(c)
            lowest_score = copy.deepcopy(score_func(c))

    # plot the result
    stop = timeit.default_timer()
    # print(f'Strings made: {len(final_configurations)}')
    print(f'Length: {length}')
    print(f'Score: {lowest_score}')
    print(f'Runtime: {stop - start}')
    plot(best_config, lowest_score)

def possible_score_func(nodes_to_visit, partial_score):
    """Calculates the best score the remaining bit of the protien can acquire."""

    possible_score = 0

    # an H atom can get at most -2 and a C atom at best -10
    for i, n in enumerate(nodes_to_visit):

        if i == len(nodes_to_visit) - 1:
            if n == 'H':
                possible_score += -5
            if n == 'C':
                possible_score += -25
            continue

        if n == 'H':
             possible_score += -4
        if n == 'C':
            possible_score += -20

    if possible_score + partial_score < min_score:
        possible_score = min_score - partial_score
        
    return possible_score

def plot(string, score):
    """Makes a graph of two lists list_x, list_y."""
    
    list_x, list_y, list_z = direction_to_xyz(string)

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
    ax.plot(red_dots_x, red_dots_y, red_dots_z, 'or', markersize=17)
    ax.plot(blue_dots_x, blue_dots_y, blue_dots_z, 'ob', markersize=17)
    ax.plot(yellow_dots_x, yellow_dots_y, yellow_dots_z, 'oy', markersize=17)
    ax.set_title(f'Folded protein, score: {score}')
    # plt.savefig("monte3dfig.png")
    plt.show()

def score_func(string):
    """Given the coordinates of a protein string, calculate the score of the shape."""

    # list to place the 'already scored' atoms into
    coordinates = []
    directions = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]
    score = 0
    list_x, list_y, list_z = direction_to_xyz(string)
    length = len(list_x)

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
def direction_to_xyz(string):
    """Converts a series of string with directions like ['L', 'R'] to lists with xy positions."""

    pos_x = [0,1]
    pos_y = [0,0]
    pos_z = [0,0]
    # go over every node
    for s in string:
    
        # previous direction is determined
        delta_x = pos_x[-1] - pos_x[-2]
        delta_y = pos_y[-1] - pos_y[-2]
        delta_z = pos_z[-1] - pos_z[-2]

        # rotation matrices used to turn into the desired direction
        if s == 'S':
            pos_x.append(pos_x[-1] + delta_x)
            pos_y.append(pos_y[-1] + delta_y)
            pos_z.append(pos_z[-1] + delta_z)
       
        elif s == 'L':
            pos_x.append(pos_x[-1] - delta_y)
            pos_y.append(pos_y[-1] + delta_x)
            pos_z.append(pos_z[-1] + delta_z)
         
        elif s == 'R':
            pos_x.append(pos_x[-1] + delta_y )
            pos_y.append(pos_y[-1] - delta_x)
            pos_z.append(pos_z[-1] + delta_z)
        
        elif s == 'U':
            pos_x.append(pos_x[-1] - delta_z)
            pos_y.append(pos_y[-1] + delta_y)
            pos_z.append(pos_z[-1] + delta_x)

        elif s == 'D':
            pos_x.append(pos_x[-1] + delta_z)
            pos_y.append(pos_y[-1] + delta_y)
            pos_z.append(pos_z[-1] - delta_x)
 
   
    return pos_x, pos_y, pos_z


def double(string):
    """Checks whether two atoms occupy the same point."""

    list_x, list_y, list_z = direction_to_xyz(string)

    coordinates = []


    # see if a coordinate is already in the list, then add that coordinate to the list
    for x, y, z in zip(list_x, list_y, list_z):
        if [x, y, z] in coordinates:
            return True
        coordinates.append([x, y, z])
    
    return False


if __name__ == "__main__":
    start = timeit.default_timer()
    main()