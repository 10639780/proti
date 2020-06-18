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
from progress.bar import Bar

# bunch of test stringsw

# protein = ['H', 'H', 'P', 'H']
protein = ['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H'] # official 8, mc -3
protein = ['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'P', 'H'] # official 14, mc -6
protein = ['H', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'H', 'P', 'P', 'H', 'P', 'H', 'H', 'P', 
  'P', 'H', 'P', 'H'] # 20, opt -9, bench 0.21, -8
# protein = ['H', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'P', 'P', 
#   'H', 'P', 'P', 'H', 'P', 'P', 'H','H'] # 24,  opt -9
# protein = ['P', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 
#   'P', 'H', 'H', 'P', 'P', 'P', 'P', 'H', 'H'] # 25, opt -8
protein = ['P', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'P', 'H', 'H', 
  'H', 'H', 'H', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'P'] # 36,  opt -14, bench 4.6, -12
# protein = ['P', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'P', 
#   'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'P', 'P', 'P', 'P', 'P', 'H', 'H', 
#       'P', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'H', 'H', 'H', 'H'] # 48, opt -23
# protein = ['P', 'P', 'H', 'P', 'P', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 
#   'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 
#       'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H']  # 50, opt -21 low -16
# # protein = ['H', 'H', 'P', 'H', 'C', 'H', 'P', 'C', 'P', 'C', 'H'] #  mc -3
# protein = ['C', 'P', 'P', 'C', 'H', 'P', 'P', 'C', 'H', 'P', 'P', 'C', 'P', 'P', 'H', 'H', 'H', 
#     'H', 'H', 'H', 'C', 'C', 'P', 'C', 'H', 'P', 'P', 'C', 'P', 'C', 'H', 'P', 'P', 'H', 'P', 'C'] # 36, mc -35
# protein = ['H', 'C', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'C', 'H', 'C', 'H', 'P', 'H', 'P', 'P', 'P', 
#     'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 
#       'H', 'H', 'H', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'H']  # 50 low -28

# protein = ['H', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'P', 'P',
#      'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 
#         'H', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H'] # 50, low -19
# protein = ['H', 'C', 'P', 'H', 'P', 'H', 'P', 'H', 'C', 'H', 'H', 'H', 'H', 'P', 'C', 'C', 'P', 'P',
#      'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'C', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P',
#      '    H', 'H', 'H', 'H', 'C', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H']  # 50 low -27
# protein = ['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'P', 'P',
#      'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 
#         'H', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 
#             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H'] # 64, opti -42
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
# protein = ['P', 'P', 'H', 'H', 'H', 'P', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'P', 'P', 'H', 'H', 'H',
#      'H', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'P', 'P', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 
#         'H', 'H', 'P', 'P', 'P', 'P', 'H', 'H', 'H', 'H', 'H', 'H', 'P', 'H', 'H', 'P', 'H', 'P'] # 60, opt -36
# protein = ['H','H', 'P','H','H','H','P','H']
# protein = ['H', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'H', 'P', 'P', 'H', 'P', 'H', 'H', 'P', 
#   'P', 'H', 'P', 'H'] # 20, opt -9, bench 0.21, -8
# protein = ['P', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'P', 'H', 'H', 
#   'H', 'H', 'H', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'P'] # 36,  opt -14, bench 4.6, -12
# protein = ['P', 'P', 'H', 'P', 'P', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 
#   'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 
#       'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H']  # 50, opt -21 low -16
# protein = ['P', 'P', 'C', 'H', 'H', 'P', 'P', 'C', 'H', 'P', 'P', 'P', 'P', 'C', 'H', 'H', 'H', 'H', 'C', 'H', 'H', 'P', 'P', 'H', 'H', 'P', 'P', 'P', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'P'] # 36, opt 0
# protein = ['C', 'P', 'P', 'C', 'H', 'P', 'P', 'C', 'H', 'P', 'P', 'C', 'P', 'P', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'P', 'C', 'H', 'P', 'P', 'C', 'P', 'C', 'H', 'P', 'P', 'H', 'P', 'C'] # 36, opt 0
# protein = ['H', 'C', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'C', 'H', 'C', 'H', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'C', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'H', 'H', 'H', 'C', 'C', 'H', 'C', 'H', 'C', 'H', 'C', 'H', 'H'] # 50, opt 0
# protein = ['H', 'C', 'P', 'H', 'P', 'H', 'P', 'H', 'C', 'H', 'H', 'H', 'H', 'P', 'C', 'C', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'C', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 'C', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H'] # 50, opt 0
protein = ['H', 'C', 'P', 'H', 'P', 'H', 'P', 'H', 'C', 'H', 'H', 'H', 'H', 'P', 'C', 'C', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'C', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 'C', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H'] # 50, opt 0
# protein = ['H', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'P', 'P', 'H', 'P', 'H', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'P', 'H', 'H'] # 50, opt 0
length = len(protein)

# set theoretical lower bound on score
even = protein[::2]
odd = protein[1::2]
H_count = protein.count('H')
C_count = protein.count('C')
min_score = 2 * max([- even.count('H') - 5 * even.count('C'), - odd.count('H') - 5 * odd.count('C')])



# list of timer to see why program so slow
possible_list = []
score_list = []
direction_list = []
double_list = []

def main():
    
    bar = Bar('Progress', max=length)
    bar.next()
    bar.next()
    k = 0

    # specifications for depth first tree building
    depth = length - 2
    q = queue.Queue()
    q.put(('',0))
    final_configurations = []

    # keep track of scores per substring
    lowest_score_k = {}
    all_scores_k = {}
    lowest_score = 0
    
    # (0,1) probabilities of pruning a path, lower is more exact but less fast
    p1 = 0.99
    p2 = 0.9

    # set inital values
    for i in range(length + 1):
        lowest_score_k[i] = 0
        all_scores_k[i] = [0]

    # create a breadth first tree
    while not q.empty():

        state = q.get()

        # if all aminos are placed, put the string in a list
        if len(state[0]) == depth and not double(state[0]):
            final_configurations.append(state)
      
        if len(state[0]) < depth:
            
            for i in ['L', 'R', 'S']:

                child = copy.deepcopy(state) 

                temp_list = list(child)
                temp_list[0] += i
                child = tuple(temp_list)

                if double(child[0]):
                    continue
                
                if len(child[0]) + 1 > k:
                    bar.next()
          
                k = len(child[0]) + 1
               
                # P's are always placed, rest have some conditions
                if not protein[k] == 'P':
                                   
                    score = child[1] + score_func(child[0])
                    
                    possible_score = possible_score_func(protein[k+1:], score)

                    if score + possible_score > lowest_score:
                        continue

                    # random number between 0 and 1
                    r = random.random()
                    
                    average_score_k = sum(all_scores_k[k]) / len(all_scores_k[k])

                    # conditions for pruning
                    if score > average_score_k and r < p1:
                        continue
                    elif (average_score_k >= score and score > lowest_score_k[k]) and r < p2:
                        continue
                    
                    temp_list = list(child)
                    temp_list[1] = score
                    child = tuple(temp_list)  

                    # add to tree
                    q.put(child) 

                    all_scores_k[k].append(score)

                    if score < lowest_score_k[k]:
                        lowest_score_k[k] = copy.deepcopy(score)
                    
                    if score < lowest_score:
                        lowest_score = copy.deepcopy(score)
                        
                else:
                    q.put(child) 

    if len(final_configurations) == 0:
        bar.finish()
        print('No conformations found.')
        return

    lowest_score = 0
    
    # weed out the best configuration from the remaining strings
    for c in final_configurations:
        if c[1] < lowest_score:
            best_config = copy.deepcopy(c[0])
            lowest_score = copy.deepcopy(c[1])

    # plot the result
    stop = timeit.default_timer()
    # print(f'Strings made: {len(final_configurations)}')
    total_time = stop - start
    bar.finish()
    print(f'Length: {length}')
    print(f'Score: {lowest_score}')
    print(f'Total runtime: {total_time}')
    print(f'Possible score function time: {round(sum(possible_list) / total_time * 100,1)}%')
    print(f'Score function time: {round(sum(score_list) / total_time * 100,1)}%')
    print(f'Direction function time: {round(sum(direction_list) / total_time * 100,1)}%')
    print(f'Double function time: {round(sum(double_list) / total_time * 100,1)}%')
    print(f'Conformation: {best_config}')
    plot(best_config, lowest_score,total_time, possible_list, score_list, direction_list, double_list)

def possible_score_func(nodes_to_visit, partial_score):
    """Calculates the best score the remaining bit of the protien can acquire."""
    possible_start = timeit.default_timer()
    possible_score = 0

    # an H atom can get at most -2 and a C atom at best -10
    for i, n in enumerate(nodes_to_visit):

        if i == len(nodes_to_visit) - 1:
            if n == 'H':
                possible_score += -3
            if n == 'C':
                possible_score += -15
            continue

        if n == 'H':
             possible_score += -2
        if n == 'C':
            possible_score += -10

    if possible_score + partial_score < min_score:
        possible_score = min_score - partial_score
    
    possible_stop = timeit.default_timer()
    possible_list.append(possible_stop-possible_start)
    return possible_score

def plot(string, score, runtime, possible, _score, direction, double):
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
    fig, (ax1) = plt.subplots(1, 1)
  
    ax1.plot(list_x, list_y, '--', color='darkgrey')
    ax1.plot(red_dots_x, red_dots_y, 'or', markersize=17)
    ax1.plot(blue_dots_x, blue_dots_y, 'ob', markersize=17)
    ax1.plot(yellow_dots_x, yellow_dots_y, 'oy', markersize=17)
    ax1.axis('equal')
    ax1.set_title(f'Folded protein of length {length}, score: {score}, total: {runtime}')

    # ax2.plot(possible)
    # ax2.set_title(f'possible score func, total: {sum(possible)}')

    # ax3.plot(_score)
    # ax3.set_title(f'score func, total: {sum(_score)}')

    # ax4.plot(direction)
    # ax4.set_title(f'direction to xy func, total: {sum(direction)}')

    # ax5.plot(double)
    # ax5.set_title(f'double func, total: {sum(double)}')

    plt.show()
    

def score_func(string):
    """Given the coordinates of a protein string, calculate the score of the shape."""
    score_start = timeit.default_timer()

    pos_x, pos_y = direction_to_xy(string)

    score = 0
    coordinates = []
    check_coordinates = []

    for x, y in zip(pos_x, pos_y):
        coordinates.append([x,y])

    delta_x = pos_x[-1] - pos_x[-2]
    delta_y = pos_y[-1] - pos_y[-2]

    left_x = pos_x[-1] - delta_y
    left_y = pos_y[-1] + delta_x
    check_coordinates.append([left_x, left_y])

    right_x = pos_x[-1] + delta_y
    right_y = pos_y[-1] - delta_x
    check_coordinates.append([right_x, right_y])

    straight_x = pos_x[-1] + delta_x
    straight_y = pos_y[-1] + delta_y
    check_coordinates.append([straight_x, straight_y])
    
    for c in check_coordinates:

        if c in coordinates:

            index = coordinates.index(c)

            if protein[len(pos_x) - 1] == 'H' and (protein[index] == 'H' or protein[index] == 'C'):
                score += -1
            elif protein[len(pos_x) - 1] == 'C' and protein[index] == 'H':
                score += -1
            elif protein[len(pos_x) - 1] == 'C' and protein[index] == 'C':
                score += -5

    score_stop = timeit.default_timer()
    score_list.append(score_stop-score_start)

    return score

# def score_func(string):
#     """Given the coordinates of a protein string, calculate the score of the shape."""
#     list_x, list_y = direction_to_xy(string)
#     # list to place the 'already scored' atoms into
#     coordinates = []
#     directions = [[-1,0],[0,1],[1,0],[0,-1]]
#     score = 0

#     for i in range(length):

#         # P's dont interact so can skip those cases
#         if not protein[i] == 'P':

#             # for every atom look around in all 4 directions
#             for d in directions:

#                 # check whether one of the previously placed atoms is in the vicinity and determine the score of the interaction with it and the current atom
#                 for j in range(len(coordinates)):

#                     if [list_x[i] + d[0], list_y[i] + d[1]] == coordinates[j] and not(list_x[i] + d[0] == list_x[i-1] and list_y[i] + d[1] == list_y[i-1]):
                        
#                         if protein[i] == 'H':
#                             if protein[j] == 'H' or protein[j] == 'C':
#                                 score += -1

#                         if protein[i] == 'C':
#                             if protein[j] == 'C':
#                                 score += -5
#                             if protein[j] == 'H':
#                                 score += -1 

#         # place in the list with coordinates
#         coordinates.append([list_x[i], list_y[i]])
    
#     return score

def direction_to_xy(string):
    """Converts a series of string with directions like ['L', 'R'] to lists with xy positions."""
    direction_start = timeit.default_timer()
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
           
    direction_stop = timeit.default_timer()
    direction_list.append(direction_stop-direction_start)
    return pos_x, pos_y


def double(string):

    """Checks whether two atoms occupy the same point."""
    double_start = timeit.default_timer()
    list_x, list_y = direction_to_xy(string)
    coordinates = []
    for x, y in zip(list_x, list_y):
        coordinates.append((x,y))
    sett = set(coordinates)
    if len(coordinates) == len(set(coordinates)):
        double_stop = timeit.default_timer()
        double_list.append(double_stop-double_start)
        return False
    
    double_stop = timeit.default_timer()
    double_list.append(double_stop-double_start)
    return True


if __name__ == "__main__":
    start = timeit.default_timer()
    main()