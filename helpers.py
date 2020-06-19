import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from anytree import Node, RenderTree, Walker, PreOrderIter

def possible_score_func(nodes_to_visit, partial_score, min_score):
    """Calculates the best score the remaining bit of the protein can acquire."""
    #possible_start = timeit.default_timer()
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

    #possible_stop = timeit.default_timer()
    #possible_list.append(possible_stop - possible_start)
    return possible_score


def xyz_plot(string, score, proti):
    # Makes a graph of three lists list_x, list_y and list_z.

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

    # search through protein and place each atom in the appropriate list
    for x, y, z, p in zip(list_x, list_y, list_z, proti.listed):

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
    fig.savefig('images/' + '{}.png'.format(string))
    plt.show()

def plot(string, score, runtime, proti):
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
    fig, (ax1) = plt.subplots(1, 1)

    ax1.plot(list_x, list_y, '--', color='darkgrey')
    ax1.plot(red_dots_x, red_dots_y, 'or', markersize=17)
    ax1.plot(blue_dots_x, blue_dots_y, 'ob', markersize=17)
    ax1.plot(yellow_dots_x, yellow_dots_y, 'oy', markersize=17)
    ax1.axis('equal')
    ax1.set_title(f'Folded protein of length {proti.length}, score: {score}, total: {runtime}')

    # ax2.plot(possible)
    # ax2.set_title(f'possible score func, total: {sum(possible)}')

    # ax3.plot(_score)
    # ax3.set_title(f'score func, total: {sum(_score)}')

    # ax4.plot(direction)
    # ax4.set_title(f'direction to xy func, total: {sum(direction)}')

    # ax5.plot(double)
    # ax5.set_title(f'double func, total: {sum(double)}')

    plt.show()


def score_func(string, proti):
    """Given the coordinates of a protein string, calculate the score of the shape."""
    #score_start = timeit.default_timer()

    pos_x, pos_y = direction_to_xy(string)

    score = 0
    coordinates = []
    check_coordinates = []

    for x, y in zip(pos_x, pos_y):
        coordinates.append([x, y])

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

            if proti.listed[len(pos_x) - 1] == 'H' and (proti.listed[index] == 'H' or proti.listed[index] == 'C'):
                score += -1
            elif proti.listed[len(pos_x) - 1] == 'C' and proti.listed[index] == 'H':
                score += -1
            elif proti.listed[len(pos_x) - 1] == 'C' and proti.listed[index] == 'C':
                score += -5

    #score_stop = timeit.default_timer()
    #score_list.append(score_stop - score_start)
    return score


def xyz_score_func(string, proti):
    """Given the coordinates of a protein string, calculate the score of the shape."""

    # list to place the 'already scored' atoms into
    coordinates = []
    directions = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
    score = 0
    list_x, list_y, list_z = direction_to_xyz(string)
    proti.length = len(list_x)

    for i in range(proti.length):

        # P's dont interact so can skip those cases
        if not proti.listed[i] == 'P':

            # for every atom look around in all 6 directions
            for d in directions:

                # check whether one of the previously placed atoms is in the vicinity and determine the score of the interaction with it and the current atom
                for j in range(len(coordinates)):

                    if [list_x[i] + d[0], list_y[i] + d[1], list_z[i] + d[2]] == coordinates[j] and not (
                            list_x[i] + d[0] == list_x[i - 1] and list_y[i] + d[1] == list_y[i - 1] and list_z[i] + d[
                        2] == list_z[i - 1]):

                        if proti.listed[i] == 'H':
                            if proti.listed[j] == 'H' or proti.listed[j] == 'C':
                                score += -1

                        if proti.listed[i] == 'C':
                            if proti.listed[j] == 'C':
                                score += -5
                            if proti.listed[j] == 'H':
                                score += -1

                                # place in the list with coordinates
        coordinates.append([list_x[i], list_y[i], list_z[i]])

    return score


def direction_to_xy(string):
    """Converts a series of string with directions like ['L', 'R'] to lists with xy positions."""
    direction_list = []
    # first two aminos already placed
    pos_x = [0, 1]
    pos_y = [0, 0]

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
            pos_x.append(pos_x[-1] + delta_y)
            pos_y.append(pos_y[-1] - delta_x)

    return pos_x, pos_y


def direction_to_xyz(string):
    """Converts a series of string with directions like ['L', 'R', 'U', 'D'] to lists with xyz positions."""

    pos_x = [0, 1]
    pos_y = [0, 0]
    pos_z = [0, 0]
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
            pos_x.append(pos_x[-1] + delta_y)
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
    list_x, list_y = direction_to_xy(string)
    coordinates = []
    for x, y in zip(list_x, list_y):
        coordinates.append((x, y))
    if len(coordinates) == len(set(coordinates)):
        return False
    return True

def xyz_double(string):
    """Checks whether two atoms occupy the same point."""

    list_x, list_y, list_z = direction_to_xyz(string)

    coordinates = []

    # see if a coordinate is already in the list, then add that coordinate to the list
    for x, y, z in zip(list_x, list_y, list_z):
        if [x, y, z] in coordinates:
            return True
        coordinates.append([x, y, z])

    return False
