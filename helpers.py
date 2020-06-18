import matplotlib.pyplot as plt

def possible_score_func(nodes_to_visit, partial_score, min_score):
    """Calculates the best score the remaining bit of the protien can acquire."""

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

    return possible_score


def plot(string, score, protein):
    """Makes a graph of two lists list_x, list_y."""

    list_x, list_y = direction_to_xy(string)

    # differentiate between types of atom
    red_dots_x = []
    red_dots_y = []
    blue_dots_x = []
    blue_dots_y = []
    yellow_dots_x = []
    yellow_dots_y = []

    # search through protein and place each atom in the appropriate list
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
    fig, ax1 = plt.subplots(1, 1, figsize=(6, 6))

    ax1.plot(list_x, list_y, '--', color='darkgrey')
    ax1.plot(red_dots_x, red_dots_y, 'or', markersize=17)
    ax1.plot(blue_dots_x, blue_dots_y, 'ob', markersize=17)
    ax1.plot(yellow_dots_x, yellow_dots_y, 'oy', markersize=17)
    ax1.set_title(f'Folded protein of length {len(protein)}, score: {score}')

    plt.show()


def score_func(string, protein):
    """Given the coordinates of a protein string, calculate the score of the shape."""
    # list to place the 'already scored' atoms into
    coordinates = []
    directions = [[-1, 0], [0, 1], [1, 0], [0, -1]]
    score = 0
    list_x, list_y = direction_to_xy(string)

    length = len(list_x)

    for i in range(length):

        # P's dont interact, so skip those cases
        if not protein[i] == 'P':

            # for every element look around in all 4 directions
            for d in directions:
                # check whether one of the previously placed atoms is in the vicinity and determine the score of the interaction with it and the current atom
                for j in range(len(coordinates)):

                    if [list_x[i] + d[0], list_y[i] + d[1]] == coordinates[j] and not (
                            list_x[i] + d[0] == list_x[i - 1] and list_y[i] + d[1] == list_y[i - 1]):

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


def direction_to_xy(string):
    """Converts a series of string with directions like ['L', 'R'] to lists with xy positions."""

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


def double(string):
    """Checks whether two atoms occupy the same point."""

    list_x, list_y = direction_to_xy(string)

    coordinates = []

    # see if a coordinate is already in the list, then add that coordinate to the list
    for x, y in zip(list_x, list_y):
        if [x, y] in coordinates:
            return True
        coordinates.append([x, y])

    return False
