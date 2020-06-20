"""
TREE.py
Minor Programmeren
Team Proti
Determines the highest possible score a protein fold can get.
"""
from anytree import Node, RenderTree, NodeMixin, LevelOrderGroupIter, PreOrderIter
import copy
import matplotlib.pyplot as plt 

class Atom(NodeMixin):
    """Creates nodes for the tree."""

    def __init__(self, name, direction=None, parent=None, children=None):
        """Every node has the atom type, and direction relative to the previous atom."""

        super(Atom, self).__init__()
        self.name = name
        self.direction = direction
        self.parent = parent
        if children:
            self.children = children

# all directions an atom can fold into
directions = [[1,0],[-1,0],[0,1],[0,-1]]

# strings of different  length, works well with string up to 11 atoms long, after that the program becomes very slow. 
# protein = ['P', 'H', 'P', 'H']
protein = ['H','H','P','H','H','H','P','H']
# protein = ['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H', 'P', 'H','H','H','P','H']
# protein = ['H', 'H', 'H', 'C', 'C', 'C', 'H', 'C', 'H']

length = len(protein)

def main():

    # create a tree with all possible folds
    root = create_tree()
    # display_tree(root)

    # extract all possible protein configurations from the tree
    route_directions = create_routes(root)

    # for r in route_directions:
    #     print(r)
    
    print(f'Total possible structures: {len(route_directions)}')

    # variables to hold the optimal solution
    best_x = []
    best_y = []
    lowest_score = 0
    same_score_counter = 0
    invalid_solutions = 0

    # loop over each of the possible configurations
    for r in route_directions:

        pos_x = []
        pos_y = []
        x = 0
        y = 0

        # lists all filled with coordinates of the atoms
        for d in r:

            x += d[0]
            y += d[1]

            pos_x.append(x)
            pos_y.append(y)

        # discard all the structures that fold into themselves
        if double(pos_x, pos_y):
            invalid_solutions += 1
            continue
        
        # calculate te score of the current structure
        current_score = score(pos_x, pos_y)

        # save the structure with the lowest score
        if current_score < lowest_score:
            best_x = copy.deepcopy(pos_x)
            best_y = copy.deepcopy(pos_y)
            lowest_score = copy.deepcopy(current_score)
            same_score_counter = 0
        if current_score == lowest_score:
            same_score_counter += 1
    
    print(f'Invalid solutions: {invalid_solutions}')
    print(f'Same score solutions: {same_score_counter}')
    print(f'Lowest score: {lowest_score}')

    # print best structure to csv file
    output(best_x, best_y, lowest_score)

    # plot the structure with the best score
    plot(best_x, best_y, lowest_score)


def output(list_x, list_y, score):
    """Prints the folded string to a csv file in the Bas Terwijn style."""

    numbers = []
   
    for i in range(len(protein)-1):
        
        # new position is compared to the old
        delta_x = list_x[i+1] - list_x[i]
        delta_y = list_y[i+1] - list_y[i]

        # conversion between coordinates  and bas terwijn numbers
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

    # write the list to a file
    f = open('output.csv', 'w')
    f.write('amino,fold\n')
    for p, n in zip(protein, numbers):
        f.write(f'{p}, {n}\n')
    f.write(f'score,{score}') 
    f.close()


def plot(list_x, list_y, score):
    """Makes a graph of two lists list_x, list_y."""
    
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

    plt.plot(list_x, list_y, '--', color='darkgrey')
    plt.plot(red_dots_x, red_dots_y, 'or')
    plt.plot(blue_dots_x, blue_dots_y, 'ob')
    plt.plot(yellow_dots_x, yellow_dots_y, 'oy')
    plt.title(f'Folded protein, score: {score}')

    plt.show()


def score(list_x, list_y):
    """Given the coordinates of a protein string, calculate the score of the shape."""

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


def double(list_x, list_y):
    """Checks whether two atoms occupy the same point."""

    coordinates = []

    # see if a coordinate is already in the list, then add that coordinate to the list
    for x, y in zip(list_x, list_y):
        if [x,y] in coordinates:
            return True
        coordinates.append([x,y])
    
    return False
 

def create_routes(root):
    """Creates a list that hold all the possible configurations a protein could have."""

    route_directions = []
    
    # list of all possible paths in the tree
    paths = [list(leaf.path) for leaf in PreOrderIter(root, filter_=lambda node: node.is_leaf)] # copied from https://stackoverflow.com/questions/59917058/how-to-get-all-possible-branch-with-python-anytree

    # add those paths to a list
    for p in paths:
        direction = []

        for node in p:
            direction.append(node.direction)

        route_directions.append(direction)
    
    return route_directions


def create_tree():
    """Create a tree structure, each atom can branch off into 4 directions, creating 4^n possible routes."""

    node_counter = 0
    parent_counter = 0
    node_list = []
    first = True
    second = True

    for p, i in zip(protein, [i for i in range(length)]):

        # first atom has no parents and is located at 0,0
        if first:
            root = Atom(p, [0,0])
            node_list.append(root)
            first = False
            continue
        
        # second atom             # first atom has no parents and is located at 0,0            # first atom has no parents and is located at 0,0            # first atom has no parents and is located at 0,0
        if second:
            node =  Atom(p, [1,0], parent=root)
            node_list.append(node)
            second = False
            continue

        for k in range((len(directions) - 1) ** (i - 2)):

            # make sure the protein doesn't walk back into itself
            previous_direction = copy.deepcopy(node_list[parent_counter + 1].direction)
            previous_direction[0] = -previous_direction[0]
            previous_direction[1] = -previous_direction[1]

            temp_direction = copy.deepcopy(directions)

            if previous_direction in temp_direction: 
                temp_direction.remove(previous_direction)

            # fill the tree and keep track of the parent child relations
            for j in range(len(directions) - 1):
            
                j = Atom(p, temp_direction[j], parent=node_list[parent_counter + 1] )
                node_list.append(j)

                node_counter += 1
                if node_counter % (len(directions) - 1) == 0:
                    parent_counter += 1

        # print statement for the long calculations
        print(f'Loading tree: {(i + 1)/length * 100}%')

    return root


def display_tree(root):
    """Print a graphical representation of the tree to the terminal."""

    for pre, fill, node in RenderTree(root):
        treestr = u"%s%s" % (pre, node.name)
        print(treestr.ljust(20), node.direction)


if __name__ == "__main__":
    main()