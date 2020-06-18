import queue
import copy

protein = ['H', 'H', 'P', 'H']
length = len(protein)
depth = length - 2


def main():

    q = queue.Queue()
    q.put(('',0))

    while not q.empty():
        state = q.get()
        print(state)

        if len(state[0]) < depth:
            for i in ['L', 'R', 'S']:

                child = copy.deepcopy(state)

                temp_list = list(child)
                temp_list[0] += i
                child = tuple(temp_list)

                score = child[1] + score_func(child[0])

                temp_list = list(child)
                temp_list[1] = score
                child = tuple(temp_list)                

                q.put((child))


def score_func(string):
    
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

            if protein[index] == 'H':
                score += -1

    return score


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

if __name__ == "__main__":
    main()