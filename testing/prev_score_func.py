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