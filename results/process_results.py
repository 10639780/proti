import statistics

def output(time_list, score_list):

    time_list, score_list = remove_invalid(time_list, score_list)


    time_avg, score_avg = statistics.mean(time_list), statistics.mean(score_list)
    time_stdev, score_stdev = statistics.stdev(time_list), statistics.stdev(score_list)

    print(f'{round(score_avg,1)}(\xb1{round(score_stdev,1)}),  {round(time_avg,1)}(\xb1{round(time_stdev,1)})')

    return  

def remove_invalid(time_list, score_list):

    to_remove = []

    for i in range(len(time_list)):

        if score_list[i] == 0:
            to_remove.append(i)
         
    to_remove.reverse()

    for i in to_remove:
        del time_list[i]
        del score_list[i]

    return time_list, score_list

time = [33.0756315, 33.283800500000005, 33.8387254, 33.203660799999994, 33.483201399999984, 33.56407949999999, 33.137783300000024, 33.19153020000002, 33.18942139999996, 33.758493199999975, 33.1253815, 33.61057390000002, 33.90959779999997, 33.28172219999999, 33.735907499999996, 33.539113400000076, 33.6643722, 34.301181999999926, 32.859732099999974, 33.37528540000005, 33.181631100000004, 32.94277390000002, 33.502696300000025, 32.81257040000003, 33.19721479999998, 33.494911900000034, 34.04892399999994, 33.724170999999956, 33.616994699999964, 35.1726936]
score = [-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6]



output(time, score)

# output: avgscre(sted), avgtime(stdev)