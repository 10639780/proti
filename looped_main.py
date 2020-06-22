"""
main.py 

Minor Programmeren
Team Proti

All the algorithms can be accessed.
In the user interface the user can enter an algorithm to run 
and they can choose a protein to fold by entering the ID of the protein.
The available proteins and their respective ID's are printed to the terminal 
for the user to see
"""

# import modules
from algorithms import BFBAB as BF
from algorithms import BFBAB3D as BF3D
from algorithms import DEE
from algorithms import FF
from algorithms import HC
from algorithms import MONTE
from algorithms import MONTE3D as M3D
from algorithms import SIMANN as SA
from algorithms import SIMANN3D as SA3D
from algorithms import SIMANNPLUS as SAPLUS
from algorithms import SIMANNPLUS3D as SAPLUS3D
from algorithms import TREE 
from algorithms import GENETIC
from algorithms import HC3D
from classes.protein import Protein
import sys

def main():

    # available algorithms
    algos = {'BFBAB':BF, 'BFBAB3D':BF3D, 'DEE':DEE, 'FF':FF, 'HC':HC, \
            'HC3D':HC3D, 'MONTE':MONTE, 'MONTE3D':M3D, 'SIMANN':SA, \
            'SIMANN3D':SA3D,'SIMANN+':SAPLUS, 'SIMANN3D+':SAPLUS3D, \
            'TREE':TREE, 'GENETIC':GENETIC} 

    proteins = ['HHPHHHPH','HHPHHHPHPHHHPH','HPHPPHHPHPPHPHHPPHPH','PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP','HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH','PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP','CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC','HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH','HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH']
    
    r = open('SIMANN+_results.txt', 'w')

    for line in proteins:
        line.strip("\n")
        time_list = []
        score_list = []
        conform_list = []

        fold_protein = line
        length = len(fold_protein)
        proti = Protein(fold_protein, length)
        print(f'line: {line}')
        algo = 'SIMANN+'
        
        for i in range(30):
            # run selected algorithm
            try:
                time, score, conform = algos[algo].run(proti, 33 * i)
            except:
                time, score, conform = 0, 0, 'not found'                                                                                                                                                                                                                                                     
            time_list.append(time)
            score_list.append(score)
            conform_list.append(conform)
        
        r.write(f'sequence = {line}\ntime = {time_list}\nscore = {score_list}\nconform = {conform_list}\n')

    r.close()

if __name__ == "__main__":
    main()
