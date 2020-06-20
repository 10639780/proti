"""
From this main landing, multiple routes can be chosen. One at a time.
"""

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
import re
import sys

if __name__ == "__main__":
    
    # read available proteins from text file
    proteins = {}
    filename = "proteins.txt"

    with open(filename, "r") as file:
        lines = file.readlines()

        print("Available proteins: ")
        # save protein string in a dictionary to choose from
        for i in range(len(lines)):
            protein = lines[i].strip("\n")             
            proteins[i] = protein
            print(f"ID: {i}, Protein: {protein}")
    
    # available algorithms
    algos = {'BFBAB':BF, 'BFBAB3D':BF3D, 'DEE':DEE, 'FF':FF, 'HC':HC, \
            'HC3D':HC3D, 'MONTE':MONTE, 'MONTE3D':M3D, 'SIMANN':SA, \
            'SIMANN3D':SA3D,'SIMANN+':SAPLUS, 'SIMANN3D+':SAPLUS3D, \
            'TREE':TREE, 'GENETIC':GENETIC}   

    # get algorithm for folding
    algo = input("\n" + "Enter algorithm to run: ")

    # get protein to fold
    choose_protein = input("Which protein would you like to fold? Enter ID: ")

    # protein information
    fold_protein = proteins[int(choose_protein)]
    length = len(fold_protein)
    proti = Protein(fold_protein, length)

    # check if algorithm name exists
    if algo not in algos:
        sys.exit("Algorithm not found. Please try again")

    # run choosen algorithm
    algos[algo].run(proti)







