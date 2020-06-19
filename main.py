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
from algorithms import SIMANNPLUS as SAPLUS
from algorithms import TREE as TREE
from classes.protein import Protein
import re
import sys

if __name__ == "__main__":
    # Input for protein variable

    algos = 'BFBAB', 'BFBAB3D', 'DEE', 'FF', 'HC', 'MONTE', 'MONTE3D', 'SIMANN', 'SIMANN+', 'TREE'
    """
    print(f'algorithms available: {algos}')
    usage = input("please input algorithm name followed by an underscore and a single protein string\n")
    usage = re.split('_+', usage)
    algo = usage[0]
    protein_list = list(usage[1])
    length = len(protein_list)
    """
    allow_chars = ['P','C','H']
    protein_list = ['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'P', 'H']
    #protein_list = ['H', 'P', 'H', 'P', 'P', 'H', 'H', 'P', 'H', 'P', 'P', 'H', 'P', 'H', 'H', 'P', 'P', 'H', 'P', 'H']
    length = len(protein_list)
    proti = Protein(protein_list, length)

    algo = 'MONTE3D'

    for i in protein_list:
        if i not in allow_chars:
            sys.exit("Usage parameters not met. Please try again")

    if algo not in algos:
        sys.exit("Algorithm not found. Please try again")



    if algo == 'BFBAB':
        BF.run(proti)

    if algo == 'BFBAB3D':
        BF3D.run(proti)

    if algo == 'DEE':
        DEE.run(proti)

    if algo =='FF':
        FF.run(proti)

    if algo =='HC':
        HC.run(proti)

    if algo =='MONTE':
        MONTE.run(proti)

    if algo =='MONTE3D':
        M3D.run(proti)

    if algo =='SIMANN':
        SA.run(proti)

    if algo =='SIMANN+':
        SAPLUS.run(proti)

    if algo =='TREE':
        TREE.run(proti)







