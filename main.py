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
from algorithms import GENETIC
from classes.protein import Protein
import re
import sys

if __name__ == "__main__":
    # Input for protein variable

    algos = 'BFBAB', 'BFBAB3D', 'DEE', 'FF', 'HC', 'MONTE', 'MONTE3D', 'SIMANN', 'SIMANN+', 'TREE', 'GENETIC'
    """
    print(f'algorithms available: {algos}')
    usage = input("please input algorithm name followed by an underscore and a single protein string\n")
    usage = re.split('_+', usage)
    algo = usage[0]
    protein_list = list(usage[1])
    length = len(protein_list)
    """
    allow_chars = ['P','C','H']

    protein_list = 'HHPHHHPHPHHHPH' # 14, opt -6

    # protein_list = 'HPHPPHHPHPPHPHHPPHPH' # 20, otp -9

    # protein_list = 'PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP' # 36, opt -14

    # protein_list = 'HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH' # 50, opt -21

    # protein_list = 'PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP' #36

    # protein_list = 'CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC' # 36

    # protein_list = 'HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH' # 50

    # protein_list = 'HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH' # 50

    length = len(protein_list)
    proti = Protein(protein_list, length)

    algo = 'GENETIC'

    for i in protein_list:
        if i not in allow_chars:
            sys.exit("Usage parameters not met. Please try again")

    if algo not in algos:
        sys.exit("Algorithm not found. Please try again")


    if algo == 'BFBAB':
        p1 = 0.99
        p2 = 0.98
        BF.run(proti, p1, p2)

    if algo == 'BFBAB3D':
        BF3D.run(proti)

    if algo == 'DEE':
        lowest_known_score = 0
        DEE.run(proti, lowest_known_score)

    if algo =='FF':
        population = 100
        N = 20
        FF.run(proti, population, N)

    if algo =='HC':
        HC.run(proti)

    if algo =='MONTE':
        N = 100
        MONTE.run(proti, N)

    if algo =='MONTE3D':
        M3D.run(proti)

    if algo =='SIMANN':
        SA.run(proti)

    if algo =='SIMANN+':
        SAPLUS.run(proti)

    if algo =='TREE':
        TREE.run(proti)

    if algo == 'GENETIC':
        population = 100
        N = 1000
        GENETIC.run(proti, population, N)







