"""
From this main landing, multiple routes can be chosen. One at a time.
"""

from algorithms import *
import re
import copy
import timeit
import random
import queue
import sys

def main():
    # Input for protein variable
    algos = 'BFBAB', 'BFBAB3D', 'DEE', 'FF', 'HC', 'MONTE', 'MONTE3D', 'SIMANN', 'SIMANN+', 'TREE'
    print(f'algorithms available: {algos}')
    usage = input("please input algorithm name followed by an underscore and a single protein string")
    usage = re.split('_+', usage)
    algo = usage[0]
    protein = list(usage[1])
    allow_chars = ['P','C','H']

    for i in protein:
        if i not in allow_chars:
            sys.exit("Usage parameters not met. Please try again")

    if algo not in algos:
        sys.exit("Usage parameters not met. Please try again")

    algo.run()

