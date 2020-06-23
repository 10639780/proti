# proti
Proteine folding tool

Once the user has cloned the repository and installed the requirements, running 'python3 main.py' will display  the available algorithms. After entering the algorithm of choice the user can select the protein string with which the program will run. Proteinstrings that are longer than twenty elements may take a long time to run. Especially TREE.py and DEE.py will become impractically slow from a string length of 14 and 20 respectively. These algorithms will guarantee an optimal folding configuration at the cost of having to consider an enormous amount of possibilities.

For the Breadth First Branch & Bound algorithm (2D & 3D) we took inspiration from the following source:
Chen, Mao, and Wen-Qi Huang. "A branch and bound algorithm for the protein folding problem in the HP lattice model." Genomics, proteomics & bioinformatics 3.4 (2005): 225-230.

For the DEE algorithm (2D) we would like to thank Okke for getting us on the right track during our techassist meetings. 

For the Firefly algorithm (2D) we took inspiration from: Zhang, Yudong, Lenan Wu, and Shuihua Wang. "Solving two-dimensional HP model by firefly algorithm and simplified energy function." Mathematical Problems in Engineering 2013 (2013).

The Genetic algorithm was based on Bui, Thang N., and Gnanasekaran Sundarraj. "An efficient genetic algorithm for predicting protein tertiary structures in the 2D HP model." Proceedings of the 7th annual conference on Genetic and evolutionary computation. 2005, supplemented by Halm, Eyal. "Genetic Algorithm for Predicting Protein Folding in the 2D HP Model." where we found the crossover funtion.

The idea to start with a full string and rotate in random directions in Monte Carlo was insprired by Ramji T. Venkatasubramanian's Computational Nanomechanics class assignment which somehow made its way to the internet. 

We hope you will enjoy our work, 

Xamanie, Jesse, Johan.
