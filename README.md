# proti
Proteine folding tool

Once the user has cloned the repository and installed the requirements, running 'python3 main.py' will display  the available algorithms. After entering the algorithm of choice the user can select the protein string with which the program will run. Protein strings that are longer than twenty elements may take a long time to run. 

-Exceptions listed here-

For the Breadth First Branch & Bound algorithm (2D & 3D) we took inspiration from the following source:
Chen, Mao, and Wen-Qi Huang. "A branch and bound algorithm for the protein folding problem in the HP lattice model." Genomics, proteomics & bioinformatics 3.4 (2005): 225-230.

For the DEE algorithm (2D) we would like to thank Okke for getting us on the right track during our techassist meetings. 

For the Firefly algorithm (2D) we took inspiration from: Zhang, Yudong, Lenan Wu, and Shuihua Wang. "Solving two-dimensional HP model by firefly algorithm and simplified energy function." Mathematical Problems in Engineering 2013 (2013).

The Genetic algorithm was based on 'An Efficient Genetic Algorithm for Predicting Protein Tertiary Structures in the 2D HP Model' (2005) by Thang N. Bui and Gnanasekaran Sundarraj, supplemented by 'Genetic Algorithm for Predicting Protein Folding in the 2D HP Model' (2007) by Eyal Halm for the section on the crossover funtion.

Monte Carlo was inspired by Ramji T. Venkatasubramanian's Computational Nanomechanics class assignment which somehow made it to the internet. 
