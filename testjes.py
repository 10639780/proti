"""
File om testjes te runnen voordat we aanpassingen doen in prot.py
"""

"""
Deze functie zouden we kunnen grebruiken om de proteins uit een file te lezen 
Waardoor we meerdere proteins in 1 keer door de code kunnen sturen (van prot.py)
Ik ga er later mee door om het in prot.py te verwerken
"""
def get_proteins():

    """Saves the proteins to analyse in a dictionary"""

    proteins = {}
    filename = "proteins.txt"

    with open(filename, "r") as file:
        lines = file.readlines()

        for i in range(len(lines)):
            protein = lines[i].strip("\n") 
            length = len(protein)
            proteins[i] = [protein, length] 

    return proteins 

def make_grid(proteins):
    """
    Creates a grid of specified size.
    Saves the grid for all proteins in the dictionary.
    """

    for key in proteins:

        grid = []
        size = 2*proteins[key][1]+1
        for i in range(size):
            row = []
            for j in range(size):
                row.append(0)
            grid.append(row)

        proteins[key].append(grid)

    return proteins

        
if __name__ == "__main__":
    p = get_proteins()
    g = make_grid(p)