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
            proteins[i] = lines[i].strip("\n")

    return proteins 

def length_protein(proteins):
    for key in proteins:
        print(proteins[key])
        print(len(proteins[key]))


        
if __name__ == "__main__":
    p = get_proteins()
    length_protein(p)