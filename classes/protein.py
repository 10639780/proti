class Protein:
    def __init__(self, listed, length):
        self.listed = listed
        self.length = length
        self.even = listed[::2]
        self.odd = listed[1::2]
        self.min_score = 2 * max([- self.even.count('H') - 5 * self.even.count('C'), - self.odd.count('H') - 5 * self.odd.count('C')])