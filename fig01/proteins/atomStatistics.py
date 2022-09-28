# Class AtomStatistics
# v. 0.1
# G. Tiana, 18 May 2020
#

from math import sqrt

class AtomStatistics:
    def __init__(self):
        self.histo = {'H':0, 'C':0, 'O':0, 'N':0, 'P':0, 'S':0 }

    def CompareHistoWith(self, otherStatistics):
        q = 0
        for a in {'H', 'N', 'C', 'O', 'S', 'P'}:
            q = q + ( self.histo[a] - otherStatistics.histo[a] )**2

        return sqrt(q)

    def HistoFromSmiles(self, smiles):
        for a in smiles:
            try:
                self.histo[string.upper(a)] = self.histo[string.upper(a)] + 1
            except:
                continue
