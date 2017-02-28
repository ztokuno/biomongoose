#!/usr/bin/env python
# Zachary Tokuno (ztokuno@ucsc.edu)

import random

class kmerGenerator:
    def __init__(self, k):
        self.k = k
        self.dna = 'ATCG'
        self.rna = 'AUCG'

    def generateSequences(self):
        '''Generate random DNA sequences of length k.'''
        seq = ''
        for pos in range(self.k):
            seq += self.dna[random.randint(0, 3)]
        print(seq)

def main():
    kmerGen = kmerGenerator(13)
    kmerGen.generateSequences()

if __name__ == '__main__':
    main()
