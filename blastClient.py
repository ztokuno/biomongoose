#!usr/bin/env python3

# Written by: Zachary Tokuno

import argparse
import numpy
import scipy
import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

class CommandLine():
    def __init__(self, forceParams = None):
        '''
        Create ArgumentParser object.
        Add arguments to parser.

        Accepts optional parameter which is a list containing parameters.
        Used to simulate command line within IDLE.
        '''

        # Change descriptions later
        self.parser = argparse.ArgumentParser(description = 'Perform a BLAST search with input fasta file as query.',
                                            add_help = True,
                                            usage = '%(prog)s -e 0.000001 -f in.fasta >output',
                                            prefix_chars = '-'
                                            )
        self.parser.add_argument('-e', '--eValThreshold', action = 'store', type = float, default = 1e-5, help = 'E-value threshold used to filter homologues.')
        self.parser.add_argument('-f', '--inFasta', action = 'store', type = str, default = None, help = 'Input fasta file to use as BLAST query.')

        # simulate command line by passing a list of parameters
        if forceParams is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(forceParams)

class BLASTClient:
    '''
    Handle BLAST IO.
    '''

    def __init__(self, eValThreshold, inFasta = None):
        '''
        Instantiate class with attributes for BLAST query.
        '''

        self.blastType = 'blastp'
        self.db = 'nr'
        self.eValThreshold = eValThreshold

        if inFasta:
            self.fa = SeqIO.read(inFasta, format = 'fasta')
        else:
            self.fa = sys.stdin

    def getBLAST(self):
        '''
        Return BLAST results as a handle with multiple records in handle.
        '''

        blastResults = NCBIWWW.qblast(self.blastType, self.db, self.fa)

        return blastResults

    def parseBLAST(self, blastResults):
        '''
        Parse through records in result handle.
        '''

        blastRecords = NCBIXML.parse(blastResults)

        numHsps = 10

        for record in blastRecords:
            for alignment in record.alignments:
                for hsp in alignment.hsps[:10]:
                    print('OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO')
                    if hsp.expect < self.eValThreshold:
                        print('*****Alignment*****')
                        print('Seq:', alignment.title)
                        print('Len:', alignment.length)
                        print('E val:', hsp.expect)
                        print(hsp.query)
                        print(hsp.match)
                        print(hsp.sbjct)
                        print()


def main(myCommandLine = None):
    if myCommandLine == None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(myCommandLine)

    blaster = BLASTClient(myCommandLine.args.eValThreshold, myCommandLine.args.inFasta)
    blastResults = blaster.getBLAST()
    blaster.parseBLAST(blastResults)






if __name__ == '__main__':
    main()
