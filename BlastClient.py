# Zachary Tokuno (ztokuno@ucsc.edu)

import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

class BlastClient:
    '''Handle BLAST IO by printing out best hsp (high scoring pair).'''
    def __init__(self):
        '''Instantiate class with attributes for BLAST query.'''
        self.blastType = 'blastp'
        self.db = 'pdb'
        self.threshold = 1e-6

        if sys.argv[1].endswith('.fasta') or sys.argv[1].endswith('.fa') or sys.argv[1].endswith('.fna'):
            self.fa = SeqIO.read(sys.argv[1], format = 'fasta')
        else:
            self.query = sys.argv[1]

    def getQuery(self):
        '''Return query sequence.'''
        if self.fa:
            return self.fa.seq
        else:
            return self.query

    def doBlast(self):
        '''Return BLAST results as a handle with multiple records in handle.'''
        if self.fa:
            blastResults = NCBIWWW.qblast(self.blastType, self.db, self.fa)
        else:
            blastResults = NCBIWWW.qblast(self.blastType, self.db, self.query)
        return blastResults

    def getHsps(self, blastResults):
        '''Parse through hsp in result handle.'''
        blastRecords = NCBIXML.parse(blastResults)
        hsps = []
        for record in blastRecords:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < self.threshold:
                        hsps.append(alignment.title)
                        hsps.append(alignment.length)
                        hsps.append(hsp.expect)
                        hsps.append(hsp.query)
                        hsps.append(hsp.match)
                        hsps.append(hsp.sbjct)
        return hsps

# blaster = BlastClient()
# blastResults = blaster.doBlast()
# hsps = blaster.getHsps(blastResults)
# print(hsps)

def main():
    # Something important to put here
    return

if __name__ == '__main__':
    main()
