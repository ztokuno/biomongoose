# Zachary Tokuno (ztokuno@ucsc.edu)

class BlastParser:
    '''Parse BLAST output from BlastClient.py.'''
    def __init__(self):
        '''Instantiate empty containers for attributes.'''
        self.pdb = None
        self.eVal = None
        self.query = None

    def setPdb(self, pdb):
        '''Set PDB id.'''
        self.pdb = pdb

    def getPdb(self):
        '''Get PDB id.'''
        return self.pdb

    def setEVal(self, eVal):
        '''Set e value.'''
        self.eVal = eVal

    def getEVal(self):
        '''Get e value.'''
        return self.eVal

    def setQuery(self, query):
        '''Set query (input amino acid) sequence.'''
        self.query = query

    def getQuery(self):
        '''Get query sequence.'''
        return self.query

    def getHspList(self, hsps):
        '''Read all lines from BlastClient into memory. Split on alignment field. Take data only from the first hsp.'''
        pdb = hsps[0].split('pdb|')[1][:4]
        eVal = hsps[2]
        query = hsps[3].strip()
        self.setPdb(pdb)
        self.setEVal(eVal)
        self.setQuery(query)

    def __str__(self):
        '''Custom __str__ method to print important attributes of class.'''
        toString = 'PDB ID: '+self.pdb+'\n'  +self.eVal+'\n'+'Query: '+self.query
        return toString

#parser = BlastParser()
#parser.getHspList()
#print(parser)

def main():
    # Something important to put here
    return

if __name__ == '__main__':
    main()
