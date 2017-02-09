# Zachary Tokuno (ztokuno@ucsc.edu)

import sys
from prody import *

class PfamClient:
    '''Handle pfam IO by searching for conserved domains with PDB ID or query sequence.'''
    def __init__(self):
        return

    def doPfamSearch(self, pdb, query):
        '''Search pfam using PDB ID or query sequence. Must give pdb and query sequence.'''
        try:
            pfamResults = searchPfam(pdb)
        except IndexError:
            sys.stderr.write('*****************************************************\n')
            sys.stderr.write('PDB ID look up failed. Searching with query sequence.\n')
            sys.stderr.write('*****************************************************\n')
            pfamResults = searchPfam(query)
        return pfamResults

    def getConservedDomains(self, pfamResults):
        '''Parse out conserved domains found in the PDB or query sequence (homolog).'''
        # python3 returns a dict_keys object, so cast into a list
        accession = list(pfamResults.keys())[0]
        conservedDomains = []
        for alignment in pfamResults[accession]['locations']:
            start = int(alignment['ali_start'])
            end = int(alignment['ali_end'])
            conservedDomains.append((start, end))
        return conservedDomains

# pfamer = PfamClient()
# pfamResults = pfamer.doPfamSearch('4xmm', 'MKRESHKHAEQARRNRLAVALHELASLIPAEWKQQNVSAAPSKATTVEAACRYIRHLQQNGST')
# conservedDomains = pfamer.getConservedDomains(pfamResults)

def main():
    # Something important to put here
    return

if __name__ == '__main__':
    main()
