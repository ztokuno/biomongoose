# Zachary Tokuno (ztokuno@ucsc.edu)

import sys
from BlastClient import *
from BlastParser import *
from PfamClient import *

def main():
    print('*****************************************************')
    print('BLAST...')
    print('*****************************************************')
    blaster = BlastClient()
    blastResults = blaster.doBlast()
    hsps = blaster.getHsps(blastResults)
    query = blaster.getQuery()

    if hsps == []:
        sys.stderr.write('No significant BLAST results found (e-value greater than 1e-6).\n')
        sys.exit(1)

    print('Parsing...')
    print('*****************************************************')
    parser = BlastParser()
    parser.getHspList(hsps)
    pdb = parser.getPdb()
    e = parser.getEVal()

    print('Pfam...')
    print('*****************************************************')
    pfamer = PfamClient()
    pfamResults = pfamer.doPfamSearch(pdb, query)
    conservedDomains = pfamer.getConservedDomains(pfamResults)

    print('*****************************************************')
    print('Query sequence most closely matched to '+pdb+' with an e-value of '+str(e))
    print('Range of conserved domains in 1-based coordinates:')
    print(conservedDomains)

if __name__ == '__main__':
    main()
