# Zachary Tokuno (ztokuno@ucsc.edu)

import sys
from BlastClient import *
from BlastParser import *
from PfamClient import *

def main():
    print('BLAST...')
    blaster = BlastClient()
    blastResults = blaster.doBlast()
    hsps = blaster.getHsps(blastResults)
    query = blaster.getQuery()

    if hsps == []:
        sys.stderr.write('No significant BLAST results found (e-value greater than 1e-6).')
        sys.exit(1)

    print('Parsing...')
    parser = BlastParser()
    parser.getHspList(hsps)
    pdb = parser.getPdb()

    print('Pfam...')
    pfamer = PfamClient()
    pfamResults = pfamer.doPfamSearch(pdb, query)
    conservedDomains = pfamer.getConservedDomains(pfamResults)
    print(conservedDomains)

if __name__ == '__main__':
    main()
