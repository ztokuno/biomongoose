# Zachary Tokuno (ztokuno@ucsc.edu)

from BlastClient import *
from BlastParser import *
from PfamClient import *

def main():
    print('BLASTing...')
    blaster = BlastClient()
    blastResults = blaster.doBlast()
    hsps = blaster.getHsps(blastResults)
    query = blaster.getQuery()

    print('Parsing...')
    parser = BlastParser()
    parser.getHspList(hsps)
    pdb = parser.getPdb()

    print('Pfaming...')
    pfamer = PfamClient()
    pfamResults = pfamer.doPfamSearch(pdb, query)
    conservedDomains = pfamer.getConservedDomains(pfamResults)
    print(conservedDomains)

if __name__ == '__main__':
    main()
