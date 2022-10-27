#!/usr/local/python/python2.7.3/bin/python
import os
import sys


def FindTLSgroups(pdbFile):
    print('\n==> XCE @ helpers: running phenix.find_tls_groups on new pdb file')
    os.system(
        'phenix.find_tls_groups {0!s} > phenix.find_tls_groups.out'.format(pdbFile))
    GroupNames = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                  'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    OutRefmac = ''
    found = 0
    i = 0
    for line in open('phenix.find_tls_groups.out'):
        if found:
            if line.startswith('}'):
                break
            temp = [line.split()]
            OutRefmac = OutRefmac + '\nTLS    {0!s}\nRANGE '.format((GroupNames)[i]) + str(temp[0][3])[:-1] + ' ' + \
                str(temp[0][6]) + '.\'  ' + str(temp[0][3])[:-1] + \
                ' ' + str(temp[0][8]) + '.\' ALL\n'
            i += 1
        if line.startswith('refinement.refine.adp {'):
            found = 1
        f = open('refmac.tls', 'w')
        f.write(OutRefmac)
        f.close()

    OutPhenix = ''
    found = 0
    for line in open('phenix.find_tls_groups.out'):
        if found:
            OutPhenix = OutPhenix+line
            if line.startswith('}'):
                break
        if line.startswith('refinement.refine.adp {'):
            found = 1
            OutPhenix = OutPhenix+line
    f = open('phenix.tls', 'w')
    f.write(OutPhenix)
    f.close()


if __name__ == '__main__':
    pdbFile = sys.argv[1]
    FindTLSgroups(pdbFile)
