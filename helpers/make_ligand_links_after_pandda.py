# last edited: 10/08/2017, 15:00

import os
import glob
import sys

# projectDir='/dls/labxchem/data/2017/lb13385-109/processing/analysis/initial_model_cov'
# panddaDir='/dls/labxchem/data/2017/lb13385-109/processing/analysis/panddas_cov'


def make_links(projectDir, panddaDir):
    os.chdir(os.path.join(panddaDir, 'processed_datasets'))
    for xtal in glob.glob('*'):
        #        for cif in glob.glob(os.path.join(xtal,'ligand_files','*.cif')):
        #            print cif
        #            os.system('/bin/rm %s' %cif)
        #        for pdb in glob.glob(os.path.join(xtal,'ligand_files','*.pdb')):
        #            print pdb
        #            os.system('/bin/rm %s' %pdb)

        for pdb in glob.glob(os.path.join(projectDir, xtal, 'compound', '*.pdb')):
            os.system('ln -s %s %s/ligand_files' % (pdb, xtal))
        for cif in glob.glob(os.path.join(projectDir, xtal, 'compound', '*.cif')):
            os.system('ln -s %s %s/ligand_files' % (cif, xtal))


if __name__ == '__main__':
    projectDir = sys.argv[1]
    panddaDir = sys.argv[2]
