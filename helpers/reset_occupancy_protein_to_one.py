import glob,os
import sys

def reset_residues(projectDir):
    for dirs in glob.glob(os.path.join(projectDir,'*')):
        xtal = dirs[dirs.rfind('/')+1:]
        os.chdir(dirs)
        for mol in toEdit:
            if os.path.isfile(mol):
                if not os.path.isfile(mol+'_original'):
                    os.system('/bin/cp %s %s' %(mol,mol+'_original'))
                print '%s: found %s' %(xtal,mol)
                out = ''
                changeResi = None
                for line in open(mol):
                    if line.startswith('ATOM'):
                        chainID=str(line[21:23]).replace(' ','')
                        resseq=str(line[23:26]).replace(' ','')
                        occupancy=str(line[56:60]).replace(' ','')
                        altLoc=str(line[16:17]).replace(' ','')
                        if altLoc != '':
                            out += line
                        else:
                            out += line[:56] + '1.00' + line[60:]
                    else:
                        out += line
                print 'saving file',mol,'...'
                newFile = open(mol,'w')
                newFile.write(out)
                newFile.close()


if __name__ == '__main__':
    projectDir = sys.argv[1]
    toEdit = ['refine.split.bound-state.pdb']
    reset_residues(projectDir)
