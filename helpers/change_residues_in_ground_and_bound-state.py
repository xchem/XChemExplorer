import glob,os
import sys

def replace_residues():
    for dirs in glob.glob(os.path.join(projectDir,'*')):
        xtal = dirs[dirs.rfind('/')+1:]
        os.chdir(dirs)
        if xtal == test:
            continue
        for mol in toEdit:
            if os.path.isfile(mol):
                if not os.path.isfile(mol+'_original'):
                    os.system('/bin/cp %s %s' %(mol,mol+'_original'))
                print '%s: found %s' %(xtal,mol)
                os.system('phenix.superpose_pdbs %s %s' %(mol,newModel))
                out = ''
                changeResi = None
                for line in open(mol):
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        chainID=str(line[21:23]).replace(' ','')
                        resseq=str(line[23:26]).replace(' ','')
                        occupancy=str(line[56:60]).replace(' ','')
                        resi = chainID + resseq
                        if changeResi != None:
                            if resi == changeResi:
                                print 'here',changeResi
                                continue
                            else:
                                changeResi = None
                        if resi in residues_to_change:
                            changeResi = resi
                            print 'changeResi',changeResi
                            for line in open(newModelName+'_fitted.pdb'):
                                if line.startswith('ATOM') or line.startswith('HETATM'):
                                    chainIDFitted=str(line[21:23]).replace(' ','')
                                    resseqFitted=str(line[23:26]).replace(' ','')
                                    residFitted = chainIDFitted + resseqFitted
                                    if residFitted == changeResi:
                                        out += line[:56] + occupancy + line[60:]
                        else:
                            out += line
                    else:
                        out+=line
                print 'saving file',mol,'...'
                newFile = open(mol,'w')
                newFile.write(out)
                newFile.close()


if __name__ == '__main__':
#    newModel = '/dls/labxchem/data/2019/lb22717-1/processing/INPP5DA/processing/reference/INPP5DA-x0080-ground-state-corrected.pdb'
#    residues_to_change = ['A567', 'A812', 'A851']
#    projectDir = '/dls/labxchem/data/2019/lb22717-1/processing/INPP5DA/processing/analysis/model_building'
#    toEdit = ['refine.split.ground-state.pdb', 'refine.split.bound-state.pdb']
#    test = 'INPP5DA-x0019'

    newModel = sys.argv[1]
    residues_to_change = sys.argv[2]
    projectDir = sys.argv[3]
    toEdit = ['refine.split.ground-state.pdb', 'refine.split.bound-state.pdb']
    newModelName = newModel[newModel.rfind('/') + 1:]

