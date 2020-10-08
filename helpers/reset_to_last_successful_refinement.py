import os,glob
import sys

def find_refine_folder(projectDir):
    for xtalDir in sorted(glob.glob(os.path.join(projectDir,'*'))):
        xtal = xtalDir[xtalDir.rfind('/')+1:]
        os.chdir(xtalDir)
        dirList = []
        for dirs in sorted(glob.glob('*')):
            if dirs.startswith('Refine_') and os.path.isdir(dirs):
                if not os.path.isfile('refine.pdb'):
                    dirList.append([dirs,os.path.getctime(dirs)])
        if dirList != []:
            find_refine_pdb(xtal, dirList)

def find_refine_pdb(xtal, dirList):
    sorted_list = sorted(dirList, key=lambda x: x[1], reverse = True)
    foundRefinePDB = False
    for d in sorted_list:
        serial = int(d[0].split('_')[1])
        if os.path.isfile(os.path.join(d[0],'refine.pdb')):
            print xtal + ': found refine.pdb in ' + d[0]
            pdbLink = os.path.join(d[0],'refine.pdb')
            mtzLink = os.path.join(d[0],'refine.mtz')
            foundRefinePDB = True
            break
        elif os.path.isfile(os.path.join(d[0],'refine_%s.pdb' %str(serial))):
            print xtal + ': found refine_%s.pdb in %s' %(str(serial),d[0])
            foundRefinePDB = True
            pdbLink = os.path.join(d[0],'refine_%s.pdb' %str(serial))
            mtzLink = os.path.join(d[0],'refine_%s.mtz' %str(serial))
            break
    if not foundRefinePDB:
        print xtal + ': ERROR -> cannot find refine.pdb'
    else:
        print os.getcwd(),pdbLink, mtzLink
        reset_links(xtal, pdbLink, mtzLink)

def find_pandda_refine_pdb(xtal, dirList):
    print 'hallp'

def reset_links(xtal, pdbLink, mtzLink):
    print xtal + ': resetting links...'
    os.system('/bin/rm refine.pdb')
    os.system('/bin/rm refine.mtz')
    os.system('/bin/rm refine.split.bound-state.pdb')
    os.system('ln -s %s refine.pdb' %pdbLink)
    os.system('ln -s %s refine.mtz' %mtzLink)
    os.system('ln -s refine.pdb refine.split.bound-state.pdb')



if __name__ == '__main__':
    projectDir = sys.argv[1]
    find_refine_folder(projectDir)