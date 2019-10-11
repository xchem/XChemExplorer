import os,sys
from rdkit import Chem
#from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

def enumerateStereoChem(compoundID,pdb,sampleDir,database):
    stereosmiles = None
    mol = Chem.MolFromPDBFile(pdb)
    Chem.AssignStereochemistry(mol,cleanIt=True,force=True,flagPossibleStereoCenters=True)
    if Chem.FindMolChiralCenters(mol,includeUnassigned=True) == []:
        print 'no chiral centres found'
    else:
        stereosmiles = Chem.MolToSmiles(mol,isomericSmiles=True)
    generateRestraints(compoundID,sampleDir,database,stereosmiles)

def generateRestraints(compoundID,sampleDir,database,stereosmiles):
    cmd = 'phenix.elbow --smiles="%s" --chiral=enumerate --id=LIG --output=newTest' %y
    os.system(cmd)

def checkFiles(compoundID,sampleDir,database,stereosmiles):
    updateDB()

def updateDB():
    # update stereo field
    # update restraintsprogram field
    # update stereosmiles

if __name__=='__main__':
    compoundID = sys.argv[1]
    pdb = sys.argv[2]
    sampleDir = sys.argv[3]
    database = sys.argv[4]
    enumerateStereoChem(compoundID,pdb,sampleDir,database)

