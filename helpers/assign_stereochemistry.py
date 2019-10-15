import os,sys
from rdkit import Chem
#from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))
import XChemDB

def enumerateStereoChem(compoundID,sampleDir,db,xtal):
    stereosmiles = None
    if os.path.isfile(os.path.join(sampleDir,'compound',compoundID+'.pdb')):
        pdb = os.path.join(sampleDir,'compound',compoundID+'.pdb')
    else:
        print 'cannot find ' + compoundID + '.pdb'
        pass
    mol = Chem.MolFromPDBFile(pdb)
    Chem.AssignStereochemistry(mol,cleanIt=True,force=True,flagPossibleStereoCenters=True)
    if Chem.FindMolChiralCenters(mol,includeUnassigned=True) == []:
        print 'no chiral centres found'
    else:
        stereosmiles = Chem.MolToSmiles(mol,isomericSmiles=True)
        generateRestraints(compoundID,sampleDir,db,stereosmiles,xtal)

def generateRestraints(compoundID,sampleDir,db,stereosmiles,xtal):
    os.chdir(os.path.join(sampleDir,'compound'))
    cmd = 'phenix.elbow --smiles="%s" --chiral=enumerate --id=LIG --output=%s' %(stereosmiles,compoundID)
    os.system(cmd)
    checkFiles(compoundID,sampleDir,db,stereosmiles,xtal)


def checkFiles(compoundID,sampleDir,db,stereosmiles,xtal):
    db_dict = {}
    db_dict['CompoundStereo'] = 'TRUE'
    db_dict['CompoundStereoSMILES'] = stereosmiles
    db_dict['CompoundStereoCIFprogram'] = 'phenix.elbow'
    updateDB(db,db_dict,xtal)

def updateDB(db,db_dict,xtal):
    # update stereo field
    # update restraintsprogram field
    # update stereosmiles
    db.update_data_source(xtal, db_dict)

if __name__=='__main__':
    compoundID = sys.argv[1]
    sampleDir = sys.argv[2]
    xtal = sampleDir[sampleDir.rfind('/')+1:]
    database = sys.argv[3]
    db=XChemDB.data_source(database)
    enumerateStereoChem(compoundID,sampleDir,db,xtal)

