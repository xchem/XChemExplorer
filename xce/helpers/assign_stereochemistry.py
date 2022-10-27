from ..lib import XChemDB
import os
import sys
import glob
from rdkit import Chem

# from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions


def enumerateStereoChem(compoundID, sampleDir, db, xtal):
    # first need to recreate the original CIF file with phenix.elbow because we cannot be sure
    # which program was used to create the initial restrains
    # this is because funny things happen to aromatic rings in case the file was made with GRADE
    os.chdir(os.path.join(sampleDir, "compound"))
    sql = "select CompoundSMILES from mainTable where CrystalName = '%s'" % xtal
    print("sqlite query:", sql)
    query = db.execute_statement(sql)
    originalSMILES = query[0][0]
    print(originalSMILES)
    cmd = 'phenix.elbow --smiles="%s" --id=LIG --output=tmp' % (originalSMILES)
    print("current working directory:", os.getcwd())
    print("runnning:", cmd)
    os.system(cmd)

    stereosmiles = None
    pdb = None
    if os.path.isfile(os.path.join(sampleDir, "compound", "tmp.pdb")):
        pdb = os.path.join(sampleDir, "compound", "tmp.pdb")
    else:
        print("cannot find tmp.pdb")
        pass
    if pdb:
        mol = Chem.MolFromPDBFile(pdb)
        Chem.AssignStereochemistry(
            mol, cleanIt=True, force=True, flagPossibleStereoCenters=True
        )
        if Chem.FindMolChiralCenters(mol, includeUnassigned=True) == []:
            print("no chiral centres found")
            db_dict = {}
            db_dict["CompoundStereo"] = "FALSE"
            updateDB(db, db_dict, xtal)
        else:
            stereosmiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            generateRestraints(compoundID, sampleDir, db, stereosmiles, xtal)


def generateRestraints(compoundID, sampleDir, db, stereosmiles, xtal):
    cmd = 'phenix.elbow --smiles="%s" --chiral=enumerate --id=LIG --output=%s' % (
        stereosmiles,
        compoundID,
    )
    os.system(cmd)
    checkFiles(compoundID, sampleDir, db, stereosmiles, xtal)


def checkFiles(compoundID, sampleDir, db, stereosmiles, xtal):
    foundCIF = False
    allCIF = ""
    for cif in glob.glob(os.path.join(compoundID + "_*.cif")):
        foundCIF = True
        allCIF += cif + ";"
        print("---", cif, "===", allCIF)
    if foundCIF:
        db_dict = {}
        db_dict["CompoundStereo"] = "TRUE"
        db_dict["CompoundStereoSMILES"] = stereosmiles
        db_dict["CompoundStereoCIFprogram"] = "phenix.elbow"
        db_dict["CompoundStereoCIFs"] = allCIF[:-1]
        print(">>>", db_dict)
        updateDB(db, db_dict, xtal)


def updateDB(db, db_dict, xtal):
    # update stereo field
    # update restraintsprogram field
    # update stereosmiles
    db.update_data_source(xtal, db_dict)


if __name__ == "__main__":
    compoundID = sys.argv[1]
    sampleDir = sys.argv[2]
    xtal = sampleDir[sampleDir.rfind("/") + 1 :]
    print("sampleID:", xtal)
    database = sys.argv[3]
    print("database:", database)
    db = XChemDB.data_source(database)
    enumerateStereoChem(compoundID, sampleDir, db, xtal)
