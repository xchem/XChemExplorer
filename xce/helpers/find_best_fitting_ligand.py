import glob
import os
import sys

import XChemDB


def parse_autofit_folder(sampleDir, compoundID):
    filenames = ["best.pdb", "ligand_fit_1.pdb"]
    pdbList = []
    for files in glob.glob(
        os.path.join(sampleDir, "autofit_ligand", compoundID + "_*", "*")
    ):
        pdb = files[files.rfind("/") + 1 :]
        if pdb in filenames:
            pdbList.append(files)
    return pdbList


def find_input_mtz(sampleDir):
    possibleMTZ = [
        "refine.mtz",
        "init.mtz",
        "dimple.mtz",
    ]
    mtzin = None
    for mtz in possibleMTZ:
        if os.path.isfile(os.path.join(sampleDir, mtz)):
            mtzin = os.path.join(sampleDir, mtz)
            break
    return mtzin


def calculate_cc(sampleDir, compoundID, db):
    xtal = sampleDir[sampleDir.rfind("/") + 1 :]
    pdbList = parse_autofit_folder(sampleDir, compoundID)
    mtzin = find_input_mtz(sampleDir)
    for pdb in pdbList:
        autofitDir = pdb[: pdb.rfind("/")]
        os.chdir(autofitDir)
        os.system(
            "phenix.get_cc_mtz_pdb %s %s > phenix.get_cc_mtz_pdb.log" % (pdb, mtzin)
        )
    new_file_links(pdbList)
    bestRun, bestCC = parse_cc_log(sampleDir, compoundID)
    changeLinks(sampleDir, pdbList, bestRun, compoundID)
    fitting_program = None
    if "phenix" in bestRun:
        fitting_program = "phenix.ligandfit"
    elif "rhofit" in bestRun:
        fitting_program = "rhofit"
    db_dict = {}
    db_dict["CompoundAutofitprogram"] = fitting_program
    db_dict["CompoundAutofitCC"] = str(bestCC)
    updateDB(db, db_dict, xtal)


def new_file_links(pdbList):
    for pdb in pdbList:
        autofitDir = pdb[: pdb.rfind("/")]
        autofitRun = pdb.split("/")[len(pdb.split("/")) - 2]
        pdbFile = pdb[pdb.rfind("/") + 1 :]
        os.chdir(autofitDir)
        if not os.path.isfile(autofitRun + ".pdb"):
            os.system("ln -s %s %s.pdb" % (pdbFile, autofitRun))


def changeLinks(sampleDir, pdbList, bestRun, compoundID):
    os.chdir(sampleDir)
    for pdb in pdbList:
        autofitDir = pdb.split("/")[len(pdb.split("/")) - 2]
        if autofitDir == bestRun:
            if os.path.isfile(pdb.replace(sampleDir, ".")):
                os.system("/bin/rm %s.pdb" % compoundID)
                os.system("ln -s %s %s.pdb" % (pdb.replace(sampleDir, "."), compoundID))


def parse_cc_log(sampleDir, compoundID):
    ccDict = {}
    ccList = []
    for ccLog in glob.glob(
        os.path.join(
            sampleDir, "autofit_ligand", compoundID + "_*", "phenix.get_cc_mtz_pdb.log"
        )
    ):
        autofitDir = ccLog.split("/")[len(ccLog.split("/")) - 2]
        for line in open(ccLog):
            if line.startswith("local CC:") and len(line.split()) > 2:
                cc = line.split()[2]
                ccList.append(float(cc))
                ccDict[autofitDir] = cc
                break
    bestRun = None
    try:
        bestCC = max(ccList)
        for ccRun in ccDict:
            print(ccRun, ccDict[ccRun], bestCC)
            if str(ccDict[ccRun]) == str(bestCC):
                bestRun = ccRun
                break
    except ValueError:
        pass
    return bestRun, bestCC


def updateDB(db, db_dict, xtal):
    db.update_data_source(xtal, db_dict)


if __name__ == "__main__":
    compoundID = sys.argv[1]
    sampleDir = sys.argv[2]
    database = sys.argv[3]
    db = XChemDB.data_source(database)

    calculate_cc(sampleDir, compoundID, db)
