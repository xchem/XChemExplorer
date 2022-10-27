# last edited: 18/11/2016, 15:00

import XChemDB
import os
import sys
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'lib'))


def update_data_source(db_file, crystalString, status):
    db = XChemDB.data_source(db_file)

    print "update mainTable set PANDDAStatus = '{0!s}' where CrystalName in ({1!s})".format(
        status, "'"+crystalString.replace(",", "','")+"'")
    db.execute_statement("update mainTable set PANDDAStatus = '{0!s}' where CrystalName in ({1!s})".format(
        status, "'"+crystalString.replace(",", "','")+"'"))


if __name__ == '__main__':
    db_file = sys.argv[1]
    crystalString = sys.argv[2]
    status = sys.argv[3]

    update_data_source(db_file, crystalString, status)
