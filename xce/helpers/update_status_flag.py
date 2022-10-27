from ..lib import XChemDB
import sys


def update_data_source(db_file, xtal, db_column, status):
    db = XChemDB.data_source(db_file)
    db_dict = {db_column: status}
    db.update_data_source(xtal, db_dict)


if __name__ == "__main__":
    db_file = sys.argv[1]
    xtal = sys.argv[2]
    db_column = sys.argv[3]
    status = sys.argv[4]

    update_data_source(db_file, xtal, db_column, status)
