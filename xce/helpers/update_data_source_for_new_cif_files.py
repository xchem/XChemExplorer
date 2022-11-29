import os
import sys

from ..lib import XChemDB

if __name__ == "__main__":
    db_file = sys.argv[1]
    xtal = sys.argv[2]
    initial_model_directory = sys.argv[3]
    compoundID = sys.argv[4]

    db = XChemDB.data_source(db_file)
    db_dict = {}
    if (
        os.path.isfile(
            os.path.join(initial_model_directory, xtal, "compound", compoundID + ".cif")
        )
        and os.path.getsize(
            os.path.join(initial_model_directory, xtal, "compound", compoundID + ".cif")
        )
        > 20
    ):
        db_dict["RefinementCIF"] = os.path.join(
            initial_model_directory, xtal, "compound", compoundID + ".cif"
        )
        db_dict["RefinementCIFStatus"] = "restraints generated"
    else:
        db_dict["RefinementCIF"] = ""
        db_dict["RefinementCIFStatus"] = "restraints failed"
        os.system(
            "/bin/rm "
            + os.path.join(
                initial_model_directory, xtal, compoundID.replace(" ", "") + ".pdb"
            )
        )
        os.system(
            "/bin/rm "
            + os.path.join(
                initial_model_directory, xtal, compoundID.replace(" ", "") + ".cif"
            )
        )
        os.system(
            "/bin/rm "
            + os.path.join(
                initial_model_directory, xtal, compoundID.replace(" ", "") + ".png"
            )
        )

    db.update_data_source(xtal, db_dict)
