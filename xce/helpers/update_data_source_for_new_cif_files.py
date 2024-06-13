import os
import sys


if __name__ == "__main__":
    sys.path.insert(
        0, os.path.join(os.environ["XChemExplorer_DIR"], "dist", "xce-1.5.0-py2.7.egg")
    )
    from xce.lib import XChemDB

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
