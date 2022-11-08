from ..lib import XChemDB
from ..lib.XChemUtils import parse
from iotbx import mtz
import os
import sys


if __name__ == "__main__":
    db_file = sys.argv[1]
    xtal = sys.argv[2]
    inital_model_directory = sys.argv[3]

    db = XChemDB.data_source(db_file)
    if os.path.isfile(os.path.join(inital_model_directory, xtal, "dimple_twin.pdb")):
        db_dict = {
            "DimpleTwinPathToPDB": os.path.join(
                inital_model_directory, xtal, "dimple_twin.pdb"
            )
        }
        dimple_ran_successfully = False
        if os.path.isfile(
            os.path.join(inital_model_directory, xtal, "dimple_twin.mtz")
        ):
            db_dict["DimpleTwinPathToMTZ"] = os.path.join(
                inital_model_directory, xtal, "dimple_twin.mtz"
            )
            dimple_ran_successfully = True
            db_dict["DataProcessingDimpleTwinSuccessful"] = "True"
            db_dict["DimpleTwinStatus"] = "finished"
        if not dimple_ran_successfully:
            db_dict["DataProcessingDimpleTwinSuccessful"] = "False"
            db_dict["DimpleTwinStatus"] = "failed"
        pdb = parse().PDBheader(
            os.path.join(inital_model_directory, xtal, "dimple_twin.pdb")
        )
        db_dict["DimpleTwinRcryst"] = pdb["Rcryst"]
        db_dict["DimpleTwinRfree"] = pdb["Rfree"]
        db_dict["RefinementOutcome"] = "1 - Analysis Pending"
        db_dict["RefinementSpaceGroup"] = pdb["SpaceGroup"]
        db_dict["DimpleTwinFraction"] = pdb["TwinFraction"]

        # setting free.mtz file
        os.chdir(os.path.join(inital_model_directory, xtal))
        os.system("/bin/rm -f %s.free.mtz" % xtal)
        mtzFree = None
        db_dict["RefinementTwinMTZfree"] = ""
        if os.path.isfile(
            os.path.join(
                inital_model_directory,
                xtal,
                "dimple_twin",
                "dimple_rerun_on_selected_file",
                "dimple_twin",
                "prepared2.mtz",
            )
        ):
            mtzFree = os.path.join(
                inital_model_directory,
                xtal,
                "dimple_twin",
                "dimple_rerun_on_selected_file",
                "dimple_twin",
                "prepared2.mtz",
            )
        elif os.path.isfile(
            os.path.join(
                inital_model_directory,
                xtal,
                "dimple_twin",
                "dimple_rerun_on_selected_file",
                "dimple_twin",
                "prepared.mtz",
            )
        ):
            mtzFree = os.path.join(
                inital_model_directory,
                xtal,
                "dimple_twin",
                "dimple_rerun_on_selected_file",
                "dimple_twin",
                "prepared.mtz",
            )
        elif os.path.isfile(
            os.path.join(
                inital_model_directory,
                xtal,
                "dimple_twin",
                "dimple_twin",
                "prepared.mtz",
            )
        ):
            mtzFree = os.path.join(
                inital_model_directory,
                xtal,
                "dimple_twin",
                "dimple_twin",
                "prepared.mtz",
            )
        elif os.path.isfile(
            os.path.join(
                inital_model_directory,
                xtal,
                "dimple_twin",
                "dimple_twin",
                "prepared2.mtz",
            )
        ):
            mtzFree = os.path.join(
                inital_model_directory,
                xtal,
                "dimple_twin",
                "dimple_twin",
                "prepared2.mtz",
            )
        elif os.path.isfile(
            os.path.join(
                inital_model_directory,
                xtal,
                "dimple_twin",
                "dimple_rerun_on_selected_file",
                "dimple_twin",
                "free.mtz",
            )
        ):
            mtzFree = os.path.join(
                inital_model_directory,
                xtal,
                "dimple_twin",
                "dimple_rerun_on_selected_file",
                "dimple_twin",
                "free.mtz",
            )
        elif os.path.isfile(
            os.path.join(
                inital_model_directory, xtal, "dimple_twin", "dimple_twin", "free.mtz"
            )
        ):
            mtzFree = os.path.join(
                inital_model_directory, xtal, "dimple_twin", "dimple_twin", "free.mtz"
            )

        if mtzFree is not None:
            if "F_unique" in mtz.object(mtzFree).column_labels():
                cmd = (
                    "cad hklin1 %s hklout %s.free.mtz << eof\n" % (mtzFree, xtal)
                    + " monitor BRIEF\n"
                    " labin file 1 E1=F E2=SIGF E3=FreeR_flag\n"
                    " labout file 1 E1=F E2=SIGF E3=FreeR_flag\n"
                    "eof\n"
                )

                os.system(cmd)
            else:
                os.symlink(mtzFree, xtal + ".free.mtz")

            db_dict["RefinementTwinMTZfree"] = xtal + ".free.mtz"

        print("==> xce: updating data source after DIMPLE run")
        db.update_data_source(xtal, db_dict)

    else:
        # the actual dimple script creates symbolic links regardless if dimple was
        # successful or not python os.path.isfile is False if symbolic link points to
        # non existing file so we remove all of them!
        os.chdir(os.path.join(inital_model_directory, xtal))
        os.system("/bin/rm dimple_twin.pdb")
        os.system("/bin/rm dimple_twin.mtz")
        os.system("/bin/rm 2fofc_twin.map")
        os.system("/bin/rm fofc_twin.map")
