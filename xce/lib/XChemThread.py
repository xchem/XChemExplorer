import csv
import glob
import math
import os
import pickle
from datetime import datetime

from PyQt4 import QtCore

from xce.lib import XChemDB
from xce.lib import XChemLog
from xce.lib import XChemMain
from xce.lib import XChemUtils
from xce.lib.cluster import slurm
from iotbx.reflection_file_reader import any_reflection_file


class synchronise_db_and_filesystem(QtCore.QThread):
    """
    - remove broken links
    - insert new samples in DB
    - update data for existing samples
    """

    def __init__(
        self, initial_model_directory, datasource, panddas_directory, xce_logfile, mode
    ):
        QtCore.QThread.__init__(self)
        self.initial_model_directory = initial_model_directory
        self.datasource = datasource
        self.db = XChemDB.data_source(self.datasource)
        self.all_samples_in_datasource = (
            self.db.get_all_samples_in_data_source_as_list()
        )
        self.panddas_directory = panddas_directory
        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.mode = mode

    def run(self):
        self.Logfile.insert("synchronising database and filesystem")
        self.Logfile.insert(
            "current project directory: " + self.initial_model_directory
        )

        XChemMain.backup_soakDB(self.datasource, self.xce_logfile)

        # get list of xtals

        self.xtal_list = []
        progress_step = 1
        progress = 0
        self.emit(QtCore.SIGNAL("update_progress_bar"), progress)
        if self.mode != "project_directory":
            # if only a single xtal is synched, then self.mode==xtalID
            self.Logfile.insert("synchronising " + self.mode + " only")
            self.xtal_list.append(self.mode)
        else:
            if len(glob.glob(os.path.join(self.initial_model_directory, "*"))) != 0:
                progress_step = 100 / float(
                    len(glob.glob(os.path.join(self.initial_model_directory, "*")))
                )
            self.Logfile.insert(
                "found "
                + str(len(glob.glob(os.path.join(self.initial_model_directory, "*"))))
                + " samples in project directory"
            )
            for directory in sorted(
                glob.glob(os.path.join(self.initial_model_directory, "*"))
            ):
                try:
                    os.chdir(directory)
                except OSError:
                    # this could happen if the user accidentaly left a file in the
                    # project directory
                    continue
                if os.listdir(directory) == []:
                    self.Logfile.warning(directory + " is empty; skipping...")
                    continue
                xtal = directory[directory.rfind("/") + 1 :]
                self.xtal_list.append(xtal)

        # go through list

        for xtal in self.xtal_list:
            self.Logfile.insert("directory name: " + xtal + " = sampleID in database")
            os.chdir(os.path.join(self.initial_model_directory, xtal))
            if xtal not in self.all_samples_in_datasource:
                self.Logfile.insert("sampleID not found in database: inserting " + xtal)
                self.db.execute_statement(
                    "insert into mainTable (CrystalName) values ('{0!s}');".format(xtal)
                )
                self.all_samples_in_datasource.append(xtal)

            db_dict = self.db.get_db_dict_for_sample(xtal)

            db_dict["ProjectDirectory"] = self.initial_model_directory

            db_dict = self.sync_data_processing(xtal, db_dict)

            db_dict = self.sync_dimple_results(xtal, db_dict)

            db_dict = self.sync_compound_information(db_dict)

            db_dict = self.sync_refinement_results(xtal, db_dict)

            if db_dict != {}:
                self.emit(
                    QtCore.SIGNAL("update_status_bar(QString)"),
                    "updating datasource for " + xtal,
                )
                self.db.update_data_source(xtal, db_dict)

            progress += progress_step
            self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        self.Logfile.insert("database mainTable update finished")
        self.Logfile.insert("updating panddaTable")
        self.emit(QtCore.SIGNAL("update_status_bar(QString)"), "updating panddaTable")
        self.sync_pandda_table_NEW()
        self.emit(
            QtCore.SIGNAL("update_status_bar(QString)"),
            "database panddaTable update finished",
        )
        self.Logfile.insert("database panddaTable update finished")

        self.emit(QtCore.SIGNAL("datasource_menu_reload_samples"))

    def sync_data_processing(self, xtal, db_dict):
        # AIMLESS logfile

        # in case the MTZ file which is used for refinement is different to the one
        # used for refinement
        if os.path.isfile("refine.mtz"):
            if os.path.isfile(xtal + ".free.mtz"):
                freeMTZ = XChemUtils.mtztools(xtal + ".free.mtz")
                nREFfree = freeMTZ.get_number_measured_reflections()
                if os.path.isfile(xtal + ".mtz"):
                    procMTZ = XChemUtils.mtztools(xtal + ".mtz")
                    nREF = procMTZ.get_number_measured_reflections()
                    (
                        CC,
                        errorMessage,
                    ) = freeMTZ.calculate_correlaton_between_intensities_in_mtzfiles(
                        xtal + ".mtz"
                    )
                    self.Logfile.insert(
                        "%s: calculating CC between %s.free.mtz"
                        " (%s refl) and %s.mtz (%s refl): %s"
                        % (xtal, xtal, str(nREFfree), xtal, str(nREF), str(CC))
                    )
                    if errorMessage != "":
                        self.Logfile.insert(
                            "pointless failed with the following error: %s"
                            % errorMessage
                        )

                    try:
                        if float(CC) < 0.999:
                            self.Logfile.insert(
                                "correlation coefficient between the two files is below"
                                " 0.999; will try to understand from dimple.log which"
                                " one was used for initial map calculation"
                            )
                            if os.path.isfile(
                                "dimple/dimple_rerun_on_selected_file/dimple/dimple.log"
                            ):
                                foundLine = False
                                mtzin = ""
                                for line in open(
                                    "dimple/dimple_rerun_on_selected_file/dimple/"
                                    "dimple.log"
                                ):
                                    if foundLine:
                                        mtzin = (
                                            line.replace(" ", "")
                                            .replace("\n", "")
                                            .replace("\r", "")
                                        )
                                        self.Logfile.insert(
                                            "%s was used for inital map calculation"
                                            % mtzin
                                        )
                                        break
                                    if line.startswith(" --no-cleanup"):
                                        foundLine = True

                                if os.path.isfile(mtzin):
                                    self.Logfile.insert(
                                        "%s: mtzfile used for refinement is not the"
                                        " same as the one chosen from autoprocessing"
                                        % xtal
                                    )
                                    self.Logfile.insert(
                                        "%s: current mtzfile after autoprocessing: %s"
                                        % (xtal, os.path.realpath(xtal + ".mtz"))
                                    )
                                    self.Logfile.insert(
                                        "%s: removing links for %s.mtz/%s.log"
                                        % (xtal, xtal, xtal)
                                    )
                                    os.system("/bin/rm %s.mtz 2> /dev/null" % xtal)
                                    os.system("/bin/rm %s.log 2> /dev/null" % xtal)
                                    self.Logfile.insert(
                                        "linking %s to %s.mtz"
                                        % (os.path.relpath(mtzin), xtal)
                                    )
                                    os.symlink(os.path.relpath(mtzin), xtal + ".mtz")
                                    for logfile in glob.glob(
                                        os.path.join(mtzin[: mtzin.rfind("/")], "*log")
                                    ):
                                        self.Logfile.insert(
                                            "linking %s to %s.log"
                                            % (os.path.relpath(logfile), xtal)
                                        )
                                        os.symlink(
                                            os.path.relpath(logfile), xtal + ".log"
                                        )
                                        break

                    except ValueError:
                        self.Logfile.insert(
                            "something went wrong: calculated CC value"
                            " does not seem to be a floating point number"
                        )

        found_logfile = False
        if os.path.isfile(xtal + ".log"):
            found_logfile = True
            db_dict["DataProcessingPathToLogfile"] = os.path.realpath(xtal + ".log")
            db_dict["DataProcessingLOGfileName"] = xtal + ".log"
            if (
                db_dict["DataCollectionOutcome"] == "None"
                or db_dict["DataCollectionOutcome"] == ""
            ):
                db_dict["DataCollectionOutcome"] = "success"
            aimless_results = XChemUtils.parse().read_aimless_logfile(
                db_dict["DataProcessingPathToLogfile"]
            )
            db_dict.update(aimless_results)
        else:
            db_dict["DataProcessingPathToLogfile"] = ""
            db_dict["DataProcessingLOGfileName"] = ""

        # MTZ file

        if os.path.isfile(xtal + ".mtz"):
            db_dict["DataProcessingPathToMTZfile"] = os.path.realpath(xtal + ".mtz")
            db_dict["DataProcessingMTZfileName"] = xtal + ".mtz"
            if not found_logfile:
                mtz_info = XChemUtils.mtztools(
                    xtal + ".mtz"
                ).get_information_for_datasource()
                db_dict.update(mtz_info)
                db_dict["DataCollectionOutcome"] = "success"
        else:
            db_dict["DataProcessingPathToMTZfile"] = ""
            db_dict["DataProcessingMTZfileName"] = ""

        return db_dict

    def sync_dimple_results(self, xtal, db_dict):
        # DIMPLE pdb

        if os.path.isfile("dimple.pdb"):
            db_dict["DimplePathToPDB"] = os.path.realpath("dimple.pdb")
            pdb_info = XChemUtils.parse().PDBheader("dimple.pdb")
            db_dict["DimpleRcryst"] = pdb_info["Rcryst"]
            db_dict["DimpleRfree"] = pdb_info["Rfree"]
            db_dict["DimpleResolutionHigh"] = pdb_info["ResolutionHigh"]
            db_dict["DimpleStatus"] = "finished"
        else:
            db_dict["DimplePathToPDB"] = ""
            db_dict["DimpleRcryst"] = ""
            db_dict["DimpleRfree"] = ""
            db_dict["DimpleResolutionHigh"] = ""
            db_dict["DimpleStatus"] = "pending"

        # DIMPLE mtz

        dimple_path = ""
        if os.path.isfile("dimple.mtz"):
            db_dict["DimplePathToMTZ"] = os.path.realpath("dimple.mtz")
            dimple_mtz = db_dict["DimplePathToMTZ"]
            dimple_path = dimple_mtz[: dimple_mtz.rfind("/")]
        else:
            db_dict["DimplePathToMTZ"] = ""
            db_dict["DimpleStatus"] = "pending"

        if os.path.isfile(
            os.path.join(
                dimple_path,
                "dimple",
                "dimple_rerun_on_selected_file",
                "dimple_run_in_progress",
            )
        ):
            db_dict["DimpleStatus"] = "running"

        # MTZ free file

        if os.path.isfile(xtal + ".free.mtz"):
            db_dict["RefinementMTZfree"] = os.path.realpath(xtal + ".free.mtz")
        else:
            db_dict["RefinementMTZfree"] = ""
            os.system("/bin/rm %s.free.mtz 2> /dev/null" % xtal)
            if os.path.isfile(os.path.join(dimple_path, "prepared2.mtz")):
                os.symlink(
                    os.path.relpath(os.path.join(dimple_path, "prepared2.mtz")),
                    xtal + ".free.mtz",
                )
                db_dict["RefinementMTZfree"] = os.path.realpath(xtal + ".free.mtz")
            elif os.path.isfile(os.path.join(dimple_path, "prepared.mtz")):
                os.symlink(
                    os.path.relpath(os.path.join(dimple_path, "prepared.mtz")),
                    xtal + ".free.mtz",
                )
                db_dict["RefinementMTZfree"] = os.path.realpath(xtal + ".free.mtz")
            elif os.path.isfile(os.path.join(dimple_path, "free.mtz")):
                os.symlink(
                    os.path.relpath(os.path.join(dimple_path, "free.mtz")),
                    xtal + ".free.mtz",
                )
                db_dict["RefinementMTZfree"] = os.path.realpath(xtal + ".free.mtz")

        return db_dict

    def sync_compound_information(self, db_dict):
        # only update database if SMILES or compoundID field is blank!

        compoundID = db_dict["CompoundCode"]
        if compoundID == "None" or compoundID == "":
            if os.path.isdir("compound"):
                for smiles in glob.glob("compound/*"):
                    if smiles.endswith("smiles"):
                        for line in open(smiles):
                            if len(line.split()) >= 1:
                                db_dict["CompoundCode"] = smiles[
                                    smiles.rfind("/") + 1 : smiles.rfind(".")
                                ]
                                compoundID = db_dict["CompoundCode"]
                                break

        if (
            os.path.isfile(compoundID + ".cif")
            and os.path.getsize(compoundID + ".cif") > 20
        ):
            db_dict["RefinementCIF"] = os.path.realpath(compoundID + ".cif").replace(
                os.getcwd() + "/", ""
            )
            db_dict["RefinementCIFStatus"] = "restraints generated"
        else:
            os.system("/bin/rm {0!s}.cif 2> /dev/null".format(compoundID))
            os.system("/bin/rm compound/{0!s}.cif 2> /dev/null".format(compoundID))
            db_dict["RefinementCIF"] = ""
            db_dict["RefinementCIFStatus"] = "pending"

        smilesDB = db_dict["CompoundSMILES"]
        smiles_found = True
        if smilesDB == "None" or smilesDB == "":
            smiles_found = False
            if os.path.isdir("compound"):
                for smiles in glob.glob("compound/*"):
                    if smiles.endswith("smiles"):
                        for line in open(smiles):
                            if len(line.split()) >= 1:
                                db_dict["CompoundSMILES"] = line.split()[0]
                                smilesDB = db_dict["CompoundSMILES"]
                                smiles_found = True
                                break

        if not smiles_found:
            db_dict["RefinementCIFStatus"] = "missing smiles"

        if (
            not os.path.isfile(compoundID + ".pdb")
            or os.path.getsize(compoundID + ".pdb") < 20
        ):
            os.system("/bin/rm {0!s}.pdb 2> /dev/null".format(compoundID))
            os.system("/bin/rm compound/{0!s}.pdb 2> /dev/null".format(compoundID))

        if (
            not os.path.isfile(compoundID + ".png")
            or os.path.getsize(compoundID + ".png") < 20
        ):
            os.system("/bin/rm {0!s}.png 2> /dev/null".format(compoundID))
            os.system("/bin/rm compound/{0!s}.png 2> /dev/null".format(compoundID))

        return db_dict

    def sync_refinement_results(self, xtal, db_dict):
        # REFINE pdb

        if os.path.isfile("refine.pdb"):
            db_dict["RefinementPDB_latest"] = os.path.realpath("refine.pdb")
            db_dict["RefinementStatus"] = "finished"
            pdb_info = XChemUtils.parse().dict_for_datasource_update("refine.pdb")
            db_dict.update(pdb_info)
            if (
                db_dict["RefinementOutcome"] == "None"
                or db_dict["RefinementOutcome"] == ""
            ):
                db_dict["RefinementOutcome"] = "3 - In Refinement"
            elif str(db_dict["RefinementOutcome"]).startswith("1"):
                db_dict["RefinementOutcome"] = "3 - In Refinement"
            elif str(db_dict["RefinementOutcome"]).startswith("2"):
                db_dict["RefinementOutcome"] = "3 - In Refinement"
        else:
            db_dict["RefinementPDB_latest"] = ""
            db_dict["RefinementStatus"] = "pending"
            db_dict["RefinementOutcome"] = "1 - Analysis Pending"
            os.system("/bin/rm refine.pdb 2> /dev/null")

        if os.path.isfile("REFINEMENT_IN_PROGRESS"):
            db_dict["RefinementStatus"] = "running"

        # REFINE bound pdb

        if os.path.isfile("refine.split.bound-state.pdb"):
            db_dict["RefinementBoundConformation"] = os.path.realpath(
                "refine.split.bound-state.pdb"
            )
        else:
            db_dict["RefinementBoundConformation"] = ""

        # REFINE mtz

        if os.path.isfile("refine.mtz"):
            db_dict["RefinementMTZ_latest"] = os.path.realpath("refine.mtz")
        else:
            db_dict["RefinementMTZ_latest"] = ""
            os.system("/bin/rm refine.mtz 2> /dev/null")

        return db_dict

    def find_apo_structures_for_PanDDA(self, panddaPATH):
        # first check if structure is already present in DB and if so if all the
        # information concur

        # need to update pandda directory for every exported structure so that
        # we know where to look for the pandda.log file that contains the relevant
        # information

        # update CrystalName_of_pandda_input in DB

        # in DB: update StructureType field accordingly

        # newer pandda versions seem to have severl copies of pandda.log with names like
        # pandda-2016-09-01-2139.log
        panddaLog = glob.glob(os.path.join(panddaPATH, "pandda*log"))
        panddaLog.sort(key=os.path.getmtime)

        readindApoStructures = False
        apoStructures = []
        apoStructureDict = {}
        apoString = ""
        for files in panddaLog:
            for line in open(files):
                if "No Statistical Maps Found:" in line:
                    readindApoStructures = True
                if readindApoStructures:
                    if "Pickling Object: processed_datasets" in line:
                        if line.split() >= 2:
                            # e.g. line.split() ->
                            # ['Pickling',
                            # 'Object:',
                            # 'processed_datasets/NUDT22A-x0055/pickles/dataset.pickle']
                            xtal = line.split()[2].split("/")[1]
                            if os.path.isfile(
                                os.path.join(
                                    panddaPATH,
                                    "processed_datasets",
                                    xtal,
                                    xtal + "-pandda-input.pdb",
                                )
                            ):
                                apoStructures.append(xtal)
                                apoString += xtal + ";"
                if (
                    "Pre-existing statistical maps (from previous runs) have been found"
                    " and will be reused:" in line
                ):
                    readindApoStructures = False
            apoStructureDict[panddaPATH] = apoStructures

        return apoString[:-1]

    def sync_pandda_table_NEW(self):
        progress_step = 1
        progress = 0
        self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        # also need to update PANDDA table...
        pandda_models = self.db.execute_statement(
            "select CrystalName,PANDDA_site_index,PANDDA_site_event_index,"
            "PANDDA_site_x,PANDDA_site_y,PANDDA_site_z,PANDDApath,ApoStructures"
            " from panddaTable"
        )
        if len(pandda_models) > 0:
            progress_step = 100 / float(len(pandda_models))
        else:
            self.Logfile.warning("panddaTable seems to be empty!")

        if pandda_models:
            for entry in pandda_models:
                db_pandda_dict = {}
                xtal = entry[0]
                site_index = entry[1]
                event_index = entry[2]
                panddaPATH = entry[6]
                apoStructures = entry[7]
                self.emit(
                    QtCore.SIGNAL("update_status_bar(QString)"),
                    "checking {0!s} -> site {1!s} -> event {2!s} ".format(
                        xtal, site_index, event_index
                    ),
                )
                self.emit(QtCore.SIGNAL("update_progress_bar"), progress)
                progress += progress_step

                try:
                    event_x = float(str(entry[3]))
                    event_y = float(str(entry[4]))
                    event_z = float(str(entry[5]))
                except ValueError:
                    pass

                # do not update pandda path since this one is updated during pandda
                # export! instead try to get apo semi-colon separated list of apo
                # structures that were used to calculate event maps; but only if field
                # is blank!
                if str(apoStructures) == "None" or apoStructures == "":
                    if panddaPATH != "None" or panddaPATH != "":
                        self.Logfile.insert(
                            "trying to find which apo structures were used to calculate"
                            " the event maps in " + panddaPATH
                        )
                        db_pandda_dict[
                            "ApoStructures"
                        ] = self.find_apo_structures_for_PanDDA(panddaPATH)
                    else:
                        self.Logfile.insert(
                            "pandda path for " + xtal + " is empty in database"
                        )

                # event map

                found_event_map = False
                db_pandda_dict["PANDDA_site_event_map"] = ""
                for file in glob.glob(
                    os.path.join(self.initial_model_directory, xtal, "*ccp4")
                ):
                    filename = file[file.rfind("/") + 1 :]
                    if filename.startswith(
                        xtal + "-event_" + event_index
                    ) and filename.endswith("map.native.ccp4"):
                        event_map = file
                        db_pandda_dict["PANDDA_site_event_map"] = os.path.realpath(
                            event_map
                        )
                        found_event_map = True
                        break
                if not found_event_map:
                    db_pandda_dict["PANDDA_site_event_map"] = ""

                db_pandda_dict["PANDDA_site_initial_model"] = ""
                for file in glob.glob(
                    os.path.join(self.initial_model_directory, xtal, "*pdb")
                ):
                    filename = file[file.rfind("/") + 1 :]
                    if filename.endswith("-ensemble-model.pdb"):
                        db_pandda_dict["PANDDA_site_initial_model"] = os.path.realpath(
                            file
                        ).replace(os.getcwd() + "/", "")
                        break

                db_pandda_dict["PANDDA_site_initial_mtz"] = ""
                for file in glob.glob(
                    os.path.join(self.initial_model_directory, xtal, "*mtz")
                ):
                    filename = file[file.rfind("/") + 1 :]
                    if filename.endswith("pandda-input.mtz"):
                        db_pandda_dict["PANDDA_site_initial_mtz"] = os.path.realpath(
                            file
                        ).replace(os.getcwd() + "/", "")
                        break

                db_pandda_dict["PANDDA_site_ligand_resname"] = ""
                db_pandda_dict["PANDDA_site_ligand_chain"] = ""
                db_pandda_dict["PANDDA_site_ligand_sequence_number"] = ""
                db_pandda_dict["PANDDA_site_ligand_altLoc"] = ""
                db_pandda_dict["PANDDA_site_ligand_placed"] = "False"
                db_pandda_dict["PANDDA_site_spider_plot"] = ""
                db_pandda_dict["PANDDA_site_ligand_id"] = ""

                db_pandda_dict["PANDDA_site_occupancy"] = ""
                db_pandda_dict["PANDDA_site_B_average"] = ""
                db_pandda_dict["PANDDA_site_B_ratio_residue_surroundings"] = ""
                db_pandda_dict["PANDDA_site_rmsd"] = ""
                db_pandda_dict["PANDDA_site_RSCC"] = ""
                db_pandda_dict["PANDDA_site_RSR"] = ""
                db_pandda_dict["PANDDA_site_RSZD"] = ""

                if os.path.isfile(
                    os.path.join(self.initial_model_directory, xtal, "refine.pdb")
                ):
                    ligands_in_file = XChemUtils.pdbtools(
                        os.path.join(self.initial_model_directory, xtal, "refine.pdb")
                    ).find_xce_ligand_details()
                    if not ligands_in_file:
                        self.Logfile.warning(
                            "{0!s}: could not find any ligands in refine.pdb".format(
                                xtal
                            )
                        )
                        continue
                    else:
                        self.Logfile.insert(
                            "{0!s}: found the following ligands in"
                            " refine.pdb: {1!s}".format(xtal, str(ligands_in_file))
                        )

                    distanceList = []
                    for ligand in ligands_in_file:
                        residue_name = ligand[0]
                        residue_chain = ligand[1]
                        residue_number = ligand[2]
                        residue_altLoc = ligand[3]
                        residue_xyz = XChemUtils.pdbtools(
                            os.path.join(
                                self.initial_model_directory, xtal, "refine.pdb"
                            )
                        ).get_center_of_gravity_of_residue_ish(
                            residue_chain, residue_number
                        )
                        distance = XChemUtils.calculate_distance_between_coordinates(
                            residue_xyz[0],
                            residue_xyz[1],
                            residue_xyz[2],
                            event_x,
                            event_y,
                            event_z,
                        )
                        distanceList.append(
                            [
                                distance,
                                residue_name,
                                residue_chain,
                                residue_number,
                                residue_altLoc,
                            ]
                        )
                        self.Logfile.insert(
                            "{0!s}: calculating distance between event and ligand"
                            " ({1!s} {2!s} {3!s}): {4!s}".format(
                                xtal,
                                residue_name,
                                residue_chain,
                                residue_number,
                                str(distance),
                            )
                        )

                    # now take the ligand that is closest to the event
                    try:
                        smallestDistance = min(distanceList, key=lambda x: x[0])
                    except ValueError:
                        self.Logfile.error(
                            "could not determine smallest distance between current"
                            " ligand pandda events"
                        )
                        continue
                    distance = smallestDistance[0]
                    residue_name = smallestDistance[1]
                    residue_chain = smallestDistance[2]
                    residue_number = smallestDistance[3]
                    residue_altLoc = smallestDistance[4]
                    self.Logfile.insert(
                        "{0!s}: ligand with the shorted distance ({1!s}A) to the"
                        " current event (id: {2!s}): {3!s} {4!s} {5!s} {6!s}".format(
                            xtal,
                            str(distance),
                            event_index,
                            residue_name,
                            residue_chain,
                            residue_number,
                            residue_altLoc,
                        )
                    )
                    db_pandda_dict["PANDDA_site_ligand_resname"] = residue_name
                    db_pandda_dict["PANDDA_site_ligand_chain"] = residue_chain
                    db_pandda_dict[
                        "PANDDA_site_ligand_sequence_number"
                    ] = residue_number
                    db_pandda_dict["PANDDA_site_ligand_altLoc"] = residue_altLoc
                    db_pandda_dict["PANDDA_site_ligand_placed"] = "True"
                    db_pandda_dict["PANDDA_site_ligand_id"] = (
                        residue_name + "-" + residue_chain + "-" + residue_number
                    )

                    if xtal + "/Refine_" in os.path.realpath(
                        os.path.join(self.initial_model_directory, xtal, "refine.pdb")
                    ):
                        tmp = os.path.realpath(
                            os.path.join(
                                self.initial_model_directory, xtal, "refine.pdb"
                            )
                        )
                        spider_plot = os.path.join(
                            tmp[: tmp.rfind("/")],
                            "residue_plots",
                            residue_chain + "-" + residue_number + ".png",
                        ).replace(" ", "")
                        if os.path.isfile(spider_plot):
                            db_pandda_dict[
                                "PANDDA_site_spider_plot"
                            ] = os.path.realpath(spider_plot)
                        if os.path.isfile(
                            os.path.join(tmp[: tmp.rfind("/")], "residue_scores.csv")
                        ):
                            with open(
                                os.path.join(
                                    tmp[: tmp.rfind("/")], "residue_scores.csv"
                                ),
                                "rb",
                            ) as csv_import:
                                csv_dict = csv.DictReader(csv_import)
                                for i, line in enumerate(csv_dict):
                                    residueNameChainNumber = line[""]
                                    if (
                                        residueNameChainNumber
                                        == residue_chain + "-" + residue_number
                                    ):
                                        db_pandda_dict["PANDDA_site_occupancy"] = line[
                                            "Occupancy"
                                        ]
                                        db_pandda_dict["PANDDA_site_B_average"] = line[
                                            "Average B-factor (Residue)"
                                        ]
                                        db_pandda_dict[
                                            "PANDDA_site_B_ratio_residue_surroundings"
                                        ] = line["Surroundings B-factor Ratio"]
                                        db_pandda_dict["PANDDA_site_rmsd"] = line[
                                            "Model RMSD"
                                        ]
                                        db_pandda_dict["PANDDA_site_RSCC"] = line[
                                            "RSCC"
                                        ]
                                        db_pandda_dict["PANDDA_site_RSR"] = line["RSR"]
                                        db_pandda_dict["PANDDA_site_RSZD"] = line[
                                            "RSZD"
                                        ]
                if db_pandda_dict != {}:
                    self.db.update_site_event_panddaTable(
                        xtal, site_index, event_index, db_pandda_dict
                    )
                    self.Logfile.insert(
                        "updating panddaTable for"
                        " xtal: {0!s}, site: {1!s}, event: {2!s}".format(
                            xtal, site_index, event_index
                        )
                    )
                    self.Logfile.insert("-> panddaDict: " + str(db_pandda_dict))


class create_png_and_cif_of_compound(QtCore.QThread):
    def __init__(
        self,
        external_software,
        initial_model_directory,
        compound_list,
        database_directory,
        data_source_file,
        todo,
        ccp4_scratch_directory,
        xce_logfile,
        max_queue_jobs,
        restraints_program,
        slurm_token,
    ):
        QtCore.QThread.__init__(self)
        self.external_software = external_software
        self.initial_model_directory = initial_model_directory
        self.compound_list = compound_list
        self.database_directory = database_directory
        self.data_source_file = data_source_file
        self.todo = todo
        self.ccp4_scratch_directory = ccp4_scratch_directory
        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.max_queue_jobs = max_queue_jobs
        self.db = XChemDB.data_source(
            os.path.join(self.database_directory, self.data_source_file)
        )
        self.restraints_program = restraints_program
        self.slurm_token = slurm_token

    def run(self):
        # first remove all ACEDRG input scripts in ccp4_scratch directory
        self.Logfile.insert(
            "removing all xce_acedrg scripts from " + self.ccp4_scratch_directory
        )
        os.chdir(self.ccp4_scratch_directory)
        os.system("/bin/rm -f xce_acedrg*")

        progress_step = 100 / float(len(self.compound_list))
        progress = 0
        counter = 1
        for item in self.compound_list:
            sampleID = item[0]
            compoundID = item[1]
            sm = item[2]
            # if counterions are present, split and take the longest string in list
            smiles = max(sm.split("."), key=len)
            self.emit(
                QtCore.SIGNAL("update_status_bar(QString)"),
                "creating ACEDRG shell script for " + sampleID,
            )
            if compoundID == "" or compoundID is None:
                compoundID = "compound"

            if not os.path.isdir(os.path.join(self.initial_model_directory, sampleID)):
                os.mkdir(os.path.join(self.initial_model_directory, sampleID))

            if self.todo == "ALL" or self.todo == "SELECTED":
                # remove symbolic links if present
                if os.path.isfile(
                    os.path.join(
                        self.initial_model_directory,
                        sampleID,
                        compoundID.replace(" ", "") + ".pdb",
                    )
                ):
                    os.system(
                        "/bin/rm "
                        + os.path.join(
                            self.initial_model_directory,
                            sampleID,
                            compoundID.replace(" ", "") + ".pdb",
                        )
                    )
                if os.path.isfile(
                    os.path.join(
                        self.initial_model_directory,
                        sampleID,
                        compoundID.replace(" ", "") + ".png",
                    )
                ):
                    os.system(
                        "/bin/rm "
                        + os.path.join(
                            self.initial_model_directory,
                            sampleID,
                            compoundID.replace(" ", "") + ".png",
                        )
                    )
                if os.path.isdir(
                    os.path.join(self.initial_model_directory, sampleID, "compound")
                ):
                    os.system(
                        "/bin/rm -fr "
                        + os.path.join(
                            self.initial_model_directory, sampleID, "compound"
                        )
                    )
                    db_dict = {"RefinementCIFStatus": "pending"}
                    self.Logfile.insert(
                        "{0!s}: removed compound directory and all its contents".format(
                            sampleID
                        )
                    )
                    self.Logfile.insert(
                        "{0!s}: setting RefinementCIFStatus flag to started".format(
                            sampleID
                        )
                    )
                    self.db.update_data_source(sampleID, db_dict)

            # create 'compound' directory if not present
            if not os.path.isdir(
                os.path.join(self.initial_model_directory, sampleID, "compound")
            ):
                os.mkdir(
                    os.path.join(self.initial_model_directory, sampleID, "compound")
                )

            # create text file which contains the smiles string
            if not os.path.isfile(
                os.path.join(
                    self.initial_model_directory,
                    sampleID,
                    "compound",
                    compoundID + ".smiles",
                )
            ):
                os.chdir(
                    os.path.join(self.initial_model_directory, sampleID, "compound")
                )
                f = open(compoundID + ".smiles", "w")
                f.write(smiles)
                f.close()

            if (
                not os.path.isfile(
                    os.path.join(
                        self.initial_model_directory,
                        sampleID,
                        compoundID.replace(" ", "") + ".cif",
                    )
                )
                or self.todo == "SELECTED"
            ):
                os.chdir(os.path.join(self.initial_model_directory, sampleID))
                os.system("/bin/rm -f %s*" % compoundID.replace(" ", ""))
                os.chdir(
                    os.path.join(self.initial_model_directory, sampleID, "compound")
                )

                XChemUtils.helpers().make_png(
                    self.initial_model_directory,
                    sampleID,
                    compoundID,
                    smiles,
                    self.external_software,
                    self.database_directory,
                    self.data_source_file,
                    self.ccp4_scratch_directory,
                    counter,
                    self.xce_logfile,
                    self.restraints_program,
                )
                counter += 1

                db_dict = {
                    "RefinementCIFprogram": self.restraints_program,
                    "RefinementCIFStatus": "started",
                }
                self.Logfile.insert(
                    "{0!s}: setting RefinementCIFStatus flag to started".format(
                        sampleID
                    )
                )
                self.db.update_data_source(sampleID, db_dict)

            progress += progress_step
            self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        # submit array job at Diamond
        self.Logfile.insert(
            "created input scripts for "
            + str(counter)
            + " ACEDRG jobs in "
            + self.ccp4_scratch_directory
        )
        os.chdir(self.ccp4_scratch_directory)
        self.Logfile.insert("changing directory to " + self.ccp4_scratch_directory)
        if counter > 1:
            Cmds = (
                "#!/bin/bash\n"
                + ". /etc/profile.d/modules.sh\n"
                + "./xce_{}_$SLURM_ARRAY_TASK_ID.sh\n".format(self.restraints_program)
            )
            f = open("%s_master.sh" % self.restraints_program, "w")
            f.write(Cmds)
            f.close()
            slurm.submit_cluster_job(
                str(self.restraints_program),
                "{!s}_master.sh".format(self.restraints_program),
                self.xce_logfile,
                self.slurm_token,
                array="0-{}".format(counter),
            )

        self.emit(QtCore.SIGNAL("datasource_menu_reload_samples"))


class fit_ligands(QtCore.QThread):
    def __init__(
        self,
        external_software,
        initial_model_directory,
        compound_list,
        database_directory,
        data_source_file,
        ccp4_scratch_directory,
        xce_logfile,
        max_queue_jobs,
        slurm_token,
    ):
        QtCore.QThread.__init__(self)
        self.external_software = external_software
        self.initial_model_directory = initial_model_directory
        self.compound_list = compound_list
        self.database_directory = database_directory
        self.data_source_file = data_source_file
        self.ccp4_scratch_directory = ccp4_scratch_directory
        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.max_queue_jobs = max_queue_jobs
        self.slurm_token = slurm_token
        self.db = XChemDB.data_source(
            os.path.join(self.database_directory, self.data_source_file)
        )
        self.n = 1

    def prepareInput(self, cmd, ligList):
        for cif in ligList:
            cmd += (
                "rhofit -m ../init.mtz -p ../init.pdb -l ../compound/%s.cif"
                " -scanchirals -d %s_rhofit\n" % (cif, cif)
            )
            cmd += (
                "phenix.ligandfit data=../init.mtz"
                " model=../init.pdb ligand=../compound/%s.cif clean_up=True\n" % cif
            )
            cmd += "/bin/mv LigandFit_run_1_ %s_phenix\n" % cif
            cmd += "/bin/rm -fr PDS\n\n"
        return cmd

    def get_header(self, sampleID):
        module = ""
        if os.path.isdir("/dls"):
            module = "module load phenix/1.20\n"
            module += "module load buster/20240123\n"

        cmd = (
            'export XChemExplorer_DIR="' + os.getenv("XChemExplorer_DIR") + '"\n'
            "\n"
            "cd %s\n"
            % os.path.join(self.initial_model_directory, sampleID, "autofit_ligand")
            + "\n"
            + module
        )

        return cmd

    def get_footer(self, cmd, sampleID, compoundID):
        cmd += (
            "$CCP4/bin/ccp4-python"
            " $XChemExplorer_DIR/xce/helpers/"
            "find_best_fitting_ligand.py {0!s} {1!s} {2!s}".format(
                compoundID.replace(" ", ""),
                os.path.join(self.initial_model_directory, sampleID),
                os.path.join(self.database_directory, self.data_source_file),
            )
        )
        return cmd

    def run(self):
        # first remove all ACEDRG input scripts in ccp4_scratch directory
        self.Logfile.insert(
            "removing all xce_fit_ligands scripts from " + self.ccp4_scratch_directory
        )
        os.chdir(self.ccp4_scratch_directory)
        os.system("/bin/rm -f xce_autofit_ligand*")

        progress_step = 100 / float(len(self.compound_list))
        progress = 0
        for item in self.compound_list:
            sampleID = item[0]
            compoundID = item[1]

            if os.path.isdir(
                os.path.join(self.initial_model_directory, sampleID, "autofit_ligand")
            ):
                os.system(
                    "/bin/rm -fr "
                    + os.path.join(
                        self.initial_model_directory, sampleID, "autofit_ligand"
                    )
                )
            os.mkdir(
                os.path.join(self.initial_model_directory, sampleID, "autofit_ligand")
            )

            # find ligands to fit (there might be more than one in case of stereoisomers
            ligList = []
            for file in glob.glob(
                os.path.join(
                    self.initial_model_directory,
                    sampleID,
                    "compound",
                    compoundID + "*.cif",
                )
            ):
                cif = file[file.rfind("/") + 1 :].replace(".cif", "")
                if "with_H" in cif:
                    continue
                else:
                    ligList.append(cif)

            if ligList != []:
                cmd = self.get_header(sampleID)
                cmd = self.prepareInput(cmd, ligList)
                cmd = self.get_footer(cmd, sampleID, compoundID)
                self.write_script(cmd)
                self.run_script()
                self.n += 1

            progress += progress_step
            self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        self.emit(QtCore.SIGNAL("datasource_menu_reload_samples"))

    def write_script(self, cmd):
        os.chdir(self.ccp4_scratch_directory)
        f = open("xce_autofit_ligand_{0!s}.sh".format(str(self.n)), "w")
        f.write(cmd)
        f.close()
        os.system("chmod +x xce_autofit_ligand_{0!s}.sh".format(str(self.n)))

    def run_script(self):
        # submit job
        self.Logfile.insert(
            "created input scripts for "
            + str(self.n)
            + " in "
            + self.ccp4_scratch_directory
        )
        os.chdir(self.ccp4_scratch_directory)
        Cmds = (
            "#!/bin/bash\n"
            + ". /etc/profile.d/modules.sh\n"
            + "./xce_autofit_ligand_$SLURM_ARRAY_TASK_ID.sh\n"
        )
        f = open("autofit_ligand_master.sh", "w")
        f.write(Cmds)
        f.close()
        slurm.submit_cluster_job(
            "xce_autofit_ligand_master",
            "autofit_ligand_master.sh",
            self.xce_logfile,
            self.slurm_token,
            array="1-{!s}".format(self.n - 1),
        )


class merge_cif_files(QtCore.QThread):
    def __init__(
        self, initial_model_directory, xce_logfile, second_cif_file, compound_list, todo
    ):
        QtCore.QThread.__init__(self)
        self.initial_model_directory = initial_model_directory
        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.second_cif_file = second_cif_file
        self.compound_list = compound_list
        self.todo = todo

    def run(self):
        progress_step = 100 / float(len(self.compound_list))
        progress = 0

        for item in self.compound_list:
            sampleID = item[0]
            compoundID = item[1]

            if os.path.isfile(
                os.path.join(
                    self.initial_model_directory,
                    sampleID,
                    "compound",
                    compoundID + ".cif",
                )
            ):
                self.Logfile.insert(
                    "%s: found %s.cif file in compound sub-directory"
                    % (sampleID, compoundID)
                )
            else:
                self.Logfile.error(
                    "%s: %s.cif file does not exist in compound sub-directory;"
                    " skipping..." % (sampleID, compoundID)
                )
                continue

            os.chdir(os.path.join(self.initial_model_directory, sampleID))
            if os.path.isfile(
                os.path.join(
                    self.initial_model_directory, sampleID, compoundID + ".cif"
                )
            ):
                self.Logfile.warning(
                    "%s: removing symbolic link to (or file) %s.cif"
                    " from sample directory" % (sampleID, compoundID)
                )
            os.system("/bin/rm %s.cif 2> /dev/null" % compoundID)

            if self.todo == "merge":
                self.emit(
                    QtCore.SIGNAL("update_status_bar(QString)"),
                    sampleID + " merging CIF files",
                )

                self.run_libcheck(sampleID, compoundID)

            elif self.todo == "restore":
                self.emit(
                    QtCore.SIGNAL("update_status_bar(QString)"),
                    sampleID + " restoring original CIF file",
                )
                self.Logfile.insert(
                    "%s: restoring symbolic link -> ln -s compound/%s.cif ."
                    % (sampleID, compoundID)
                )
                os.system("ln -s compound/%s.cif ." % compoundID)

            progress += progress_step
            self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        self.emit(QtCore.SIGNAL("finished()"))

    def run_libcheck(self, sampleID, compoundID):
        cmd = (
            "#!/bin/bash\n"
            "\n"
            ". /etc/profile.d/modules.sh\n"
            "$CCP4/bin/libcheck << eof \n"
            "_Y\n"
            "_FILE_L compound/%s.cif\n" % compoundID
            + "_FILE_L2 "
            + self.second_cif_file
            + "\n"
            "_FILE_O " + compoundID + ".cif\n"
            "_END\n"
            "eof\n"
        )
        self.Logfile.insert(
            "%s: running libcheck with the following input:\n%s" % (sampleID, cmd)
        )
        os.system(cmd)
        if os.path.isfile(compoundID + ".cif.lib"):
            self.Logfile.insert("%s: merged CIF file successfully created" % sampleID)
            os.system("/bin/mv %s %s" % (compoundID + ".cif.lib", compoundID + ".cif"))
        else:
            self.Logfile.error("%s: could not create merged CIF file" % sampleID)
            self.Logfile.warning(
                "%s: will re-create symbolic links to original restraints file"
                % sampleID
            )
            os.system("ln -s compound/%s.cif ." % compoundID)


class run_dimple_on_all_autoprocessing_files_new(QtCore.QThread):
    def __init__(
        self,
        sample_list,
        initial_model_directory,
        external_software,
        ccp4_scratch_directory,
        database_directory,
        data_source_file,
        max_queue_jobs,
        xce_logfile,
        dimple_twin_mode,
        pipeline,
        slurm_token,
    ):
        QtCore.QThread.__init__(self)
        self.sample_list = sample_list
        self.initial_model_directory = initial_model_directory
        self.external_software = external_software
        self.ccp4_scratch_directory = ccp4_scratch_directory
        self.database_directory = database_directory
        self.data_source_file = data_source_file
        self.db = XChemDB.data_source(
            os.path.join(self.database_directory, self.data_source_file)
        )
        self.max_queue_jobs = max_queue_jobs
        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.pipeline = pipeline
        self.dimple_twin_mode = dimple_twin_mode
        self.slurm_token = slurm_token

        self.n = 1

        self.Logfile.insert(
            "running initial refinement with the following pipeline: " + self.pipeline
        )

    def run(self):
        progress_step = 1
        if len(self.sample_list) != 0:
            progress_step = 100 / float(len(self.sample_list))
        progress = 0
        self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        os.chdir(self.ccp4_scratch_directory)
        os.system("/bin/rm -f xce_{0!s}*sh".format(self.pipeline))

        for item in sorted(self.sample_list):
            xtal = item[0]
            visit_run_autoproc = item[1]
            mtzin = item[2]
            ref_pdb = item[3]
            ref_mtz = item[4]
            ref_cif = item[5]

            if "dimple_rerun_on_selected_file" in visit_run_autoproc:
                if self.pipeline == "dimple":
                    twin = self.prepare_dimple_shell_script(
                        xtal, visit_run_autoproc, mtzin, ref_pdb, ref_mtz, ref_cif
                    )
                elif self.pipeline == "pipedream":
                    twin = self.prepare_pipedream_shell_script(
                        xtal, visit_run_autoproc, mtzin, ref_pdb, ref_mtz, ref_cif
                    )
                elif self.pipeline == "phenix.ligand_pipeline":
                    twin = self.prepare_phenix_ligand_pipeline_shell_script(
                        xtal, visit_run_autoproc, mtzin, ref_pdb, ref_mtz, ref_cif
                    )
            else:
                twin = self.prepare_dimple_shell_script(
                    xtal, visit_run_autoproc, mtzin, ref_pdb, ref_mtz, ref_cif
                )

            progress += progress_step
            self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        self.run_script(twin)

        self.emit(QtCore.SIGNAL("datasource_menu_reload_samples"))

    def prepare_phenix_ligand_pipeline_shell_script(
        self, xtal, visit_run_autoproc, mtzin, ref_pdb, ref_mtz, ref_cif
    ):
        # check if reference mtzfile has an Rfree column; if not, then ignore
        # DIMPLE assumes an Rfree column and barfs if it is not present
        # note: ref_mtz looks like this: ref mtz  -R reference.mtz

        if not os.path.isdir(os.path.join(self.initial_model_directory, xtal)):
            os.mkdir(os.path.join(self.initial_model_directory, xtal))
        os.chdir(os.path.join(self.initial_model_directory, xtal))
        if os.path.isdir(
            os.path.join(self.initial_model_directory, xtal, "phenix.ligand_pipeline")
        ):
            os.system("/bin/rm -fr phenix.ligand_pipeline")
        os.mkdir(
            os.path.join(self.initial_model_directory, xtal, "phenix.ligand_pipeline")
        )
        os.system("touch phenix.ligand_pipeline_run_in_progress")

        if "bash" in os.getenv("SHELL"):
            ccp4_scratch = "export CCP4_SCR=" + self.ccp4_scratch_directory + "\n"
        else:
            ccp4_scratch = ""

        if os.path.isdir("/dls"):
            ccp4_scratch += "module load phenix/1.20\n"
            ccp4_scratch += "module load ccp4/7.1.018\n"

        mtz_column_list = XChemUtils.mtztools(mtzin).get_all_columns_as_list()
        rfree = ""
        if "FreeR_flag" in mtz_column_list:
            rfree = ' xray_data.r_free_flags.label="FreeR_flag"'

        Cmds = (
            'export XChemExplorer_DIR="' + os.getenv("XChemExplorer_DIR") + '"\n'
            "\n"
            "cd %s\n"
            % os.path.join(self.initial_model_directory, xtal, "phenix.ligand_pipeline")
            + "\n"
            "module load global/cluster\n"
            "module load phenix/1.20\n"
            "module load buster/20240123\n"
            "\n" + ccp4_scratch + "\n"
            "$CCP4/bin/ccp4-python"
            " $XChemExplorer_DIR/xce/helpers/update_status_flag.py %s %s %s %s\n"
            % (
                os.path.join(self.database_directory, self.data_source_file),
                xtal,
                "DimpleStatus",
                "running",
            )
            + "\n"
            "phenix.ligand_pipeline %s %s" % (ref_pdb, mtzin) + " mr=False"
            " ligand_copies=0"
            " build=False"
            " prune=False"
            " remove_waters=False"
            " stop_if_r_free_greater_than=0.4"
            " update_waters=False" + rfree + " build_hydrogens=False\n"
            "\n"
            "fft hklin pipeline_1/refine_final.mtz mapout 2fofc.map << EOF\n"
            " labin F1=2FOFCWT_filled PHI=PH2FOFCWT\n"
            "EOF\n"
            "\n"
            "fft hklin pipeline_1/refine_final.mtz mapout fofc.map << EOF\n"
            " labin F1=FOFCWT PHI=PHFOFCWT\n"
            "EOF\n"
            "\n"
            "cd %s\n" % os.path.join(self.initial_model_directory, xtal) + "\n"
            "/bin/rm phenix.ligand_pipeline.pdb\n"
            "/bin/rm phenix.ligand_pipeline.mtz\n"
            "\n"
            "ln -s phenix.ligand_pipeline/pipeline_1/refine_final.pdb"
            " phenix.ligand_pipeline.pdb\n"
            "ln -s phenix.ligand_pipeline/pipeline_1/refine_final.mtz"
            " phenix.ligand_pipeline.mtz\n"
            "\n"
            "/bin/rm init.pdb\n"
            "/bin/rm init.mtz\n"
            "\n"
            "ln -s phenix.ligand_pipeline.pdb init.pdb\n"
            "ln -s phenix.ligand_pipeline.mtz init.mtz\n"
            "\n"
            "/bin/rm 2fofc.map\n"
            "/bin/rm fofc.map\n"
            "\n"
            "ln -s phenix.ligand_pipeline/2fofc.map .\n"
            "ln -s phenix.ligand_pipeline/fofc.map .\n"
            "\n"
            "$CCP4/bin/ccp4-python "
            + os.path.join(
                os.getenv("XChemExplorer_DIR"),
                "xce",
                "helpers",
                "update_data_source_for_new_dimple_pdb.py",
            )
            + " {0!s} {1!s} {2!s}\n".format(
                os.path.join(self.database_directory, self.data_source_file),
                xtal,
                self.initial_model_directory,
            )
            + "\n"
            "/bin/rm phenix.ligand_pipeline_run_in_progress\n"
        )

        os.chdir(self.ccp4_scratch_directory)
        f = open("xce_{0!s}_{1!s}.sh".format(self.pipeline, str(self.n)), "w")
        f.write(Cmds)
        f.close()
        os.system("chmod +x xce_{0!s}_{1!s}.sh".format(self.pipeline, str(self.n)))
        self.n += 1
        db_dict = {"DimpleStatus": "started"}
        self.Logfile.insert(
            "{0!s}: setting DataProcessingStatus flag to started".format(xtal)
        )
        self.db.update_data_source(xtal, db_dict)
        twin = ""
        return twin

    def prepare_pipedream_shell_script(
        self, xtal, visit_run_autoproc, mtzin, ref_pdb, ref_mtz, ref_cif
    ):
        if not os.path.isdir(os.path.join(self.initial_model_directory, xtal)):
            os.mkdir(os.path.join(self.initial_model_directory, xtal))
        if os.path.isdir(os.path.join(self.initial_model_directory, xtal, "pipedream")):
            os.chdir(os.path.join(self.initial_model_directory, xtal))
            os.system("/bin/rm -fr pipedream")
        os.mkdir(os.path.join(self.initial_model_directory, xtal, "pipedream"))
        os.system("touch pipedream_run_in_progress")

        if "bash" in os.getenv("SHELL"):
            ccp4_scratch = "export CCP4_SCR=" + self.ccp4_scratch_directory + "\n"
        else:
            ccp4_scratch = ""

        if os.path.isdir("/dls"):
            ccp4_scratch += "module load buster/20240123\n"
            ccp4_scratch += "module load ccp4/7.1.018\n"

        if os.path.isfile(ref_mtz):
            hklref_line = " -hklref {0!s}".format(ref_mtz)
        else:
            hklref_line = " -nofreeref"

        Cmds = (
            'export XChemExplorer_DIR="' + os.getenv("XChemExplorer_DIR") + '"\n'
            "\n"
            "cd %s\n" % os.path.join(self.initial_model_directory, xtal, "pipedream")
            + "\n"
            "module load global/cluster\n"
            "module load phenix/1.20\n"
            "module load buster/20240123\n"
            "\n" + ccp4_scratch + "\n"
            "$CCP4/bin/ccp4-python $XChemExplorer_DIR/xce/helpers/update_status_flag.py"
            " %s %s %s %s\n"
            % (
                os.path.join(self.database_directory, self.data_source_file),
                xtal,
                "DimpleStatus",
                "running",
            )
            + "\n"
            "pointless hklin {0!s} xyzin {1!s} hklout pointless.mtz >"
            " pointless.log\n".format(mtzin, ref_pdb) + "\n"
            "pipedream "
            " -d pipedreamDir"
            " -xyzin %s" % ref_pdb + hklref_line + " -hklin pointless.mtz"
            " -keepwater\n"
            "\n"
            "fft hklin pipedreamDir/refine/refine.mtz mapout 2fofc.map << EOF\n"
            " labin F1=2FOFCWT PHI=PH2FOFCWT\n"
            "EOF\n"
            "\n"
            "fft hklin pipedreamDir/refine/refine.mtz mapout fofc.map << EOF\n"
            " labin F1=FOFCWT PHI=PHFOFCWT\n"
            "EOF\n"
            "\n"
            "cd %s\n" % os.path.join(self.initial_model_directory, xtal) + "\n"
            "/bin/rm pipedream.pdb\n"
            "/bin/rm pipedream.mtz\n"
            "\n"
            "ln -s pipedream/pipedreamDir/refine/refine.pdb pipedream.pdb\n"
            "ln -s pipedream/pipedreamDir/refine/refine.mtz pipedream.mtz\n"
            "\n"
            "/bin/rm init.pdb\n"
            "/bin/rm init.mtz\n"
            "\n"
            "ln -s pipedream.pdb init.pdb\n"
            "ln -s pipedream.mtz init.mtz\n"
            "\n"
            "/bin/rm 2fofc.map\n"
            "/bin/rm fofc.map\n"
            "\n"
            "ln -s pipedream/2fofc.map .\n"
            "ln -s pipedream/fofc.map .\n"
            "\n"
            "$CCP4/libexec/python "
            + os.path.join(
                os.getenv("XChemExplorer_DIR"),
                "xce",
                "helpers",
                "update_data_source_for_new_dimple_pdb.py",
            )
            + " {0!s} {1!s} {2!s}\n".format(
                os.path.join(self.database_directory, self.data_source_file),
                xtal,
                self.initial_model_directory,
            )
            + "\n"
            "/bin/rm pipedream_run_in_progress\n"
        )

        os.chdir(self.ccp4_scratch_directory)
        f = open("xce_{0!s}_{1!s}.sh".format(self.pipeline, str(self.n)), "w")
        f.write(Cmds)
        f.close()
        os.system("chmod +x xce_{0!s}_{1!s}.sh".format(self.pipeline, str(self.n)))
        self.n += 1
        db_dict = {"DimpleStatus": "started"}
        self.Logfile.insert(
            "{0!s}: setting DataProcessingStatus flag to started".format(xtal)
        )
        self.db.update_data_source(xtal, db_dict)
        twin = ""
        return twin

    def prepare_dimple_shell_script(
        self, xtal, visit_run_autoproc, mtzin, ref_pdb, ref_mtz, ref_cif
    ):
        # check if reference mtzfile has an Rfree column; if not, then ignore
        # DIMPLE assumes an Rfree column and barfs if it is not present
        # note: ref_mtz looks like this: ref mtz  -R reference.mtz
        if os.path.isfile(ref_mtz):
            mtz_column_dict = XChemUtils.mtztools(ref_mtz).get_all_columns_as_dict()
            if "FreeR_flag" not in mtz_column_dict["RFREE"]:
                self.Logfile.insert(
                    "cannot find FreeR_flag in reference mtz file: {0!s} ->"
                    " ignoring reference mtzfile!!!".format(ref_mtz)
                )
                ref_mtz = ""
                if mtz_column_dict["RFREE"]:
                    self.Logfile.insert(
                        "found Rfree set with other column name though: {0!s}".format(
                            str(mtz_column_dict["RFREE"])
                        )
                    )
                    self.Logfile.insert(
                        "try renaming Rfree column to FreeR_flag with CAD!"
                    )

        db_dict = {"DimpleReferencePDB": ref_pdb}
        self.db.update_data_source(xtal, db_dict)

        self.emit(
            QtCore.SIGNAL("update_status_bar(QString)"),
            "creating input script for " + xtal + " in " + visit_run_autoproc,
        )

        if not os.path.isdir(os.path.join(self.initial_model_directory, xtal)):
            os.mkdir(os.path.join(self.initial_model_directory, xtal))
        os.chdir(os.path.join(self.initial_model_directory, xtal))

        twin = ""
        twinRefmac = ""
        if self.dimple_twin_mode:
            twinRefmac = "--refmac-key 'TWIN'"
            twin = "_twin"
            if os.path.isdir(
                os.path.join(self.initial_model_directory, xtal, "dimple_twin")
            ):
                os.system("/bin/rm -fr dimple_twin")
            os.mkdir(os.path.join(self.initial_model_directory, xtal, "dimple_twin"))
            os.system("touch dimple_twin_run_in_progress")
        else:
            if os.path.isdir(
                os.path.join(self.initial_model_directory, xtal, "dimple")
            ):
                os.system("/bin/rm -fr dimple")
            os.mkdir(os.path.join(self.initial_model_directory, xtal, "dimple"))
            os.system("touch dimple_run_in_progress")
        os.system("/bin/rm final.mtz 2> /dev/null")
        os.system("/bin/rm final.pdb 2> /dev/null")

        if "bash" in os.getenv("SHELL"):
            ccp4_scratch = "export CCP4_SCR=" + self.ccp4_scratch_directory + "\n"
        else:
            ccp4_scratch = ""

        if os.path.isdir("/dls"):
            ccp4_scratch += "module load ccp4/7.1.018\n"

        hkl = any_reflection_file(file_name=mtzin)
        miller_arrays = hkl.as_miller_arrays()
        mtzFile = miller_arrays[0]

        if mtzFile.space_group_info().symbol_and_number() == "R 3 :H (No. 146)":
            symNoAbsence = "H3"
        elif mtzFile.space_group_info().symbol_and_number() == "R 3 2 :H (No. 155)":
            symNoAbsence = "H32"
        else:
            symNoAbsence = (
                str(
                    [
                        x[0]
                        for x in str(
                            mtzFile.space_group_info().symbol_and_number().split("(")[0]
                        ).split()
                    ]
                )
                .replace("[", "")
                .replace("]", "")
                .replace("'", "")
                .replace(",", "")
                .replace(" ", "")
            )

        Cmds = (
            'export XChemExplorer_DIR="' + os.getenv("XChemExplorer_DIR") + '"\n'
            "\n"
            "cd %s\n"
            % os.path.join(self.initial_model_directory, xtal, "dimple%s" % twin)
            + "\n"
            "module load global/cluster\n"
            "module load phenix/1.20\n"
            "module load buster/20240123\n"
            "\n" + ccp4_scratch + "\n"
            "$CCP4/bin/ccp4-python $XChemExplorer_DIR/xce/helpers/update_status_flag.py"
            " %s %s %s %s\n"
            % (
                os.path.join(self.database_directory, self.data_source_file),
                xtal,
                "DimpleStatus",
                "running",
            )
            + "\n"
            "cd %s\n"
            % os.path.join(self.initial_model_directory, xtal, "dimple%s" % twin)
            + "\n"
            "unique hklout unique.mtz << eof\n"
            " cell %s\n"
            % str([round(float(i), 2) for i in mtzFile.unit_cell().parameters()])
            .replace("[", "")
            .replace("]", "")
            + " symmetry %s\n" % symNoAbsence
            + " resolution %s\n" % str(round(float(mtzFile.d_min()), 3))
            + "eof\n"
            "\n"
            "\n"
            "sftools << eof > sftools.log\n"
            " read unique.mtz\n"
            " calc col F = 10.0\n"
            " calc col SIGF = 1.0\n"
            " write sftools.mtz\n"
            "eof\n"
            "\n"
            "cad hklin1 sftools.mtz hklin2 %s hklout %s.999A.mtz << eof\n"
            % (mtzin, xtal)
            + " monitor BRIEF\n"
            " labin file 1 E1=F E2=SIGF\n"
            " labout file 1 E1=F_unique E2=SIGF_unique\n"
            " labin file 2 ALL\n"
            " resolution file 1 999.0 %s\n" % str(round(float(mtzFile.d_min()), 2))
            + "eof\n"
            "\n"
            "pointless hklin %s.999A.mtz hklout %s.999A.reind.mtz xyzin %s << eof >"
            " pointless.reind.log\n" % (xtal, xtal, ref_pdb) + " tolerance 5\n"
            "eof\n"
            "\n"
            "dimple --no-cleanup %s.999A.reind.mtz %s %s %s %s dimple%s\n"
            % (xtal, ref_pdb, ref_mtz, ref_cif, twinRefmac, twin)
            + "\n"
            "fft hklin dimple%s/final.mtz mapout 2fofc%s.map << EOF\n" % (twin, twin)
            + " labin F1=FWT PHI=PHWT\n"
            "EOF\n"
            "\n"
            "fft hklin dimple%s/final.mtz mapout fofc%s.map << EOF\n" % (twin, twin)
            + " labin F1=DELFWT PHI=PHDELWT\n"
            "EOF\n"
            "\n"
            "cd %s\n" % os.path.join(self.initial_model_directory, xtal) + "\n"
            "/bin/rm dimple%s.pdb\n" % twin + "/bin/rm dimple%s.mtz\n" % twin + "\n"
            "ln -s dimple%s/dimple%s/final.pdb dimple%s.pdb\n" % (twin, twin, twin)
            + "ln -s dimple%s/dimple%s/final.mtz dimple%s.mtz\n" % (twin, twin, twin)
            + "\n"
            "/bin/rm init%s.pdb\n" % twin + "/bin/rm init%s.mtz\n" % twin + "\n"
            "ln -s dimple%s.pdb init%s.pdb\n" % (twin, twin)
            + "ln -s dimple%s.mtz init%s.mtz\n" % (twin, twin)
            + "\n"
            "/bin/rm 2fofc%s.map\n" % twin + "/bin/rm fofc%s.map\n" % twin + "\n"
            "ln -s dimple%s/2fofc%s.map .\n" % (twin, twin)
            + "ln -s dimple%s/fofc%s.map .\n" % (twin, twin)
            + "\n"
            "$CCP4/bin/ccp4-python "
            + os.path.join(
                os.getenv("XChemExplorer_DIR"),
                "xce",
                "helpers",
                "update_data_source_for_new_dimple%s_pdb.py" % twin,
            )
            + " {0!s} {1!s} {2!s}\n".format(
                os.path.join(self.database_directory, self.data_source_file),
                xtal,
                self.initial_model_directory,
            )
            + "\n"
            "/bin/rm dimple_run_in_progress\n"
            "\n"
        )

        os.chdir(self.ccp4_scratch_directory)
        f = open(
            "xce_{0!s}{1!s}_{2!s}.sh".format(self.pipeline, twin, str(self.n)), "w"
        )
        f.write(Cmds)
        f.close()
        os.system(
            "chmod +x xce_{0!s}{1!s}_{2!s}.sh".format(self.pipeline, twin, str(self.n))
        )
        self.n += 1
        db_dict = {"DimpleStatus": "started"}
        self.Logfile.insert(
            "{0!s}: setting DataProcessingStatus flag to started".format(xtal)
        )
        self.db.update_data_source(xtal, db_dict)
        return twin

    def run_script(self, twin):
        # submit job
        self.Logfile.insert(
            "created input scripts for "
            + str(self.n)
            + " in "
            + self.ccp4_scratch_directory
        )
        os.chdir(self.ccp4_scratch_directory)
        Cmds = (
            "#!/bin/bash\n"
            + ". /etc/profile.d/modules.sh\n"
            + "./xce_{!s}{!s}_$SLURM_ARRAY_TASK_ID.sh\n".format(self.pipeline, twin)
        )
        f = open("{!s}{!s}_master.sh".format(self.pipeline, twin), "w")
        f.write(Cmds)
        f.close()
        slurm.submit_cluster_job(
            "xce_{!s}{!s}_master".format(self.pipeline, twin),
            "{!s}{!s}_master.sh".format(self.pipeline, twin),
            self.xce_logfile,
            self.slurm_token,
            array="1-{!s}".format(self.n - 1),
        )


class remove_selected_dimple_files(QtCore.QThread):
    def __init__(
        self,
        sample_list,
        initial_model_directory,
        xce_logfile,
        database_directory,
        data_source_file,
        pipeline,
    ):
        QtCore.QThread.__init__(self)
        self.sample_list = sample_list
        self.initial_model_directory = initial_model_directory
        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.db = XChemDB.data_source(
            os.path.join(database_directory, data_source_file)
        )
        self.pipeline = pipeline

    def run(self):
        progress_step = 1
        if len(self.sample_list) != 0:
            progress_step = 100 / float(len(self.sample_list))
        progress = 0
        self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        for n, xtal in enumerate(self.sample_list):
            db_dict = {}
            if not os.path.isdir(os.path.join(self.initial_model_directory, xtal)):
                self.Logfile.insert("{0!s}: directory does not exist".format(xtal))
                continue
            os.chdir(os.path.join(self.initial_model_directory, xtal))

            if self.pipeline == "dimple":
                if os.path.isfile("init.pdb"):
                    if "dimple" in os.path.realpath("init.pdb"):
                        self.Logfile.warning(
                            "{0!s}: init.pdb & init.mtz is linked to"
                            " dimple outcome".format(xtal)
                        )
                        self.Logfile.warning(
                            "{0!s}: removing init.pdb & init.mtz & (2)fofc maps".format(
                                xtal
                            )
                        )
                        db_dict = self.remove_init(db_dict)
                else:
                    db_dict = self.remove_init(db_dict)
                self.Logfile.warning(
                    "{0!s}: removing dimple folder & dimple.pdb/dimple.mtz".format(xtal)
                )
                os.system("/bin/rm dimple_run_in_progress 2> /dev/null")
                os.system("/bin/rm dimple.pdb 2> /dev/null")
                os.system("/bin/rm dimple.mtz 2> /dev/null")
                os.system("/bin/rm -fr dimple")
            elif self.pipeline == "pipedream":
                if os.path.isfile("init.pdb"):
                    if "dimple" in os.path.realpath("init.pdb"):
                        self.Logfile.warning(
                            "{0!s}: init.pdb & init.mtz is linked to"
                            " pipedream outcome".format(xtal)
                        )
                        self.Logfile.warning(
                            "{0!s}: removing init.pdb & init.mtz & (2)fofc maps".format(
                                xtal
                            )
                        )
                        db_dict = self.remove_init(db_dict)
                else:
                    db_dict = self.remove_init(db_dict)
                self.Logfile.warning(
                    "{0!s}: removing pipedream folder &"
                    " pipedream.pdb/pipedream.mtz".format(xtal)
                )
                os.system("/bin/rm pipedream_run_in_progress 2> /dev/null")
                os.system("/bin/rm pipedream.pdb 2> /dev/null")
                os.system("/bin/rm pipedream.mtz 2> /dev/null")
                os.system("/bin/rm -fr pipedream")
            elif self.pipeline == "phenix.ligand_pipeline":
                if os.path.isfile("init.pdb"):
                    if "dimple" in os.path.realpath("init.pdb"):
                        self.Logfile.warning(
                            "{0!s}: init.pdb & init.mtz is linked to"
                            " phenix.ligand_pipeline outcome".format(xtal)
                        )
                        self.Logfile.warning(
                            "{0!s}: removing init.pdb & init.mtz & (2)fofc maps".format(
                                xtal
                            )
                        )
                        db_dict = self.remove_init(db_dict)
                else:
                    db_dict = self.remove_init(db_dict)
                self.Logfile.warning(
                    "{0!s}: removing phenix.ligand_pipeline folder &"
                    " phenix.ligand_pipeline.pdb/phenix.ligand_pipeline.mtz".format(
                        xtal
                    )
                )
                os.system("/bin/rm phenix.ligand_pipeline_run_in_progress 2> /dev/null")
                os.system("/bin/rm phenix.ligand_pipeline.pdb 2> /dev/null")
                os.system("/bin/rm phenix.ligand_pipeline.mtz 2> /dev/null")
                os.system("/bin/rm -fr phenix.ligand_pipeline")

            if db_dict != {}:
                self.Logfile.insert("{0!s}: updating database".format(xtal))
                self.db.update_data_source(xtal, db_dict)

            progress += progress_step
            self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        self.emit(QtCore.SIGNAL("datasource_menu_reload_samples"))

    def remove_init(self, db_dict):
        os.system("/bin/rm init.pdb")
        os.system("/bin/rm init.mtz")
        os.system("/bin/rm 2fofc.map")
        os.system("/bin/rm fofc.map")
        db_dict["DimpleResolutionHigh"] = ""
        db_dict["DimpleRcryst"] = ""
        db_dict["DimpleRfree"] = ""
        db_dict["DimplePathToPDB"] = ""
        db_dict["DimplePathToMTZ"] = ""
        db_dict["DimpleReferencePDB"] = ""
        db_dict["DimplePANDDAwasRun"] = "False"
        db_dict["DimplePANDDAhit"] = "False"
        db_dict["DimplePANDDAreject"] = "False"
        db_dict["DimplePANDDApath"] = ""
        db_dict["DimpleStatus"] = "pending"
        return db_dict


class set_results_from_selected_pipeline(QtCore.QThread):
    def __init__(
        self,
        sample_list,
        initial_model_directory,
        xce_logfile,
        database_directory,
        data_source_file,
        pipeline,
    ):
        QtCore.QThread.__init__(self)
        self.sample_list = sample_list
        self.initial_model_directory = initial_model_directory
        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.db = XChemDB.data_source(
            os.path.join(database_directory, data_source_file)
        )
        self.pipeline = pipeline

    def run(self):
        progress_step = 1
        if len(self.sample_list) != 0:
            progress_step = 100 / float(len(self.sample_list))
        progress = 0
        self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        for n, xtal in enumerate(self.sample_list):
            db_dict = {}
            if not os.path.isdir(os.path.join(self.initial_model_directory, xtal)):
                self.Logfile.insert("{0!s}: directory does not exist".format(xtal))
                continue
            os.chdir(os.path.join(self.initial_model_directory, xtal))

            if os.path.isfile("init.pdb"):
                self.Logfile.warning("{0!s}: init.pdb & init.mtz exist".format(xtal))
                self.Logfile.warning(
                    "{0!s}: removing init.pdb & init.mtz & (2)fofc maps".format(xtal)
                )
                os.system("/bin/rm init.pdb")
                os.system("/bin/rm init.mtz")
                os.system("/bin/rm 2fofc.map")
                os.system("/bin/rm fofc.map")
                db_dict["DimpleResolutionHigh"] = ""
                db_dict["DimpleRcryst"] = ""
                db_dict["DimpleRfree"] = ""
                db_dict["DimplePathToPDB"] = ""
                db_dict["DimplePathToMTZ"] = ""
                db_dict["DimpleReferencePDB"] = ""
                db_dict["DimplePANDDAwasRun"] = "False"
                db_dict["DimplePANDDAhit"] = "False"
                db_dict["DimplePANDDAreject"] = "False"
                db_dict["DimplePANDDApath"] = ""
                db_dict["DimpleStatus"] = "pending"

            if self.pipeline == "dimple":
                if os.path.isfile("dimple.pdb"):
                    self.Logfile.insert("%s: selecting output from dimple" % xtal)
                    os.system("ln -s dimple.pdb init.pdb")
                    pdb_info = XChemUtils.parse().PDBheader("dimple.pdb")
                    db_dict["DimpleRcryst"] = pdb_info["Rcryst"]
                    db_dict["DimpleRfree"] = pdb_info["Rfree"]
                    db_dict["DimpleResolutionHigh"] = pdb_info["ResolutionHigh"]
                    db_dict["DimpleStatus"] = "finished"
                    db_dict["DimplePathToPDB"] = os.path.realpath("dimple.pdb")
                if os.path.isfile("dimple.mtz"):
                    os.system("ln -s dimple.mtz init.mtz")
                    db_dict["DimplePathToMTZ"] = os.path.realpath("dimple.mtz")
                if os.path.isfile("dimple/2fofc.map"):
                    os.system("ln -s dimple/2fofc.map .")
                if os.path.isfile("dimple/fofc.map"):
                    os.system("ln -s dimple/fofc.map .")
            elif self.pipeline == "pipedream":
                if os.path.isfile("pipedream.pdb"):
                    self.Logfile.insert("%s: selecting output from pipedream" % xtal)
                    os.system("ln -s pipedream.pdb init.pdb")
                    pdb_info = XChemUtils.parse().PDBheader("pipedream.pdb")
                    db_dict["DimpleRcryst"] = pdb_info["Rcryst"]
                    db_dict["DimpleRfree"] = pdb_info["Rfree"]
                    db_dict["DimpleResolutionHigh"] = pdb_info["ResolutionHigh"]
                    db_dict["DimpleStatus"] = "finished"
                    db_dict["DimplePathToPDB"] = os.path.realpath("pipedream.pdb")
                if os.path.isfile("pipedream.mtz"):
                    os.system("ln -s pipedream.mtz init.mtz")
                    db_dict["DimplePathToMTZ"] = os.path.realpath("pipedream.mtz")
                if os.path.isfile("pipedream/2fofc.map"):
                    os.system("ln -s pipedream/2fofc.map .")
                if os.path.isfile("pipedream/fofc.map"):
                    os.system("ln -s pipedream/fofc.map .")
            elif self.pipeline == "phenix.ligand_pipeline":
                if os.path.isfile("phenix.ligand_pipeline.pdb"):
                    self.Logfile.insert(
                        "%s: selecting output from phenix.ligand_pipeline" % xtal
                    )
                    os.system("ln -s phenix.ligand_pipeline.pdb init.pdb")
                    pdb_info = XChemUtils.parse().PDBheader(
                        "phenix.ligand_pipeline.pdb"
                    )
                    db_dict["DimpleRcryst"] = pdb_info["Rcryst"]
                    db_dict["DimpleRfree"] = pdb_info["Rfree"]
                    db_dict["DimpleResolutionHigh"] = pdb_info["ResolutionHigh"]
                    db_dict["DimpleStatus"] = "finished"
                    db_dict["DimplePathToPDB"] = os.path.realpath(
                        "phenix.ligand_pipeline.pdb"
                    )
                if os.path.isfile("phenix.ligand_pipeline.mtz"):
                    os.system("ln -s phenix.ligand_pipeline.mtz init.mtz")
                    db_dict["DimplePathToMTZ"] = os.path.realpath(
                        "phenix.ligand_pipeline.mtz"
                    )
                if os.path.isfile("phenix.ligand_pipeline/2fofc.map"):
                    os.system("ln -s phenix.ligand_pipeline/2fofc.map .")
                if os.path.isfile("phenix.ligand_pipeline/fofc.map"):
                    os.system("ln -s phenix.ligand_pipeline/fofc.map .")

            self.Logfile.insert("{0!s}: updating database".format(xtal))
            self.db.update_data_source(xtal, db_dict)

            progress += progress_step
            self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        self.emit(QtCore.SIGNAL("datasource_menu_reload_samples"))


class start_COOT(QtCore.QThread):
    def __init__(self, settings, interface):
        QtCore.QThread.__init__(self)
        self.settings = settings
        if interface == "test":
            self.pylib = "XChemCoot.py"
        elif interface == "new":
            self.pylib = "XChemCootNew.py"
        elif interface == "panddaV1":
            self.pylib = "XChemCootOld.py"
        elif interface == "reference":
            self.pylib = "XChemCootReference.py"
        elif interface == "buster":
            self.pylib = "XChemCootBuster.py"
        elif interface == "dimple_twin":
            self.pylib = "XChemCootTwin.py"

    def run(self):
        # coot at Diamond always or sometimes at least open in home directory
        # so then it won't find the .pkl file
        pickle.dump(
            self.settings,
            open(os.path.join(os.getenv("HOME"), ".xce_settings.pkl"), "wb"),
        )
        os.system(
            "cd {0!s}\ncoot --no-guano --no-state-script --script {1!s}".format(
                os.getenv("HOME"),
                os.path.join(os.getenv("XChemExplorer_DIR"), "coot", self.pylib),
            )
        )


class start_pandda_inspect(QtCore.QThread):
    def __init__(self, settings, xce_logfile):
        QtCore.QThread.__init__(self)
        self.panddas_directory = settings["panddas_directory"]
        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)

    def run(self):
        Cmds = (
            "#!" + os.getenv("SHELL") + "\n"
            "unset PYTHONPATH\n"
            "module load buster/20240123\n"
            "source /dls/science/groups/i04-1/software/pandda_0.2.12"
            "/ccp4/ccp4-7.0/bin/ccp4.setup-sh\n"
            "cd " + self.panddas_directory + "\n"
            "pandda.inspect\n"
        )

        self.Logfile.insert(
            "starting pandda.inspect with the following command:\n" + Cmds
        )
        os.system(Cmds)

class start_pandda_2_inspect(QtCore.QThread):
    def __init__(self, settings, xce_logfile):
        QtCore.QThread.__init__(self)
        self.panddas_directory = settings["panddas_directory"]
        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)

    def run(self):

        # If running at DLS get the appropriate path, otherwise rely on the user
        # having pandda2.inspect in path
        # Maybe XCE has a DLS settings file which would be a much nicer
        # solution to this? start_pandda_inspect calls the module system so this
        # doesn't seem any less portable
        dls_pandda_2_inspect_path = '/dls_sw/i04-1/software/PanDDA2Inspect/out/Moorhen-linux-x64/Moorhen'

        if os.path.exists(dls_pandda_2_inspect_path):
            Cmds = (
               "#!" + os.getenv("SHELL") + "\n"
               "" + str(dls_pandda_2_inspect_path) + " " + self.panddas_directory + "\n"
            )
        else:
            Cmds = (
                "#!" + os.getenv("SHELL") + "\n"
                "pandda2.inspect " + self.panddas_directory + "\n"
            )

        self.Logfile.insert(
            "starting pandda2.inspect with the following command:\n" + Cmds
        )
        os.system(Cmds)


# --- new module from hell -------------------------------------------------------------


class read_pinIDs_from_gda_logs(QtCore.QThread):
    def __init__(self, beamline, visit, database, gdaLogInstructions, xce_logfile):
        QtCore.QThread.__init__(self)
        self.beamline = beamline
        self.visit = visit

        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)

        self.db = XChemDB.data_source(database)
        self.allSamples = self.db.collected_xtals_during_visit_for_scoring(visit)

        self.gdaLogInstructions = gdaLogInstructions
        self.gda_log_start_line = gdaLogInstructions[0]
        self.gzipped_logs_parsed = gdaLogInstructions[1]

    def run(self):
        if self.gzipped_logs_parsed:
            self.Logfile.insert("parsed gzipped gda logfiles before, won't do again...")
        self.Logfile.insert(
            "will start parsing of current gda logfile at line {0!s}".format(
                str(self.gda_log_start_line)
            )
        )
        self.emit(
            QtCore.SIGNAL("update_status_bar(QString)"),
            "checking GDA logiles for pinID details",
        )
        pinDict, self.gda_log_start_line = XChemMain.get_gda_barcodes(
            self.allSamples,
            self.gzipped_logs_parsed,
            self.gda_log_start_line,
            self.beamline,
            self.xce_logfile,
        )

        self.update_database(pinDict)

        self.gdaLogInstructions = [self.gda_log_start_line, True]
        self.Logfile.insert("====== finished checking GDA logfiles ======")
        self.emit(
            QtCore.SIGNAL("update_gdaLog_parsing_instructions_and_score"),
            self.gdaLogInstructions,
        )
        self.emit(QtCore.SIGNAL("finished()"))

    def update_database(self, pinDict):
        self.Logfile.insert("updating database with pinDIs from GDA logfiles")

        progress = 0
        progress_step = XChemMain.getProgressSteps(len(pinDict))

        for sample in pinDict:
            self.emit(
                QtCore.SIGNAL("update_status_bar(QString)"),
                "updating pinID in DB for " + sample,
            )
            dbDict = {}
            dbDict["DataCollectionPinBarcode"] = pinDict[sample]
            self.db.update_specified_table(sample, dbDict, "collectionTable")
            progress += progress_step
            self.emit(QtCore.SIGNAL("update_progress_bar"), progress)


class choose_autoprocessing_outcome(QtCore.QThread):
    def __init__(
        self,
        database,
        visit,
        reference_file_list,
        preferences,
        projectDir,
        rescore,
        xce_logfile,
        agamemnon,
    ):
        QtCore.QThread.__init__(self)
        self.visit = visit
        self.projectDir = projectDir
        self.reference_file_list = reference_file_list
        self.rescore = rescore
        self.selection_mechanism = preferences["dataset_selection_mechanism"]
        self.acceptable_unitcell_volume_difference = preferences[
            "allowed_unitcell_difference_percent"
        ]
        self.acceptable_low_resolution_limit_for_data = preferences[
            "acceptable_low_resolution_limit_for_data"
        ]
        self.acceptable_low_resolution_Rmerge = preferences[
            "acceptable_low_resolution_Rmerge"
        ]
        self.agamemnon = agamemnon
        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)

        self.db = XChemDB.data_source(os.path.join(database))
        if self.agamemnon:
            self.allSamples = []
            for v in self.visit:
                x = self.db.collected_xtals_during_visit_for_scoring(v)
                for e in x:
                    self.allSamples.append(e)
        else:
            self.allSamples = self.db.collected_xtals_during_visit_for_scoring(visit)

    def run(self):
        progress = 0
        progress_step = XChemMain.getProgressSteps(len(self.allSamples))

        for sample in sorted(self.allSamples):
            if self.db.autoprocessing_result_user_assigned(sample) and not self.rescore:
                if os.path.isfile(
                    os.path.join(self.projectDir, sample, sample + ".mtz")
                ):
                    self.Logfile.warning(
                        "{0!s}: user has manually selected auto-processing result;"
                        " will NOT auto-select!".format(sample)
                    )
                    continue
                else:
                    self.Logfile.warning(
                        "{0!s}: user has manually selected auto-processing result"
                        " before, but {1!s}.mtz does not exist".format(sample, sample)
                    )
                    self.Logfile.insert("%s: selecting autoprocessing result" % sample)
            elif self.rescore:
                self.Logfile.warning(
                    "{0!s}: rescore selected -> might overwrite user selection!".format(
                        sample
                    )
                )
            else:
                self.Logfile.insert("%s: selecting autoprocessing result" % sample)
            dbList = self.db.all_autoprocessing_results_for_xtal_as_dict(sample)
            self.Logfile.insert(
                "%s: found %s different autoprocessing results"
                % (sample, str(len(dbList)))
            )

            # 0.) first check for which results files actually exist
            #
            dbList = self.checkExistingFiles(dbList, sample)
            if not dbList:
                self.Logfile.error(
                    sample + ": cannot find any MTZ & LOG files; skipping..."
                )
                continue

            # 1.) if posssible, only carry forward samples with similar UCvolume
            # and same point group
            dbList = self.selectResultsSimilarToReference(dbList)

            # 2.) if possible, only carry forward samples with low resolution Rmerge
            # smaller than specified in perferences
            dbList = self.selectResultsWithAcceptableLowResoRmerge(dbList)

            # 3.) Make selection based on speified selection mechanism
            if self.selection_mechanism == "IsigI*Comp*UniqueRefl":
                dbDict = self.selectHighestScore(dbList)
            elif self.selection_mechanism == "highest_resolution":
                dbDict = self.selectHighestResolution(dbList)
            elif "only" in self.selection_mechanism:
                dbDict = self.selectSpecificPipelineOnly(dbList)

            # 4.) Set new symbolic links in project directory
            XChemMain.linkAutoProcessingResult(
                sample, dbDict, self.projectDir, self.xce_logfile
            )

            # 5.) Determine DataProcessing Outcome
            dbDict["DataCollectionOutcome"] = self.determine_processing_outcome(dbDict)

            # 6.) update database
            dbDict["DataProcessingAutoAssigned"] = "True"
            self.updateDB(sample, dbDict)

            progress += progress_step
            self.emit(
                QtCore.SIGNAL("update_status_bar(QString)"),
                "scoring auto-processing results for " + sample,
            )
            self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        self.Logfile.insert("====== finished scoring data processing results ======")
        self.emit(QtCore.SIGNAL("populate_datasets_summary_table_NEW"))
        self.emit(QtCore.SIGNAL("finished()"))

    def report_forward_carried_pipelines(self, dbListOut, dbList):
        if not dbListOut:
            dbListOut = dbList
            self.Logfile.warning(
                "none of the MTZ files fulfilled criteria;"
                " will carry forward all results:"
            )
        else:
            self.Logfile.insert(
                "will carry forward the MTZ files from the following auto-processing"
                " pipelines:"
            )
        self.Logfile.insert(
            "{0:30} {1:10} {2:10} {3:10}".format(
                "pipeline", "Rmerge(Low)", "PG", "Score"
            )
        )
        self.Logfile.insert(
            "----------------------------------------------------------------------"
        )
        for resultDict in dbListOut:
            self.Logfile.insert(
                "{0:30} {1:10} {2:10} {3:10}".format(
                    resultDict["DataProcessingProgram"],
                    resultDict["DataProcessingRmergeLow"],
                    resultDict["DataProcessingPointGroup"],
                    resultDict["DataProcessingScore"],
                )
            )
        return dbListOut

    def checkExistingFiles(self, dbList, xtal):
        self.Logfile.insert("checking if MTZ & LOG files exisit")
        os.chdir(os.path.join(self.projectDir, xtal))
        self.Logfile.insert(
            xtal + ": changing directory " + os.path.join(self.projectDir, xtal)
        )
        dbListOut = []
        for resultDict in dbList:
            try:
                run = resultDict["DataCollectionRun"]
                subDir = resultDict["DataCollectionSubdir"]
                if subDir != "":
                    procCode = "_" + subDir
                else:
                    procCode = ""
                visit = resultDict["DataCollectionVisit"]
                autoproc = resultDict["DataProcessingProgram"]
                mtzFileAbs = resultDict["DataProcessingPathToMTZfile"]
                mtzfileName = mtzFileAbs[mtzFileAbs.rfind("/") + 1 :]
                logFileAbs = resultDict["DataProcessingPathToLogfile"]
                logfileName = logFileAbs[logFileAbs.rfind("/") + 1 :]

                mtzfile = os.path.join(
                    "autoprocessing",
                    visit + "-" + run + autoproc + procCode,
                    mtzfileName,
                )
                logfile = os.path.join(
                    "autoprocessing",
                    visit + "-" + run + autoproc + procCode,
                    logfileName,
                )

                if os.path.isfile(mtzfile):
                    self.Logfile.insert(xtal + ": found " + mtzfile)
                    if os.path.isfile(logfile):
                        self.Logfile.insert(xtal + ": found " + logfile)
                        dbListOut.append(resultDict)
                    else:
                        self.Logfile.error(
                            xtal + ": cannot find " + logfile + " ; skipping..."
                        )
                else:
                    self.Logfile.error(
                        xtal + ": cannot find " + mtzfile + " ; skipping..."
                    )

            except ValueError:
                pass
        return dbListOut

    def selectResultsSimilarToReference(self, dbList):
        self.Logfile.insert(
            "checking if MTZ files are similar to supplied reference files"
        )
        dbListOut = []
        for resultDict in dbList:
            try:
                if isinstance(float(resultDict["DataProcessingUnitCellVolume"]), float):
                    self.Logfile.insert(
                        "checking unit cell volume difference and point group:"
                    )
                    for reference_file in self.reference_file_list:
                        if not reference_file[4] == 0:
                            self.Logfile.insert(
                                "unitcell volume reference:" + str(reference_file[4])
                            )
                            self.Logfile.insert(
                                "unitcell volume dataset:  "
                                + str(resultDict["DataProcessingUnitCellVolume"])
                            )
                            unitcell_difference = round(
                                (
                                    math.fabs(
                                        reference_file[4]
                                        - float(
                                            resultDict["DataProcessingUnitCellVolume"]
                                        )
                                    )
                                    / reference_file[4]
                                )
                                * 100,
                                1,
                            )
                            self.Logfile.insert(
                                resultDict["DataProcessingProgram"]
                                + ": "
                                + str(unitcell_difference)
                                + "% difference -> pg(ref): "
                                + reference_file[5]
                                + " -> pg(mtz): "
                                + resultDict["DataProcessingPointGroup"]
                            )
                            if (
                                unitcell_difference
                                < self.acceptable_unitcell_volume_difference
                                and reference_file[5]
                                == resultDict["DataProcessingPointGroup"]
                            ):
                                self.Logfile.insert(
                                    "=> passed -> mtz file has same point group"
                                    " as reference file and similar unit cell volume"
                                )
                                dbListOut.append(resultDict)
                            else:
                                self.Logfile.warning(
                                    "mtz file has different point group/"
                                    " unit cell volume as reference file"
                                )
            except ValueError:
                pass
        dbListOut = self.report_forward_carried_pipelines(dbListOut, dbList)
        return dbListOut

    def selectResultsWithAcceptableLowResoRmerge(self, dbList):
        self.Logfile.insert(
            "checking if MTZ files have acceptable low resolution Rmerge values"
            " (currently set to %s)" % str(self.acceptable_low_resolution_Rmerge)
        )
        dbListOut = []
        for resultDict in dbList:
            try:
                if (
                    float(resultDict["DataProcessingRmergeLow"])
                    < self.acceptable_low_resolution_Rmerge
                ):
                    self.Logfile.insert(
                        resultDict["DataProcessingProgram"]
                        + ": Rmerge(low) of MTZ file is below threshold: "
                        + str(resultDict["DataProcessingRmergeLow"])
                    )
                    dbListOut.append(resultDict)
                else:
                    self.Logfile.warning(
                        resultDict["DataProcessingProgram"]
                        + ": Rmerge(low) of MTZ file is ABOVE threshold: "
                        + str(resultDict["DataProcessingRmergeLow"])
                    )
            except ValueError:
                pass
        dbListOut = self.report_forward_carried_pipelines(dbListOut, dbList)
        return dbListOut

    def selectHighestScore(self, dbList):
        tmp = []
        for resultDict in dbList:
            try:
                tmp.append(float(resultDict["DataProcessingScore"]))
            except ValueError:
                tmp.append(0.0)
        highestScoreDict = dbList[tmp.index(max(tmp))]
        return highestScoreDict

    def selectHighestResolution(self, dbList):
        tmp = []
        for resultDict in dbList:
            try:
                tmp.append(float(resultDict["DataProcessingResolutionHigh"]))
            except ValueError:
                tmp.append(100.0)
        highestResoDict = dbList[tmp.index(min(tmp))]
        return highestResoDict

    def selectSpecificPipelineOnly(self, dbList):
        tmp = []
        self.Logfile.insert(
            "selecting datasets by auto-processing pipeline: "
            + self.selection_mechanism
        )
        for resultDict in dbList:
            if (
                self.selection_mechanism == "dials - only"
                and "dials" in resultDict["DataProcessingProgram"]
            ):
                tmp.append(resultDict)
            if (
                self.selection_mechanism == "xia2 3d - only"
                and "3d-" in resultDict["DataProcessingProgram"]
            ):
                tmp.append(resultDict)
            if (
                self.selection_mechanism == "xia2 3dii - only"
                and "3dii" in resultDict["DataProcessingProgram"]
            ):
                tmp.append(resultDict)
            if (
                self.selection_mechanism == "autoProc - only"
                and "autoPROC" in resultDict["DataProcessingProgram"]
                and "staraniso" not in resultDict["DataProcessingProgram"]
            ):
                tmp.append(resultDict)
            if (
                self.selection_mechanism == "autoProc_staraniso - only"
                and "autoPROC" in resultDict["DataProcessingProgram"]
                and "staraniso" in resultDict["DataProcessingProgram"]
            ):
                tmp.append(resultDict)
        if not tmp:
            tmp = dbList
        pipelineDict = self.selectHighestScore(tmp)
        return pipelineDict

    def determine_processing_outcome(self, db_dict):
        outcome = "Failed - unknown"
        try:
            if (
                float(db_dict["DataProcessingResolutionHigh"])
                < self.acceptable_low_resolution_limit_for_data
            ):
                outcome = "success"
            else:
                outcome = "Failed - low resolution"
        except ValueError:
            pass
        return outcome

    def updateDB(self, sample, dbDict):
        self.Logfile.insert("{0!s}: updating database".format(sample))
        self.db.update_insert_data_source(sample, dbDict)


class read_write_autoprocessing_results_from_to_disc(QtCore.QThread):
    """
    major changes:
    - pkl file is obsolete
    - results for every autoprocessing result is recorded in new DB table
    - crystal centring images are copied into project directory
    - beamline directory in project directory will not be used anymore
    - users need to actively select the actual data collection visit as Data Collection
      Directory in the settings tab (e.g. /dls/i04-1/data/2017/mx15433-50
    - results from fast_dp are not copied over and included in analysis

    - DB mainTable gets flag if user updated autoprocessing selection
    - checking of reprocessed files needs to be explicit
    - all dictionaries used to store information are retired
    - at the moment one can only review/ rescore crystals collected during the selected
      visit
    - parsing of pinIDs in gda logfiles is still missing
    """

    def __init__(
        self, processedDir, database, projectDir, xce_logfile, target, agamemnon
    ):
        QtCore.QThread.__init__(self)
        self.processedDir = processedDir
        self.visit, self.beamline = XChemMain.getVisitAndBeamline(self.processedDir)
        self.projectDir = projectDir
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.target = target
        self.agamemnon = agamemnon

        self.db = XChemDB.data_source(os.path.join(database))
        self.exisitingSamples = self.getExistingSamples()

        self.toParse = [
            [
                os.path.join("*"),
                os.path.join("LogFiles", "*aimless.log"),
                os.path.join("DataFiles", "*free.mtz"),
            ],
            [
                os.path.join("*"),
                os.path.join("LogFiles", "*merging-statistics.json"),
                os.path.join("DataFiles", "*free.mtz"),
            ],
            [
                os.path.join("multi-xia2", "*"),
                os.path.join("LogFiles", "*aimless.log"),
                os.path.join("DataFiles", "*free.mtz"),
            ],
            [os.path.join("autoPROC"), "*aimless.log", "*truncate-unique.mtz"],
            [
                os.path.join("autoPROC"),
                # staraniso_alldata-unique.table1 only available in tar archive
                "*summary.tar.gz",
                "*staraniso_alldata-unique.mtz",
            ],
            [os.path.join("autoPROC-*"), "*aimless.log", "*truncate-unique.mtz"],
            [
                os.path.join("autoPROC-*"),
                "*summary.tar.gz",
                "*staraniso_alldata-unique.mtz",
            ],
        ]

    def run(self):
        self.parse_file_system()

    def getExistingSamples(self):
        existingSamples = {}
        self.Logfile.insert("reading existing samples from collectionTable")
        allEntries = self.db.execute_statement(
            "select CrystalName,DataCollectionVisit,DataCollectionRun,"
            "DataProcessingProgram, DataCollectionSubdir from collectionTable where"
            ' DataCollectionOutcome = "success"'
        )
        for item in allEntries:
            if str(item[0]) not in existingSamples:
                existingSamples[str(item[0])] = []
                self.Logfile.insert("%s: adding %s" % (str(item[0]), str(item[1])))
                # visit-runautoproc-subdir
            existingSamples[str(item[0])].append(
                str(item[1]) + "-" + str(item[2]) + str(item[3]) + "-" + str(item[4])
            )
        return existingSamples

    def createSampleDir(self, xtal):
        if not os.path.isdir(os.path.join(self.projectDir, xtal)):
            os.mkdir(os.path.join(self.projectDir, xtal))

    def createAutoprocessingDir(self, xtal, run, autoproc, proc_code):
        # create all the directories if necessary
        self.Logfile.insert(
            "%s: checking if new directory needs to be created for %s"
            % (
                xtal,
                os.path.join(
                    self.projectDir,
                    xtal,
                    "autoprocessing",
                    self.visit + "-" + run + autoproc + "_" + proc_code,
                ),
            )
        )
        if not os.path.isdir(os.path.join(self.projectDir, xtal, "autoprocessing")):
            os.mkdir(os.path.join(self.projectDir, xtal, "autoprocessing"))
        if not os.path.isdir(
            os.path.join(
                self.projectDir,
                xtal,
                "autoprocessing",
                self.visit + "-" + run + autoproc + "_" + proc_code,
            )
        ):
            self.Logfile.insert(
                "%s: making directory %s"
                % (
                    xtal,
                    os.path.join(
                        self.projectDir,
                        xtal,
                        "autoprocessing",
                        self.visit + "-" + run + autoproc + "_" + proc_code,
                    ),
                )
            )
            os.mkdir(
                os.path.join(
                    self.projectDir,
                    xtal,
                    "autoprocessing",
                    self.visit + "-" + run + autoproc + "_" + proc_code,
                )
            )
        else:
            self.Logfile.warning("%s: directory exists; skipping..." % xtal)

    def cleanUpDir(self, xtal, run, autoproc, mtzfile, logfile, proc_code):
        toKeep = [
            "staraniso_alldata-unique.mtz",
            "staraniso_alldata-unique.table1",
            "staraniso_alldata.log",
            xtal + ".mtz",
            xtal + ".log",
        ]
        os.chdir(
            os.path.join(
                self.projectDir,
                xtal,
                "autoprocessing",
                self.visit + "-" + run + autoproc + "_" + proc_code,
            )
        )
        for files in glob.glob("*"):
            if files not in toKeep:
                os.system("/bin/rm -f " + files)

    def copyMTZandLOGfiles(self, xtal, run, autoproc, mtzfile, logfile, proc_code):
        mtzNew = ""
        logNew = ""
        os.chdir(
            os.path.join(
                self.projectDir,
                xtal,
                "autoprocessing",
                self.visit + "-" + run + autoproc + "_" + proc_code,
            )
        )
        # MTZ file
        if not os.path.isfile(mtzfile[mtzfile.rfind("/") + 1 :]):
            self.Logfile.insert("%s: copying %s" % (xtal, mtzfile))
            os.system("/bin/cp " + mtzfile + " .")
            for mmcif in glob.glob(os.path.join(mtzfile[: mtzfile.rfind("/")], "*")):
                if mmcif.endswith(".mmcif"):
                    self.Logfile.insert("%s: copying %s" % (xtal, mmcif))
                    os.system("/bin/cp " + mmcif + " .")
                elif mmcif.endswith(".mmcif.bz2"):
                    self.Logfile.insert(
                        "%s: copying and decompressing %s" % (xtal, mmcif)
                    )
                    os.system("/bin/cp " + mmcif + " .")
                    os.system("bzip2 -d ./" + mmcif)
        if os.path.isfile(mtzfile[mtzfile.rfind("/") + 1 :]) and not os.path.isfile(
            xtal + ".mtz"
        ):
            os.symlink(mtzfile[mtzfile.rfind("/") + 1 :], xtal + ".mtz")
        if os.path.isfile(mtzfile[mtzfile.rfind("/") + 1 :]):
            mtzNew = os.path.join(
                self.projectDir,
                xtal,
                "autoprocessing",
                self.visit + "-" + run + autoproc + "_" + proc_code,
                mtzfile[mtzfile.rfind("/") + 1 :],
            )
        # MTZ file unmerged
        if os.path.isfile(mtzfile.replace("_free.mtz", "_scaled_unmerged.mtz")):
            self.Logfile.insert("%s: found unmerged mtz file" % xtal)
            if not os.path.isfile(xtal + "_unmerged.mtz"):
                self.Logfile.insert(
                    "%s: copying %s"
                    % (xtal, mtzfile.replace("_free.mtz", "_scaled_unmerged.mtz"))
                )
                os.system(
                    "/bin/cp "
                    + mtzfile.replace("_free.mtz", "_scaled_unmerged.mtz")
                    + " ."
                )
        # AIMLESS logfile
        if not os.path.isfile(logfile[logfile.rfind("/") + 1 :]):
            self.Logfile.insert("%s: copying %s" % (xtal, logfile))
            os.system("/bin/cp " + logfile + " .")
            if logfile.endswith("summary.tar.gz"):
                self.Logfile.insert("unpacking summary.tar.gz")
                os.system("tar -xzvf summary.tar.gz")
                logfile = logfile.replace(
                    "summary.tar.gz", "staraniso_alldata-unique.table1"
                )
                self.cleanUpDir(xtal, run, autoproc, mtzfile, logfile, proc_code)
        if os.path.isfile(logfile[logfile.rfind("/") + 1 :]) and not os.path.isfile(
            xtal + ".log"
        ):
            os.symlink(logfile[logfile.rfind("/") + 1 :], xtal + ".log")
        if os.path.isfile(logfile[logfile.rfind("/") + 1 :]):
            logNew = os.path.join(
                self.projectDir,
                xtal,
                "autoprocessing",
                self.visit + "-" + run + autoproc + "_" + proc_code,
                logfile[logfile.rfind("/") + 1 :],
            )
        # September 2021: xia2 does not also output a xia2.mmcif and json merging
        # statistics file, even if aimless was used for scaling, however,
        # the xia2.mmcif file is differently formatted than the one from dials
        # hence, if xtal.log file already exists, use this one
        if "aimless" in os.readlink(xtal + ".log"):
            logNew = os.path.join(
                self.projectDir,
                xtal,
                "autoprocessing",
                self.visit + "-" + run + autoproc + "_" + proc_code,
                xtal + ".log",
            )
        return mtzNew, logNew

    def makeJPGdir(self, xtal, run):
        if not os.path.isdir(os.path.join(self.projectDir, xtal, "jpg")):
            self.Logfile.insert("making jpg directory in " + xtal)
            os.mkdir(os.path.join(self.projectDir, xtal, "jpg"))
        if not os.path.isdir(
            os.path.join(self.projectDir, xtal, "jpg", self.visit + "-" + run)
        ):
            os.mkdir(os.path.join(self.projectDir, xtal, "jpg", self.visit + "-" + run))

    def copyJPGs(self, xtal, run, auto):
        self.Logfile.insert("%s: trying to copy crystal snapshots..." % xtal)
        found = False
        if self.agamemnon:
            proposal = self.visit.split("-")[0]
            self.Logfile.insert(
                "looking for images in "
                + os.path.join(
                    self.processedDir.replace(proposal, self.visit),
                    "jpegs",
                    auto,
                    self.target,
                    xtal,
                    run + "*.0.png",
                )
            )
            for img in glob.glob(
                os.path.join(
                    self.processedDir.replace(proposal, self.visit),
                    "jpegs",
                    auto,
                    self.target,
                    xtal,
                    run + "*.0.png",
                )
            ):
                if not os.path.isfile(
                    os.path.join(
                        self.projectDir,
                        xtal,
                        "jpg",
                        self.visit + "-" + run,
                        img[img.rfind("/") + 1 :],
                    )
                ):
                    self.Logfile.insert("%s: copying %s" % (xtal, img))
                    os.system(
                        "/bin/cp %s %s"
                        % (
                            img,
                            os.path.join(
                                self.projectDir,
                                xtal,
                                "jpg",
                                self.visit + "-" + auto + "_" + run,
                            ),
                        )
                    )
        else:
            for img in glob.glob(
                os.path.join(
                    self.processedDir.replace("processed", "jpegs"), run + "*.0.png"
                )
            ):
                found = True
                if not os.path.isfile(
                    os.path.join(
                        self.projectDir,
                        xtal,
                        "jpg",
                        self.visit + "-" + run,
                        img[img.rfind("/") + 1 :],
                    )
                ):
                    self.Logfile.insert("%s: copying %s" % (xtal, img))
                    os.system(
                        "/bin/cp %s %s"
                        % (
                            img,
                            os.path.join(
                                self.projectDir, xtal, "jpg", self.visit + "-" + run
                            ),
                        )
                    )
            if not found:
                for img in glob.glob(
                    os.path.join(
                        self.processedDir.replace("processed", "jpegs"),
                        xtal,
                        run + "*.0.png",
                    )
                ):
                    if not os.path.isfile(
                        os.path.join(
                            self.projectDir,
                            xtal,
                            "jpg",
                            self.visit + "-" + run,
                            img[img.rfind("/") + 1 :],
                        )
                    ):
                        self.Logfile.insert("%s: copying %s" % (xtal, img))
                        os.system(
                            "/bin/cp %s %s"
                            % (
                                img,
                                os.path.join(
                                    self.projectDir, xtal, "jpg", self.visit + "-" + run
                                ),
                            )
                        )

    def findJPGs(self, xtal, run):
        jpgDict = {}
        for n, img in enumerate(
            glob.glob(
                os.path.join(
                    self.projectDir, xtal, "jpg", self.visit + "-" + run, "*.0.png"
                )
            )
        ):
            if n <= 3:
                jpgDict["DataCollectionCrystalImage" + str(n + 1)] = img
        return jpgDict

    def readProcessingUpdateResults(
        self, xtal, folder, log, mtz, timestamp, current_run, autoproc, proc_code
    ):
        db_dict = {}
        self.Logfile.warning("%s: looking for %s" % (xtal, os.path.join(folder, mtz)))
        for mtzfile in glob.glob(os.path.join(folder, mtz)):
            self.Logfile.insert("%s: found %s" % (xtal, mtzfile))
            self.Logfile.warning(
                "%s: looking for %s" % (xtal, os.path.join(folder, log))
            )
            for logfile in glob.glob(os.path.join(folder, log)):
                self.Logfile.insert("%s: found %s" % (xtal, logfile))
                self.createAutoprocessingDir(xtal, current_run, autoproc, proc_code)
                mtzNew, logNew = self.copyMTZandLOGfiles(
                    xtal, current_run, autoproc, mtzfile, logfile, proc_code
                )
                if self.target == "=== project directory ===":
                    target = "unknown"
                else:
                    target = self.target
                db_dict = {
                    "DataCollectionDate": timestamp,
                    "DataProcessingPathToLogfile": logNew,
                    "DataProcessingPathToMTZfile": mtzNew,
                    "DataProcessingDirectoryOriginal": folder,
                    # success in collection Table only means that a logfile was found
                    "DataCollectionOutcome": "success",
                    "ProteinName": target,
                }
                db_dict.update(XChemUtils.parse().read_aimless_logfile(logNew))
                # image exist even if data processing failed
                db_dict.update(self.findJPGs(xtal, current_run))
                db_dict["DataCollectionBeamline"] = self.beamline
                self.update_data_collection_table(
                    xtal, current_run, autoproc, db_dict, proc_code
                )

    def getAutoProc(self, folder_rel, staraniso):
        self.Logfile.insert("checking name of auto-processing pipeline...")
        folder = os.path.realpath(folder_rel)
        self.Logfile.insert("folder: " + folder)
        autoproc = "unknown"
        if "ap-run" in folder:
            autoproc = "autoPROC"
        else:
            for f in folder.split("/"):
                if "autoPROC" in f:
                    autoproc = f + staraniso
                    break
                elif "xia2" in f:
                    autoproc = f
                    break
                elif "dials" in f:
                    autoproc = f
                    break
        self.Logfile.insert("name of auto-processing pipeline: %s" % autoproc)
        return autoproc

    def update_data_collection_table(
        self, xtal, current_run, autoproc, db_dict, proc_code
    ):
        condition_dict = {
            "CrystalName": xtal,
            "DataCollectionVisit": self.visit,
            "DataCollectionRun": current_run,
            "DataCollectionSubdir": proc_code,
            "DataProcessingProgram": autoproc,
        }
        self.db.update_insert_any_table("collectionTable", db_dict, condition_dict)

    def alreadyParsed(self, xtal, current_run, proc_code, autoproc):
        parsed = False
        self.Logfile.insert("checking if this processed directory was already parsed")
        if xtal in self.exisitingSamples:
            self.Logfile.insert(xtal + " exists in collectionTable")
            self.Logfile.insert("current identifier:")
            self.Logfile.insert(
                self.visit + "-" + current_run + autoproc + "-" + proc_code
            )
            self.Logfile.insert("indentifier in collectionTable")
            for x in self.exisitingSamples[xtal]:
                self.Logfile.insert(x)
            # visit-runautoproc-subdir
            if (
                self.visit + "-" + current_run + autoproc + "-" + proc_code
                in self.exisitingSamples[xtal]
            ):
                self.Logfile.warning(
                    "%s: results from %s already parsed; skipping..."
                    % (xtal, self.visit + "-" + current_run + autoproc)
                )
                parsed = True
        return parsed

    def empty_folder(self, xtal, folder):
        empty = True
        stuff = []
        for x in glob.glob(os.path.join(folder, "*")):
            stuff.append(x)
        if not stuff:
            self.Logfile.warning(
                "{0!s}: {1!s} is empty; probably waiting for autoprocessing to finish;"
                " try later!".format(xtal, folder)
            )
        else:
            empty = False
        return empty

    def junk(self, folder):
        do_not_parse = False
        if not os.path.isdir(folder):
            do_not_parse = True
        elif "dimple" in folder:
            do_not_parse = True
        return do_not_parse

    def parse_file_system(self):
        if self.agamemnon:
            d = self.processedDir[: self.processedDir.rfind("/")]
            t = self.processedDir[self.processedDir.rfind("/") + 1 :]
            autoDir = os.path.join(d, "auto", t)
        else:
            autoDir = self.processedDir

        self.Logfile.insert("checking for new data processing results in " + autoDir)
        progress = 0
        progress_step = XChemMain.getProgressSteps(
            len(glob.glob(os.path.join(autoDir, "*")))
        )

        runList = []
        self.Logfile.insert("--> " + os.path.join(autoDir, "*"))
        for nx, collected_xtals in enumerate(
            sorted(glob.glob(os.path.join(autoDir, "*")))
        ):
            self.Logfile.insert("%s: %s" % (nx, collected_xtals))
            self.visit = collected_xtals.split("/")[5]
            if (
                "tmp" in collected_xtals
                or "results" in collected_xtals
                or "scre" in collected_xtals
            ):
                continue
            if not os.path.isdir(collected_xtals):
                continue

            if (
                "tmp" in collected_xtals
                or "results" in collected_xtals
                or "scre" in collected_xtals
            ):
                continue
            if not os.path.isdir(collected_xtals):
                self.Logfile.warning(collected_xtals + " is not a directory")
                continue

            xtal = collected_xtals[collected_xtals.rfind("/") + 1 :]

            self.Logfile.insert("%s: checking auto-processing results" % xtal)
            self.createSampleDir(xtal)

            if self.target == "=== project directory ===":
                runDir = os.path.join(collected_xtals, "processed", "*")
            else:
                runDir = os.path.join(collected_xtals, "*")
            self.Logfile.insert("current runDir: " + runDir)

            for run in sorted(glob.glob(runDir)):
                current_run = run[run.rfind("/") + 1 :]
                for code in glob.glob(os.path.join(run, "*")):
                    if os.path.islink(code):
                        continue
                    proc_code = code.split("/")[len(code.split("/")) - 1]
                    self.Logfile.insert(xtal + ": processed directory -> " + proc_code)
                    if current_run + proc_code in runList:
                        continue
                    self.Logfile.insert(
                        "%s -> run: %s -> current run: %s -> %s"
                        % (xtal, run, current_run, proc_code)
                    )
                    timestamp = datetime.fromtimestamp(os.path.getmtime(run)).strftime(
                        "%Y-%m-%d %H:%M:%S"
                    )
                    # create directory for crystal aligment images in projectDir
                    self.makeJPGdir(xtal, current_run)
                    # 'non-auto' is irrelevant here
                    self.copyJPGs(xtal, current_run, "non-auto")

                    for item in self.toParse:
                        procDir = os.path.join(code, item[0])
                        logfile = item[1]
                        mtzfile = item[2]
                        self.Logfile.insert(
                            "%s: search template: procDir - logfile - mtzfile" % xtal
                        )
                        self.Logfile.insert("%s: procDir = %s" % (xtal, procDir))
                        self.Logfile.insert("%s: logfile = %s" % (xtal, logfile))
                        self.Logfile.insert("%s: mtzfile = %s" % (xtal, mtzfile))

                        for folder in glob.glob(procDir):
                            staraniso = ""
                            self.Logfile.insert("%s: searching %s" % (xtal, folder))
                            if self.junk(folder):
                                continue
                            if self.empty_folder(xtal, folder):
                                continue
                            if "staraniso" in logfile or "summary.tar.gz" in logfile:
                                staraniso = "_staraniso"
                            autoproc = self.getAutoProc(folder, staraniso)
                            if self.alreadyParsed(
                                xtal, current_run, proc_code, autoproc
                            ):
                                continue
                            self.readProcessingUpdateResults(
                                xtal,
                                folder,
                                logfile,
                                mtzfile,
                                timestamp,
                                current_run,
                                autoproc,
                                proc_code,
                            )
                    runList.append(current_run + proc_code)
            progress += progress_step
            self.emit(
                QtCore.SIGNAL("update_status_bar(QString)"),
                "parsing auto-processing results for " + xtal,
            )
            self.emit(QtCore.SIGNAL("update_progress_bar"), progress)

        self.Logfile.insert("====== finished parsing beamline directory ======")
        self.emit(QtCore.SIGNAL("read_pinIDs_from_gda_logs"))
        self.emit(QtCore.SIGNAL("finished()"))
