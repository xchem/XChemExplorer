import XChemLog
import XChemDB
import os
import glob
import subprocess
import getpass
import gzip
from datetime import datetime


def space_group_list():
    space_group_list = [
        "P1",
        "P2",
        "P21",
        "C121",
        "P1211",
        "P121",
        "I2",
        "I121",
        "P222",
        "P2122",
        "P2212",
        "P2221",
        "P21212",
        "P21221",
        "P22121",
        "P212121",
        "C222",
        "C2221",
        "F222",
        "I222",
        "I212121",
        "P4",
        "P41",
        "P42",
        "P43",
        "I4",
        "I41",
        "P422",
        "P4212",
        "P4122",
        "P41212",
        "P4222",
        "P42212",
        "P4322",
        "P43212",
        "I422",
        "I4122",
        "P3",
        "P31",
        "P32",
        "P312",
        "P321",
        "P3112",
        "P3121",
        "P3212",
        "P3221",
        "P6",
        "P61",
        "P65",
        "P62",
        "P64",
        "P63",
        "P622",
        "P6122",
        "P6522",
        "P6222",
        "P6422",
        "P6322",
        "H3",
        "H32",
        "P23",
        "F23",
        "I23",
        "P213",
        "I213",
        "P432",
        "P4232",
        "F432",
        "F4132",
        "I432",
        "P4332",
        "P4132",
        "I4132",
    ]
    return space_group_list


def get_target_and_visit_list(beamline_directory, agamemnon):
    target_list = ["=== SELECT TARGET ===", "=== project directory ==="]
    visit_list = []
    # the beamline directory could be a the real directory or
    # a directory where the visits are linked into
    if (
        len(beamline_directory.split("/"))
        and beamline_directory.split("/")[1] == "dls"
        and beamline_directory.split("/")[3] == "data"
        and "labxchem" not in beamline_directory
    ):
        visit_list.append(beamline_directory)
    else:
        visit_list.append(os.path.realpath(beamline_directory))

    for visit in visit_list:
        print("-->", os.path.join(visit, "processed", "*"))
        if agamemnon:
            for target in glob.glob(os.path.join(visit, "processed", "auto", "*")):
                print(target)
                if target[target.rfind("/") + 1 :] not in [
                    "results",
                    "README-log",
                    "edna-latest.html",
                ]:
                    if target[target.rfind("/") + 1 :] not in target_list:
                        target_list.append(target[target.rfind("/") + 1 :])
        else:
            for target in glob.glob(os.path.join(visit, "processed", "*")):
                print(target)
                if target[target.rfind("/") + 1 :] not in [
                    "results",
                    "README-log",
                    "edna-latest.html",
                ]:
                    if target[target.rfind("/") + 1 :] not in target_list:
                        target_list.append(target[target.rfind("/") + 1 :])
    return target_list, visit_list


def get_jobs_running_on_cluster():
    out_dict = {}

    dimple_jobs = []
    acedrg_jobs = []
    pandda_jobs = []
    refmac_jobs = []
    xia2_jobs = []
    others_jobs = []

    # note: each job_details list contains a list with
    # [job_ID, status, run_time]
    out = subprocess.Popen(["qstat"], stdout=subprocess.PIPE)
    for n, line in enumerate(iter(out.stdout.readline, "")):
        if len(line.split()) >= 7:
            if line.split()[3] == getpass.getuser():

                job_id = line.split()[0]
                job_name = line.split()[2]
                job_status = line.split()[4]

                ##########################################################
                # determine run time of each job in minutes
                start_date = ""
                start_time = ""
                run_time_minutes = ""
                start_date = line.split()[5]
                if len(start_date.split("/")) == 3:
                    month_start = start_date.split("/")[0]
                    day_start = start_date.split("/")[1]
                    year_start = start_date.split("/")[2]

                start_time = line.split()[6]
                if len(start_time.split(":")) == 3:
                    hour_start = start_time.split(":")[0]
                    minute_start = start_time.split(":")[1]
                    second_start = start_time.split(":")[2]

                if start_time != "" and start_date != "":
                    start = "{0!s}-{1!s}-{2!s} {3!s}:{4!s}:{5!s}".format(
                        year_start,
                        month_start,
                        day_start,
                        hour_start,
                        minute_start,
                        second_start,
                    )
                    run_time = datetime.now() - datetime.strptime(
                        start, "%Y-%m-%d %H:%M:%S"
                    )
                    run_time_minutes = int(run_time.total_seconds() / 60)

                ##########################################################
                # determine run time of each job in minutes
                if "dimple" in job_name:
                    dimple_jobs.append([job_id, job_status, run_time_minutes])
                elif "acedrg" in job_name:
                    acedrg_jobs.append([job_id, job_status, run_time_minutes])
                elif "pandda" in job_name:
                    pandda_jobs.append([job_id, job_status, run_time_minutes])
                elif "refmac" in job_name:
                    refmac_jobs.append([job_id, job_status, run_time_minutes])
                elif "xia2" in job_name:
                    xia2_jobs.append([job_id, job_status, run_time_minutes])
                else:
                    others_jobs.append([job_id, job_status, run_time_minutes])

    out_dict["dimple"] = dimple_jobs
    out_dict["acedrg"] = acedrg_jobs
    out_dict["pandda"] = pandda_jobs
    out_dict["refmac"] = refmac_jobs
    out_dict["xia2"] = xia2_jobs
    out_dict["others"] = others_jobs

    return out_dict


def print_acedrg_status(xce_logfile, xtal_db_dict):
    Logfile = XChemLog.updateLog(xce_logfile)
    Logfile.insert("compound restraints summary:")
    pending = 0
    started = 0
    running = 0
    missing_smiles = 0
    failed = 0
    success = 0
    unknown = 0
    for xtal in xtal_db_dict:
        db_dict = xtal_db_dict[xtal]
        status = db_dict["RefinementCIFStatus"]
        if "pending" in status:
            pending += 1
        elif "started" in status:
            started += 1
        elif "running" in status:
            running += 1
        elif "missing" in status:
            missing_smiles += 1
        elif "failed" in status:
            failed += 1
        elif "generated" in status:
            success += 1
        else:
            unknown += 1
    Logfile.insert("restraint generation pending: ...... {0!s}".format(str(pending)))
    Logfile.insert("restraint generation started: ...... {0!s}".format(str(started)))
    Logfile.insert("restraint generation running: ...... {0!s}".format(str(running)))
    Logfile.insert(
        "missing smiles string: ............. {0!s}".format(str(missing_smiles))
    )
    Logfile.insert("restraint generation failed: ....... {0!s}".format(str(failed)))
    Logfile.insert("restraints successfully created: ... {0!s}".format(str(success)))
    Logfile.insert("unknown status: .................... {0!s}".format(str(unknown)))


def print_cluster_status_message(program, cluster_dict, xce_logfile):
    Logfile = XChemLog.updateLog(xce_logfile)
    Logfile.insert("cluster status summary:")
    Logfile.insert(
        "{0!s} {1!s} jobs are running on the cluster".format(
            len(cluster_dict[program]), program
        )
    )
    if len(cluster_dict[program]) > 0:
        cumulative_runtime = 0
        job_ids = []
        for n, item in enumerate(cluster_dict[program]):
            cumulative_runtime += item[2]
            if not item[0] in job_ids:
                job_ids.append(item[0])
        average_runtime = round(float(cumulative_runtime) / float(n + 1), 0)
        Logfile.insert("average run time " + str(average_runtime) + " minutes")
        if job_ids:
            Logfile.insert(
                "you can kill them by pasting the following line"
                " into a new terminal window:"
            )
            out = "qdel "
            for job in job_ids:
                out += str(job) + " "
            Logfile.insert(out)


def get_datasource_summary(db_file):
    db = XChemDB.data_source(db_file)

    out_dict = {}

    out_dict["nr_samples"] = len(
        db.execute_statement(
            "select CrystalName from mainTable where CrystalName is not NULL;"
        )
    )
    out_dict["nr_samples_failed_to_mount"] = len(
        db.execute_statement(
            "select HarvestStatus from mainTable where HarvestStatus is 'fail';"
        )
    )

    out_dict["nr_smiles_for_samples"] = len(
        db.execute_statement(
            "select compoundSMILES from mainTable"
            " where compoundSMILES is not (NULL or '')"
        )
    )

    out_dict["nr_data_collection_success"] = len(
        db.execute_statement(
            "select DataCollectionOutcome from mainTable"
            " where DataCollectionOutcome is 'success';"
        )
    )
    out_dict["nr_data_collection_centring_fail"] = len(
        db.execute_statement(
            "select DataCollectionOutcome from mainTable"
            " where DataCollectionOutcome is 'Failed - centring failed';"
        )
    )
    out_dict["nr_data_collection_no-diffraction"] = len(
        db.execute_statement(
            "select DataCollectionOutcome from mainTable"
            " where DataCollectionOutcome is 'Failed - no diffraction';"
        )
    )
    out_dict["nr_data_collection_processing_fail"] = len(
        db.execute_statement(
            "select DataCollectionOutcome from mainTable"
            " where DataCollectionOutcome is 'Failed - processing';"
        )
    )
    out_dict["nr_data_collection_loop-empty"] = len(
        db.execute_statement(
            "select DataCollectionOutcome from mainTable"
            " where DataCollectionOutcome is 'Failed - loop empty';"
        )
    )
    out_dict["nr_data_collection_loop-broken"] = len(
        db.execute_statement(
            "select DataCollectionOutcome from mainTable where"
            " DataCollectionOutcome is 'Failed - loop broken';"
        )
    )
    out_dict["nr_data_collection_low-resolution"] = len(
        db.execute_statement(
            "select DataCollectionOutcome from mainTable where"
            " DataCollectionOutcome is 'Failed - low resolution';"
        )
    )
    out_dict["nr_data_collection_no-X-rays"] = len(
        db.execute_statement(
            "select DataCollectionOutcome from mainTable where"
            " DataCollectionOutcome is 'Failed - no X-rays';"
        )
    )
    out_dict["nr_data_collection_unknown"] = len(
        db.execute_statement(
            "select DataCollectionOutcome from mainTable where"
            " DataCollectionOutcome is 'Failed - unknown';"
        )
    )

    out_dict["nr_data_collection_failed"] = (
        out_dict["nr_data_collection_centring_fail"]
        + out_dict["nr_data_collection_no-diffraction"]
        + out_dict["nr_data_collection_processing_fail"]
        + out_dict["nr_data_collection_loop-empty"]
        + out_dict["nr_data_collection_loop-broken"]
        + out_dict["nr_data_collection_low-resolution"]
        + out_dict["nr_data_collection_no-X-rays"]
        + out_dict["nr_data_collection_unknown"]
    )

    out_dict["nr_data_collection_pending"] = (
        out_dict["nr_samples"]
        - out_dict["nr_data_collection_success"]
        - out_dict["nr_data_collection_centring_fail"]
        - out_dict["nr_data_collection_no-diffraction"]
        - out_dict["nr_data_collection_processing_fail"]
        - out_dict["nr_data_collection_loop-empty"]
        - out_dict["nr_data_collection_loop-broken"]
        - out_dict["nr_data_collection_low-resolution"]
        - out_dict["nr_data_collection_no-X-rays"]
        - out_dict["nr_data_collection_unknown"]
    )

    out_dict["nr_initial_maps_available"] = len(
        db.execute_statement(
            "select DimplePathToMTZ from mainTable where DimplePathToMTZ is not '';"
        )
    )
    out_dict["nr_initial_maps_fail"] = len(
        db.execute_statement(
            "select DataProcessingDimpleSuccessful from mainTable"
            " where DataProcessingDimpleSuccessful = 'False';"
        )
    )
    out_dict["nr_initial_maps_pending"] = (
        out_dict["nr_data_collection_success"]
        - out_dict["nr_initial_maps_available"]
        - out_dict["nr_initial_maps_fail"]
    )

    out_dict["nr_pandda_hits"] = len(
        db.execute_statement(
            "select DimplePANDDAhit from mainTable where DimplePANDDAhit = 'True';"
        )
    )
    out_dict["nr_pandda_reject"] = len(
        db.execute_statement(
            "select DimplePANDDAreject from mainTable"
            " where DimplePANDDAreject = 'True';"
        )
    )
    out_dict["nr_pandda_processed"] = (
        len(
            db.execute_statement(
                "select DimplePANDDAwasRun from mainTable"
                " where DimplePANDDAwasRun = 'True';"
            )
        )
        - out_dict["nr_pandda_hits"]
        - out_dict["nr_pandda_reject"]
    )
    out_dict["nr_pandda_pending"] = (
        out_dict["nr_initial_maps_available"]
        - out_dict["nr_pandda_hits"]
        - out_dict["nr_pandda_reject"]
        - out_dict["nr_pandda_processed"]
    )

    out_dict["nr_cif_files"] = len(
        db.execute_statement(
            "select RefinementCIF from mainTable"
            " where RefinementCIF is not (Null or '');"
        )
    )

    out_dict["nr_analysis-pending"] = len(
        db.execute_statement(
            "select RefinementOutcome from mainTable"
            " where RefinementOutcome is '1 - Analysis Pending';"
        )
    )
    out_dict["nr_pandda-models"] = len(
        db.execute_statement(
            "select RefinementOutcome from mainTable"
            " where RefinementOutcome is '2 - PANDDA model';"
        )
    )
    out_dict["nr_in-refinement"] = len(
        db.execute_statement(
            "select RefinementOutcome from mainTable"
            " where RefinementOutcome is '3 - In Refinement';"
        )
    )
    out_dict["nr_comp-chem-ready"] = len(
        db.execute_statement(
            "select RefinementOutcome from mainTable"
            " where RefinementOutcome is '4 - ComChem ready';"
        )
    )
    out_dict["nr_deposition-ready"] = len(
        db.execute_statement(
            "select RefinementOutcome from mainTable"
            " where RefinementOutcome is '5 - Deposition ready';"
        )
    )

    return out_dict


def change_links_to_selected_data_collection_outcome(
    sample,
    data_collection_dict,
    data_collection_column_three_dict,
    dataset_outcome_dict,
    initial_model_directory,
    data_source_file,
    xce_logfile,
):
    Logfile = XChemLog.updateLog(xce_logfile)
    # find out which row was selected in respective data collection table
    selected_processing_result = "n/a"
    indexes = (
        data_collection_column_three_dict[sample][0].selectionModel().selectedRows()
    )
    if indexes:  # i.e. logfile exists
        for index in sorted(indexes):
            selected_processing_result = index.row()

    for n, entry in enumerate(data_collection_dict[sample]):
        if entry[0] == "logfile":
            if entry[7] == selected_processing_result:
                visit = entry[1]
                run = entry[2]
                autoproc = entry[4]
                db_dict = entry[6]
                path_to_logfile = db_dict["DataProcessingPathToLogfile"]
                path_to_mtzfile = db_dict["DataProcessingPathToMTZfile"]
                mtz_filename = db_dict["DataProcessingMTZfileName"]
                log_filename = db_dict["DataProcessingLOGfileName"]
                #                relative_path_to_mtzfile='./'+path_to_mtzfile.replace(initial_model_directory,'')
                relative_path_to_mtzfile = "./" + path_to_mtzfile.replace(
                    os.path.join(initial_model_directory, sample), ""
                )
                if relative_path_to_mtzfile.startswith(".//"):
                    relative_path_to_mtzfile = relative_path_to_mtzfile.replace(
                        ".//", "./"
                    )
                relative_path_to_logfile = "./" + path_to_logfile.replace(
                    os.path.join(initial_model_directory, sample), ""
                )
                if relative_path_to_logfile.startswith(".//"):
                    relative_path_to_logfile = relative_path_to_logfile.replace(
                        ".//", "./"
                    )

                # first check if folders and files exist
                # since user might do this before data are actually copied over

                if os.path.isdir(
                    os.path.join(
                        initial_model_directory,
                        sample,
                        "autoprocessing",
                        visit + "-" + run + autoproc,
                    )
                ):
                    db_dict["DataProcessingAutoAssigned"] = "False"
                    Logfile.insert(
                        "changing directory to: "
                        + os.path.join(initial_model_directory, sample)
                    )
                    os.chdir(os.path.join(initial_model_directory, sample))
                    # first remove old links
                    os.system("/bin/rm " + sample + ".mtz 2> /dev/null")
                    os.system("/bin/rm " + sample + ".log 2> /dev/null")
                    # make new links
                    Logfile.insert(
                        "setting relative symlink: "
                        + os.path.join(relative_path_to_logfile, log_filename)
                        + " -> "
                        + sample
                        + ".log"
                    )
                    os.symlink(
                        os.path.join(relative_path_to_logfile, log_filename),
                        sample + ".log",
                    )
                    Logfile.insert(
                        "setting relative symlink: "
                        + os.path.join(relative_path_to_mtzfile, mtz_filename)
                        + " -> "
                        + sample
                        + ".mtz"
                    )
                    os.symlink(
                        os.path.join(relative_path_to_mtzfile, mtz_filename),
                        sample + ".mtz",
                    )

                    # update data source
                    data_source = XChemDB.data_source(data_source_file)
                    data_source.update_insert_data_source(sample, db_dict)

                else:
                    Logfile.insert("please copy data to PROJECT DIRECTORY first!")


def get_gda_barcodes(
    sampleList, gzipped_logs_parsed, gda_log_start_line, beamline, xce_logfile
):
    Logfile = XChemLog.updateLog(xce_logfile)
    Logfile.insert(
        "checking GDA logfile in {0!s}".format(
            os.path.join("/dls_sw", beamline, "logs")
        )
    )
    pinDict = {}
    found_barcode_entry = False
    for gdaLogFile in glob.glob(
        os.path.join("/dls_sw", beamline, "logs", "gda-server*log*")
    ):
        if gdaLogFile.endswith("tmp"):
            Logfile.warning("ignoring temporary file " + gdaLogFile)
            continue
        Logfile.insert("parsing {0!s}".format(gdaLogFile))
        if gzipped_logs_parsed and gdaLogFile.endswith(".gz"):
            Logfile.insert(
                "{0!s} was already parsed during this visit".format(gdaLogFile)
            )
            continue
        if gdaLogFile.endswith(".gz"):
            try:
                with gzip.open(gdaLogFile, "r") as f:
                    for line in f:
                        if "BART SampleChanger - getBarcode() returning" in line:
                            barcode = line.split()[len(line.split()) - 1]
                            found_barcode_entry = True
                        if found_barcode_entry:
                            if "Snapshots will be saved" in line:
                                sampleID = line.split()[len(line.split()) - 1].split(
                                    "/"
                                )[-1]
                                if sampleID in sampleList:
                                    pinDict[sampleID] = barcode
                                    Logfile.insert(
                                        "found: sample={0!s}, barcode={1!s},"
                                        " file={2!s}".format(
                                            sampleID, barcode, gdaLogFile
                                        )
                                    )
                                found_barcode_entry = False
            except IOError:
                Logfile.warning("cannot open file %s" % gdaLogFile)
        else:
            try:
                for n, line in enumerate(
                    open(gdaLogFile).readlines()[gda_log_start_line:]
                ):
                    if "BART SampleChanger - getBarcode() returning" in line:
                        barcode = line.split()[len(line.split()) - 1]
                        found_barcode_entry = True
                    if found_barcode_entry:
                        if "Snapshots will be saved" in line:
                            sampleID = line.split()[len(line.split()) - 1].split("/")[
                                -1
                            ]
                            if sampleID in sampleList:
                                pinDict[sampleID] = barcode
                                Logfile.insert(
                                    "found: sample={0!s}, barcode={1!s},"
                                    " file={2!s}".format(sampleID, barcode, gdaLogFile)
                                )
                            found_barcode_entry = False
            except IOError:
                Logfile.error("IOError -> " + gdaLogFile)

            try:
                gda_log_start_line = gda_log_start_line + n - 1
            except UnboundLocalError:
                gda_log_start_line = gda_log_start_line

    return pinDict, gda_log_start_line


def linkAutoProcessingResult(xtal, dbDict, projectDir, xce_logfile):
    Logfile = XChemLog.updateLog(xce_logfile)

    run = dbDict["DataCollectionRun"]
    subDir = dbDict["DataCollectionSubdir"]
    if subDir != "":
        procCode = "_" + subDir
    else:
        procCode = ""
    visit = dbDict["DataCollectionVisit"]
    autoproc = dbDict["DataProcessingProgram"]
    mtzFileAbs = dbDict["DataProcessingPathToMTZfile"]
    mtzfile = mtzFileAbs[mtzFileAbs.rfind("/") + 1 :]
    logFileAbs = dbDict["DataProcessingPathToLogfile"]
    logfile = logFileAbs[logFileAbs.rfind("/") + 1 :]

    Logfile.insert("changing directory to " + os.path.join(projectDir, xtal))
    os.chdir(os.path.join(projectDir, xtal))

    # MTZ file
    Logfile.warning("removing %s.mtz" % xtal)
    os.system("/bin/rm %s.mtz" % xtal)
    Logfile.insert(
        xtal
        + ": looking for "
        + os.path.join(
            "autoprocessing", visit + "-" + run + autoproc + procCode, mtzfile
        )
    )
    if os.path.isfile(
        os.path.join("autoprocessing", visit + "-" + run + autoproc + procCode, mtzfile)
    ):
        os.symlink(
            os.path.join(
                "autoprocessing", visit + "-" + run + autoproc + procCode, mtzfile
            ),
            xtal + ".mtz",
        )
        Logfile.insert("linking MTZ file from different auto-processing pipeline:")
        Logfile.insert(
            "ln -s "
            + os.path.join(
                "autoprocessing", visit + "-" + run + autoproc + procCode, mtzfile
            )
            + " "
            + xtal
            + ".mtz"
        )
    # LOG file
    Logfile.warning("removing %s.log" % xtal)
    os.system("/bin/rm %s.log" % xtal)
    Logfile.insert(
        xtal
        + ": looking for "
        + os.path.join(
            "autoprocessing", visit + "-" + run + autoproc + procCode, logfile
        )
    )
    if os.path.isfile(
        os.path.join("autoprocessing", visit + "-" + run + autoproc + procCode, logfile)
    ):
        os.symlink(
            os.path.join(
                "autoprocessing", visit + "-" + run + autoproc + procCode, logfile
            ),
            xtal + ".log",
        )
        Logfile.insert("linking LOG file from different auto-processing pipeline:")
        Logfile.insert(
            "ln -s "
            + os.path.join(
                "autoprocessing", visit + "-" + run + autoproc + procCode, logfile
            )
            + " "
            + xtal
            + ".log"
        )


def getProgressSteps(iterations):
    if iterations == 0:
        progress_step = 1
    else:
        progress_step = 100 / float(iterations)
    return progress_step


def getVisitAndBeamline(visitDirectory):
    visit = "unknown"
    beamline = "unknown"
    if "attic" in visitDirectory:
        try:
            visit = visitDirectory.split("/")[6]
            beamline = visitDirectory.split("/")[3]
        except IndexError:
            pass
    else:
        try:
            visit = visitDirectory.split("/")[5]
            beamline = visitDirectory.split("/")[2]
        except IndexError:
            pass
    if not visitDirectory.startswith("/dls"):
        # this is all a bit of a fudge in case someone transfers a DLS visit directory
        # back home does most certainly not catch all possible scenarios
        if visitDirectory.split("/")[len(visitDirectory.split("/")) - 2] == "processed":
            visit = visitDirectory.split("/")[len(visitDirectory.split("/")) - 3]
        else:
            visit = visitDirectory.split("/")[len(visitDirectory.split("/")) - 1]
        beamline = "unknown"
    return visit, beamline


def crystal_growth_methods():

    methods = [
        "VAPOR DIFFUSION, SITTING DROP",
        "VAPOR DIFFUSION, HANGING DROP",
        "BATCH MODE",
        "LIPIDIC CUBIC PHASE",
        "MICROBATCH",
        "MICROFLUIDIC",
    ]

    return methods


def wwBeamlines():

    beamlines = [
        "DIAMOND BEAMLINE I02",
        "DIAMOND BEAMLINE I03",
        "DIAMOND BEAMLINE I04",
        "DIAMOND BEAMLINE I04-1",
        "DIAMOND BEAMLINE I23",
        "DIAMOND BEAMLINE I24",
    ]

    return beamlines


def radiationSource():

    source = ["SYNCHROTRON", "ROTATING ANODE", "SEALED TUBE"]

    return source


def detector():

    detectorPrinciple = ["PIXEL", "CCD", "IMAGE PLATE", "CMOS"]

    return detectorPrinciple


def detectorType():

    detector = [
        "DECTRIS PILATUS 2M",
        "DECTRIS PILATUS 2M-F",
        "DECTRIS PILATUS 6M",
        "DECTRIS PILATUS 6M-F",
        "DECTRIS PILATUS 12M",
        "DECTRIS PILATUS3 2M",
        "DECTRIS PILATUS3 6M",
        "DECTRIS EIGER X 9M",
        "DECTRIS EIGER X 16M",
        "ADSC QUANTUM 315",
        "ADSC QUANTUM 315r",
    ]

    return detector


def NCBI_taxonomy_ID():

    taxonomy_dict = {
        "9606": "Homo sapiens",
        "10090": "Mus musculus",
        "7108": "SPODOPTERA FRUGIPERDA",
        "5693 ": "Trypanosoma cruzi",
        "1508227": "BAT SARS-LIKE CORONAVIRUS",
        "2697049": "SARS-CoV-2",
        "562": "Escherichia coli",
        "837": "Porphyromonas gingivalis",
    }

    return taxonomy_dict


def data_integration_software():

    software = ["XDS", "HKL", "DENZO", "DTREK", "MOSFLM"]

    return software


def phasing_software():

    software = [
        "REFMAC",
        "PHENIX",
        "SOLVE",
        "PHASER",
        "CNS",
        "XPLOR",
        "MLPHARE",
        "SHELX",
        "SNB",
        "BnP",
        "BP3",
        "SHARP",
        "PHASES",
        "WARP",
    ]

    return software


def pdbx_keywords():

    keywords = [
        "",
        "ALLERGEN",
        "ANTIBIOTIC",
        "ANTIFREEZE PROTEIN",
        "ANTIFUNGAL PROTEIN",
        "ANTIMICROBIAL PROTEIN",
        "ANTITOXIN",
        "ANTITUMOR PROTEIN",
        "ANTIVIRAL PROTEIN",
        "APOPTOSIS",
        "ATTRACTANT",
        "BIOSYNTHETIC PROTEIN",
        "BLOOD CLOTTING",
        "CARBOHYDRATE",
        "CELL ADHESION",
        "CELL CYCLE",
        "CELL INVASION",
        "CHAPERONE",
        "CHOLINE - BINDING PROTEIN",
        "CIRCADIAN CLOCK PROTEIN",
        "CONTRACTILE PROTEIN",
        "CYTOKINE",
        "CYTOSOLIC PROTEIN",
        "DE NOVO PROTEIN",
        "DNA",
        "DNA BINDING PROTEIN",
        "DNA - RNA HYBRID",
        "ELECTRON TRANSPORT",
        "ENDOCYTOSIS",
        "EXOCYTOSIS",
        "FLAVOPROTEIN",
        "FLUORESCENT PROTEIN",
        "GENE REGULATION",
        "HORMONE",
        "HYDROLASE",
        "IMMUNE SYSTEM",
        "IMMUNOSUPPRESSANT",
        "ISOMERASE",
        "LIGASE",
        "LIPID BINDING PROTEIN",
        "LIPID TRANSPORT",
        "LUMINESCENT PROTEIN",
        "LYASE",
        "MEMBRANE PROTEIN",
        "METAL BINDING PROTEIN",
        "METAL TRANSPORT",
        "MOTOR PROTEIN",
        "NEUROPEPTIDE",
        "NUCLEAR PROTEIN",
        "ONCOPROTEIN",
        "OXIDOREDUCTASE",
        "OXYGEN BINDING",
        "OXYGEN STORAGE",
        "OXYGEN TRANSPORT",
        "PEPTIDE BINDING PROTEIN",
        "PHOTOSYNTHESIS",
        "PLANT PROTEIN",
        "PROTEIN BINDING",
        "PROTEIN FIBRIL",
        "PROTEIN TRANSPORT",
        "PROTON TRANSPORT",
        "RECOMBINATION",
        "REPLICATION",
        "RIBOSOMAL PROTEIN",
        "RIBOSOME",
        "RNA",
        "RNA BINDING PROTEIN",
        "SIGNALING PROTEIN",
        "SPLICING",
        "STRUCTURAL GENOMICS",
        "STRUCTURAL PROTEIN",
        "SUGAR BINDING PROTEIN",
        "SURFACTANT PROTEIN",
        "TOXIN",
        "TRANSCRIPTION",
        "TRANSFERASE",
        "TRANSLATION",
        "TRANSLOCASE",
        "TRANSPORT PROTEIN",
        "UNKNOWN FUNCTION",
        "VIRAL PROTEIN",
        "VIRUS",
        "VIRUS LIKE PARTICLE",
    ]

    return keywords


def pdbx_country():

    countries = [
        "United Kingdom",
        "United States",
        "Japan",
        "Albania",
        "Andorra",
        "Argentina",
        "Armenia",
        "Australia",
        "Austria",
        "Azerbaijan",
        "Bahamas",
        "Bangladesh",
        "Barbados",
        "Belarus",
        "Belgium",
        "Brazil",
        "Bulgaria",
        "Canada",
        "Chile",
        "China",
        "Croatia",
        "Cuba",
        "Cyprus",
        "Czech Republic",
        "Denmark",
        "Estonia",
        "Finland",
        "France",
        "Germany",
        "Greece",
        "Hungary",
        "Iceland",
        "India",
        "Indonesia",
        "Iran, Islamic Republic Of",
        "Iraq",
        "Ireland",
        "Israel",
        "Italy",
        "Jamaica",
        "Jordan",
        "Kazakhstan",
        "Kenya",
        "Kiribati",
        "Korea, Republic Of",
        "Latvia",
        "Lithuania",
        "Luxembourg",
        "Malta",
        "Mexico",
        "Netherlands",
        "New Zealand",
        "Norway",
        "Pakistan",
        "Paraguay",
        "Peru",
        "Philippines",
        "Poland",
        "Portugal",
        "Romania",
        "Russian Federation",
        "Serbia",
        "Singapore",
        "Slovakia",
        "Slovenia",
        "South Africa",
        "Spain",
        "Sweden",
        "Switzerland",
        "Taiwan",
        "Thailand",
        "Turkey",
        "Ukraine",
        "United Arab Emirates",
        "Uruguay",
    ]
    return countries


def read_html(file_name):
    source_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "html_fragments",
        file_name + ".html",
    )
    with open(source_path) as file:
        return file.read()


def html_header():
    return read_html("header")


def html_ngl(firstPDB, firstEvent, firstMap, firstDiffMap, ligID):
    ligChain = ligID.split("-")[1]
    ligResid = ligID.split("-")[2]

    return read_html("ngl") % (
        firstPDB,
        firstEvent,
        firstMap,
        firstDiffMap,
        ligChain,
        ligResid,
    )


def html_download(protein_name):
    return read_html("download") % (
        protein_name,
        protein_name,
        protein_name,
        protein_name,
        protein_name,
    )


def html_guide():
    return read_html("guide")


def html_table_header():
    return read_html("table_header")


def html_table_row(
    xtalID,
    pdbID,
    ligID,
    compoundImage,
    residuePlot,
    pdb,
    event,
    thumbNail,
    resoHigh,
    spg,
    unitCell,
    FWT,
    DELFWT,
    ligConfidence,
    modelStatus,
):
    ligChain = ligID.split("-")[1]
    ligResid = ligID.split("-")[2]

    return read_html("table_row") % (
        xtalID,
        pdbID,
        pdbID,
        ligID,
        ligConfidence,
        modelStatus,
        compoundImage,
        residuePlot,
        pdbID,
        pdb,
        event,
        FWT,
        DELFWT,
        ligChain,
        ligResid,
        thumbNail,
        resoHigh,
        spg,
        unitCell,
        pdb.replace(".pdb", ""),
        ligID,
    )


def html_footer():
    return read_html("footer")


def coot_prepare_input(x, y, z, ligID, sampleDir, eventMap):

    os.chdir(sampleDir)
    cmd = (
        "# !/usr/bin/env coot\n"
        "# python script for coot - generated by dimple\n"
        "import coot\n"
        'set_nomenclature_errors_on_read("ignore")\n'
        'molecule = read_pdb("refine.split.bound-state.pdb")\n'
        "set_rotation_centre(%s, %s, %s)\n" % (x, y, z) + "set_zoom(30.)\n"
        "set_view_quaternion(-0.180532, -0.678828, 0, 0.711759)\n"
        'coot.handle_read_ccp4_map(("%s"),0)\n' % eventMap
        + 'coot.raster3d("%s.r3d")\n' % ligID
        + "coot_real_exit(0)\n"
    )
    f = open(ligID + ".py", "w")
    f.write(cmd)
    f.close()


def coot_write_raster_file(ligID, sampleDir):
    os.chdir(sampleDir)
    os.system("coot --no-graphics --no-guano --script %s.py" % ligID)


def render_scene(xtal, ligID, sampleDir):
    os.chdir(sampleDir)
    os.system("render < %s.r3d -png %s_%s.png" % (ligID, xtal, ligID))


def make_thumbnail(xtal, ligID, sampleDir):
    os.chdir(sampleDir)
    os.system(
        "convert -thumbnail 150x150 %s_%s.png %s_%s_thumb.png"
        % (xtal, ligID, xtal, ligID)
    )


def backup_soakDB(database, xce_logfile):
    Logfile = XChemLog.updateLog(xce_logfile)
    Logfile.insert(
        "backing up soakDB: " + database + str(datetime.now()).replace(" ", "_")
    )
    os.system(
        "/bin/cp %s %s"
        % (
            database,
            database + "." + str(datetime.now()).replace(" ", "_").replace(":", "-"),
        )
    )
