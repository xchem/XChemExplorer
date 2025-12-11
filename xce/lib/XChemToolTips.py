import getpass
import os


def dataset_task_tip():
    tip = "describes what you can do"
    return tip


def dataset_task_run_button_tip():
    tip = "Run dataset"
    return tip


def dataset_task_status_button_tip():
    tip = "Status dataset"
    return tip


def map_cif_file_task_tip():
    tip = "describes what you can do"
    return tip


def map_cif_file_task_run_button_tip():
    tip = "Run map_cif_file"
    return tip


def map_cif_file_task_status_button_tip():
    tip = "Status map_cif_file"
    return tip


def panddas_file_task_tip():
    tip = "describes what you can do"
    return tip


def panddas_file_task_run_button_tip():
    tip = "Run panddas_file"
    return tip


def panddas_file_task_status_button_tip():
    tip = "Status panddas_file"
    return tip


def refine_file_task_tip():
    tip = "describes what you can do"
    return tip


def refine_file_task_run_button_tip():
    tip = "Run refine_file"
    return tip


def refine_file_task_status_button_tip():
    tip = "Status refine_file"
    return tip


def update_from_datasource_button_tip():
    tip = "Status validation_file"
    return tip


def run_pandda_inspect_at_home(pandda_directory):
    instruction = (
        "\n\n"
        "Be sure to have pandda installed at home, and go to a clean subdirectory.\n"
        "From that directory, do the steps below.\n"
        "This moves the relevant files to your site so you can do the model building"
        " locally, and then moves the files back to Diamond.\n"
        "1.	run:  rsync -av %s@nx.diamond.ac.uk:%s .\n"
        % (getpass.getuser(), pandda_directory)
        + '2.	run: "pandda.inspect", and build all relevant models, etc.\n'
        "3.	run:  rsync -av * %s@nx.diamond.ac.uk:%s\n"
        % (getpass.getuser(), pandda_directory)
        + "Now proceed within XChemExplorer as before.\n"
    )

    print(instruction)


def deposition_interface_note():
    note = (
        "\n\n"
        "Note: you can use this mask to save identical information for ALL structures"
        " to be deposited.\n"
        "However, this may not be suitable in cases where the information is different"
        " for certain samples.\n"
        "In such cases, please use for example SQLiteBrowser to edit the relevant"
        " fields in the depositTable."
    )

    return note


def pandda_pre_run(reference_directory):
    msg = (
        "The aim of the pre-run is NOT to identify bound ligands,\n"
        "but to create mean ground state maps.\n"
        "Hence, the pre-run will only comprise 100 datasets.\n"
        "After the pre-run is finished use the resulting ground state mean maps\n"
        "to build a suitable reference model for the subsequent PanDDA production"
        " run.\n"
        "The appendix will determine the name of the folder where the results from\n"
        "the pre-run will be stored. It is used by the COOT plugin to distinguish\n"
        "between different pre-runs.\n"
        "The result from the pre-run will be stored in:\n"
        "%s/pannda_<appenddix\n"
        "The bullet points below highlight the next steps after the pre-run is"
        " finished:\n"
        '- PanDDA tab: run "Build ground state model" \n'
        "- MAPS tab: select ALL datasets \n"
        '- MAPS tab: press "Refresh reference file list"\n'
        "- MAPS tab: select ground state model and set as new reference\n"
        '- MAPS tab: press "Run DIMPLE on selected MTZ files"\n'
        '- PanDDA tab: run "pandda.analyse\n'
        '- run "pandda.analyse\n'
    )

    return msg


def deposition_introduction():
    msg = (
        "<b>DISCLAIMER:</b> The group deposition process in XCE is now partially deprecated. "
        "You may still be able to assemble the required files using the steps below with some "
        "additional patching, but please note that we are actively working on replacing this "
        "part of the pipeline.<br><br>"
        '<h3 style="margin-top: 10px; margin-bottom: 5px;">Uploading to Fragalysis</h3>'
        "Instructions for uploading and sharing your data via Fragalysis can be found "
        '<a href="https://xchem-align.readthedocs.io/en/latest/" style="color: blue;">here</a>.'
    )
    return msg


def deposition_introduction_link():
    # Deprecated - no longer used
    lnk = ""
    return lnk


def deposition_bound_state_prerequisites():
    msg = (
        "1. Event map to MTZ conversion.\n"
        "     All pandda event maps need to be converted to MTZ format."
        " This is currently not done automatically.\n"
        '     If you have not done so, select and run "Event Map -> SF"'
        " from the Hit Identification action box.\n"
        "     This may take a while, but you only need to run this once.\n"
        "2. Select datasets to be deposited.\n"
        '     Set crystals to "5-ready for deposition" in XCE refinement tab\n'
        "3. Enter additional data required for PDB deposition.\n"
        '     In the Deposition menu, select "Edit information".'
        ' Fill out all the required items and press "Save to Database".\n'
        "     Note: this needs to be done after the datasets have been selected"
        " for deposition."
    )
    return msg


def deposition_bound_state_preparation_step_one_text():
    msg = (
        "1. Press the button below to generate structure and structure factor mmcif"
        " files of all selected datasets.\n"
        "     Note: all previously generated mmcif files will be overwritten."
    )
    return msg


def deposition_bound_state_preparation_step_two_text():
    msg = (
        "2. Press the button below to copy all the mmcif files into the"
        ' "group deposition directory" (see settings tab).\n'
        "All mmcif files will be bundled into a single, bzipped tar archive"
        " which can be uploaded into via the PDB group deposition website."
    )
    return msg


def pdb_group_deposition_instruction_one():
    msg = (
        "3. Go to the group deposition website, create a session and upload the"
        " ligand-bound.tar.bz2 file from the group deposition directory."
    )
    return msg


def pdb_group_deposition_link():
    lnk = (
        "<a "
        'href="https://deposit-group.rcsb.rutgers.edu/groupdeposit/"'
        ">'https://deposit-group.rcsb.rutgers.edu/groupdeposit/'</a>"
    )
    return lnk


def pdb_group_deposition_instruction_two():
    msg = "     user: grouptester\n" "     password: !2016rcsbpdb "
    return msg


def deposition_ground_state_prerequisites():
    msg = (
        "1. Convert all apo MTZ files to mmcif format.\n"
        "     Swtich to the pandda tab, select the PanDDA directory which contains the"
        " analysis you want to deposit.\n"
        '     Then select and run "apo -> mmcif" from the Hit Identification'
        " action box.\n"
        "     Note: you only need to do this once."
    )
    return msg


def deposition_ground_state_preparation_step_one_text():
    msg = (
        "1. Select ground-state PDB file.\n"
        "     Note: the file is usually in the reference directory."
    )
    return msg


def deposition_ground_state_preparation_step_three_text():
    msg = (
        "1. Please check the settings tab that you have selected the correct"
        " pandda directory.\n"
        "   (Note: we will take all the apo mmcif files from this directory)\n"
        "   Current PanDDA directory:"
    )
    return msg


def deposition_ground_state_preparation_step_four_text():
    msg = "2. Add the ground-state entry to the database."
    return msg


def deposition_ground_state_preparation_step_five_text():
    msg = (
        "3. Enter meta-data for ground-state model:\n"
        '     - Open "Deposition -> Edit information"\n'
        "     - Fill out form or load .deposit file\n"
        '     - Press "Save to Database"\n'
        '     - Press "OK"'
    )
    return msg


def deposition_ground_state_preparation_step_six_text():
    msg = (
        "4. Prepare the ground-state mmcif file.\n"
        "     Note: the mmcif files are saved into the selected pandda directory"
    )
    return msg


def deposition_ground_state_preparation_step_seven_text():
    msg = (
        "5. Press the button below to copy the structire and structure factor mmcif"
        ' files into the "group deposition directory" (see settings tab).\n'
        "Both mmcif files will be bundled into a single, bzipped tar archive which can"
        " be uploaded into via the PDB group deposition website."
    )
    return msg


def deposition_ground_state_preparation_step_eight_text():
    msg = (
        "6. Go to the group deposition website, create a session and upload the"
        " ligand-bound.tar.bz2 file from the group deposition directory."
    )
    return msg


def after_deposition_step_one_text():
    msg = (
        "After you have successfully submitted the ligand-bound structures via the PDB"
        " group deposition interface, you will immediately get\n"
        "an email with the PDB codes. There will be a single line for each PDB"
        ' submission. Highlight and copy the text! Then go to the "Deposition" menu\n'
        'and select "Update DB with PDB codes". A pop-up window will appear, paste the'
        ' text into the window and press "Update Database".'
    )
    return msg


def deposition_bounnd_state_preparation_ignore_event_map():
    msg = (
        "do NOT include PanDDA event maps"
        " (ONLY USE IN CASE IF DATA WERE NOT ANALYSED WITH PANDDA!)"
    )
    return msg


def second_cif_file_info(cif_file):
    second_cif_file = str(cif_file)
    if os.path.isfile(second_cif_file):
        msg = (
            "You have selected the following restraints file\n"
            "for merging into the the ligand CIF files:\n" + second_cif_file + "\n"
            "Press NO in case you do not want to continue\n"
            "In case you want to use another file, open the preferences menu"
            " and set it there\n"
            "Please check the XCE logfile in case of unexpected behaviour"
        )
    elif second_cif_file.replace(" ", "") == "" or second_cif_file.lower() == "none":
        msg = "No restraints file was selected!"
    else:
        msg = (
            "The selected restraints file does not exist:\n" + second_cif_file + "\n"
            "Open the preferences menu and set it there\n"
        )
    return msg


def second_cif_file_not_exists():
    msg = (
        "The CIF file for the non-standard ligand does not exist! "
        "It was either not selected or the file was deleted/ moved in the meantime. "
        "Please check Menu -> Preferences."
    )
    return msg


def notification_about_changes_to_apo_deposition():
    msg = (
        "The ground-state deposition procedure has changes in XCE v1.5.0!\n\n"
        "Previously, users had to select a specific reference PDB file,"
        " which had to have a ground-state mean-map and logfile associated with it.\n"
        "However, there were two main problems: (i) the map to mtz conversion"
        " did sometimes lead to strange maps or maps in space group P1 only\n"
        "and (2) occasionally the resulting MTZ files had missing columns."
        " This lead to errors and delays during deposition, while inclusion of\n"
        "ground-state mean maps has little benefit since the mean map"
        " can be calculated with the deposited apo datasets. Hence, there is really\n"
        "no need to include it in the deposition.\n\n"
        "All you need to do is to select the relevant PanDDA directory"
        " and XCE will arbitrarily select a high resolution structure with low Rfree\n"
        "as the model for the deposition bundle and then put all structure factor"
        " MMCIF files into a single file. Note that it does not matter\n"
        "which PDB file gets chosen at this point since the positional modifications"
        " during refinement with REFMAC will be minimal.\n"
        "If you still want to make the ground-state mean map available,"
        " we would recommend to deposit ground-state mean map and \n"
        "everything else you think is relevant on ZENODO and ask the PDB annotator to"
        " include the respective DOI in the deposition\n"
    )
    return msg


def pandda_export_ligand_bound_models_only_disclaimer():
    msg = (
        "You have chosen to export and refine the ligand-bound model only. "
        "This is different to the original procedure described by"
        " Pearce et al. (doi 10.1107/S2059798317003412), "
        "which is based on refinement of an ensemble model consisting of ligand-bound"
        " state and confounding ground-state."
        "Please note that working with the ligand-bound model only usually works well"
        " for reasonably well-defined "
        "ligands and conformational changes, but it may not be ideal for weakly bound"
        " fragments and corresponding "
        "conformational changes."
    )
    return msg
