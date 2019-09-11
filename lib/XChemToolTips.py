# last edited: 03/08/2017, 17:00

import os
import getpass

def dataset_task_tip():
    tip =   (   'describes what you can do'   )
    return tip

def dataset_task_run_button_tip():
    tip =   (   'Run dataset'   )
    return tip

def dataset_task_status_button_tip():
    tip =   (   'Status dataset'   )
    return tip

def map_cif_file_task_tip():
    tip =   (   'describes what you can do'   )
    return tip

def map_cif_file_task_run_button_tip():
    tip =   (   'Run map_cif_file'   )
    return tip

def map_cif_file_task_status_button_tip():
    tip =   (   'Status map_cif_file'   )
    return tip

def panddas_file_task_tip():
    tip =   (   'describes what you can do'   )
    return tip

def panddas_file_task_run_button_tip():
    tip =   (   'Run panddas_file'   )
    return tip

def panddas_file_task_status_button_tip():
    tip =   (   'Status panddas_file'   )
    return tip

def refine_file_task_tip():
    tip =   (   'describes what you can do'   )
    return tip

def refine_file_task_run_button_tip():
    tip =   (   'Run refine_file'   )
    return tip

def refine_file_task_status_button_tip():
    tip =   (   'Status refine_file'   )
    return tip

def validation_file_task_tip():
    tip =   (   'describes what you can do'   )
    return tip

def validation_file_task_run_button_tip():
    tip =   (   'Run validation_file'   )
    return tip

def validation_file_task_status_button_tip():
    tip =   (   'Status validation_file'   )
    return tip

def update_from_datasource_button_tip():
    tip =   (   'Status validation_file'   )
    return tip




def load_samples_from_datasource():
    tip =   (   'loads all information\n'
                'which is currently in the\n'
                'data source'   )
    return tip

def save_samples_to_datasource():
    tip =   (   'saves all changes you made\n'
                'to any of the cells to the\n'
                'data source'   )
    return tip

def create_pdb_cif_png_files():
    tip =   (   'uses ACEDRG to create pdb,cif and png files for all compounds\n'
                'that have a smiles string assigned. The compounds ID\n'
                'serves as name if available.\n'
                'The files are written to initial_model/<sample ID>/compound\n'
                'data source'   )
    return tip

def check_status_create_pdb_cif_png_files():
    tip =   (   'the jobs are sent to a computer cluster if available\n'
                'gives you an overview about how it progresses')
    return tip

def run_pandda_inspect_at_home(pandda_directory):
#    instruction =   (   '\n\n'
#                        'Please follow the instructions below to run pandda.inspect at your home institution:\n\n'
#                        ' 1. Make sure that the pandda package is installed at your home institution.\n'
#                        '    For more details see: http://pandda.bitbucket.org\n\n'
#                        ' 2. At home: go to the folder where you want the files from pandda.analyse to be;\n'
#                        '    e.g. cd /home/tkrojer/fragment_screening\n\n'
#                        ' 3. At home: run the following command:\n'
#                        '    rsync -av %s@nx.diamond.ac.uk:%s .\n\n' %(getpass.getuser(),pandda_directory)+
#                        ' 4. At home: go into the pandda directory\n'
#                        '    cd %s\n\n' %pandda_directory[pandda_directory.rfind('/')+1:]+
#                        ' 5. At home: run pandda.inspect\n'
#                        '    pandda.inspect\n\n'
#                        ' 6. Once you finished inspecting your models, close pandda.inspect and copy the data back to DLS\n'
#                        '    with the following command:\n'
#                        '    rsync -av * %s@nx.diamond.ac.uk:%s\n\n' %(getpass.getuser(),pandda_directory)+
#                        ' 7. At DLS: continue using XChemExplorer; go to the PANDDA tab and run:\n'
#                        '    "Export PANDDA models"\n'
#                        '    This will trigger the following steps:\n'
#                        '    - transform the pandda models back into the crystallographic unit cell\n'
#                        '    - copy the transformed models into the respective folder of the project directory\n'
#                        '    - launch an initial round of refinement with REFMAC\n'
#                        '      Note: All refinement jobs run on the cluster at DLS.\n'
#                        '            You can check the status of the jobs by typing:\n'
#                        '            module load global/cluster\n'
#                        '            qstat\n\n'
#                        ' 8. At DLS: go to the REFINEMENT tab and run "Open COOT" to review and further refine your models if necessary\n' )

    instruction =   (   '\n\n'
                        'Be sure to have pandda installed at home, and go to a clean subdirectory.\n'
                        'From that directory, do the steps below.\n'
                        'This moves the relevant files to your site so you can do the model building locally, and then moves the files back to Diamond.\n'
                        '1.	run:  rsync -av %s@nx.diamond.ac.uk:%s .\n' %(getpass.getuser(),pandda_directory)+
                        '2.	run: "pandda.inspect", and build all relevant models, etc.\n'
                        '3.	run:  rsync -av * %s@nx.diamond.ac.uk:%s\n' %(getpass.getuser(),pandda_directory)+
                        'Now proceed within XChemExplorer as before.\n' )

    print instruction

def deposition_interface_note():

    note = (    '\n\n'
                'Note: you can use this mask to save identical information for ALL structures to be deposited.\n'
                'However, this may not be suitable in cases where the information is different for certain samples.\n'
                'In such cases, please use for example SQLiteBrowser to edit the relevant fields in the depositTable.'   )

    return note

def html_summary_introduction():

    text = (    'This section describes how to generate an HTML summary page of the current project. \n'
                'All you need to do is select the HTML output directory where the html and additional auxiliary files \n'
                'and will be written to. '
                'Then press the "Export to HTML" button and \nall bound structures that are labeled as "4 - compchem ready"'
                'or higher will be included in the final document. \n'
                'Please checkout https://www.thesgc.org/fragment-screening to see what the final pages look like. ')

    return text

def html_export_directory_background():

    text = (    'The HTML export directory will contain all files necessary to generate the summary file.\n'
                'It is also makes it easy to share the results with collaborators.')

    return text

def html_export_step():

    text = (    'This step will create an index.html file and it will copy pdb, mtz and cif files\n'
                'as well as PanDDA event maps into the html export directory.'  )

    return text


def icb_file_background():

    text = (    'This step will use ICM-pro to create interactive structure thumbnails in the html file (http://www.molsoft.com/activeicmjs.html)')

    return text



def prepare_ICB_files():

    instruction = (     '- Open file browser and navigate to HTML export directory.\n'
                        '- Open ICM (at DLS: select SERVER and enter diamvicmpro.diamond.ac.uk in case ICM asks).\n'
                        '- Drag-and-drop dsEvent_sqlite.icm file into ICM main window.\n'
                        '- In ICM workspace panel: right-click on dsEvent_sqlite.icm and select "RUN".')
    return instruction

def zenodo_upload_start(html_export_directory):

    instruction = (     'This will copy all files required for ZENODO upload to\n'
                        '%s'   %os.path.join(html_export_directory,'zenodo')    )
    return instruction

def zenodo_upload_part_one(html_export_directory):

    instruction = (     '- Register with ZENODO (www.zenodo.org)\n'
                        '- Select Upload -> New Upload\n'
                        '- upload all files from %s/zenodo, BUT NOT any of the .html files!!!' %html_export_directory  )
    return instruction

def zenodo_upload_part_two():

    instruction = (     'You might have noticed that once you started uploading files, the address bar in your browser will change to something like this:')

    return instruction

def zenodo_upload_part_three():

    instruction = (     'Please enter the upload ID, i.e. the number at the end of the line, into the field below')

    return instruction

def zenodo_upload_part_four(html_export_directory):

    instruction = (     'Upload ALL html files in {0!s} to ZENODO and publish the page.'.format(os.path.join(html_export_directory,'zenodo')))

    return instruction

def deposition_pandda_site_not_ready(xtal):
    msg = (     'Please make sure that all placed ligands in '+xtal+' are ready for deposition.\n'
                'Go to the Refinement tab -> Open COOT and either delete the ligand or set the site '
                ' to "5 - Deposition ready"'    )
    return msg

def pandda_pre_run(reference_directory):
    msg = (    'The aim of the pre-run is NOT to identify bound ligands,\n'
                'but to create mean ground state maps.\n'
                'Hence, the pre-run will only comprise 100 datasets.\n'
                'After the pre-run is finished use the resulting ground state mean maps\n'
                'to build a suitable reference model for the subsequent PanDDA production run.\n'
                'The appendix will determine the name of the folder where the results from\n'
                'the pre-run will be stored. It is used by the COOT plugin to distinguish\n'
                'between different pre-runs.\n'
                'The result from the pre-run will be stored in:\n'
                '%s/pannda_<appenddix\n'
                'The bullet points below highlight the next steps after the pre-run is finished:\n'
                        '- PanDDA tab: run "Build ground state model" \n'
                        '- MAPS tab: select ALL datasets \n'
                        '- MAPS tab: press "Refresh reference file list"\n'
                        '- MAPS tab: select ground state model and set as new reference\n'
                        '- MAPS tab: press "Run DIMPLE on selected MTZ files"\n'
                        '- PanDDA tab: run "pandda.analyse\n'
                        '- run "pandda.analyse\n'    )

    return msg

def deposition_introduction():
    msg = ( 'Some background about batch deposition of PanDDA models can be found here: '
    )
    return msg

def deposition_introduction_link():
    lnk = (
        "<a href=\"https://openlabnotebooks.org/update-on-batch-deposition-of-xchem-structures\">'https://openlabnotebooks.org/update-on-batch-deposition-of-xchem-structures'</a>"
    )
    return lnk

def deposition_bound_state_prerequisites():
    msg = (
        '1. Event map to MTZ conversion.\n'
        '     All pandda event maps need to be converted to MTZ format. This is currently not done automatically.\n'
        '     If you have not done so, select and run "Event Map -> SF" from the Hit Identification action box.\n'
        '     This may take a while, but you only need to run this once.\n'
        '2. Select datasets to be deposited.\n'
        '     Set crystals to "5-ready for deposition" in XCE refinement tab\n'
        '3. Enter additional data required for PDB deposition.\n'
        '     In the Deposition menu, select "Edit information". Fill out all the required items and press "Save to Database".\n'
        '     Note: this needs to be done after the datasets have been selected for deposition.'
    )
    return msg

def deposition_bound_state_preparation_step_one_text():
    msg = (
        '1. Press the button below to generate structure and structure factor mmcif files of all selected datasets.\n'
        '     Note: all previously generated mmcif files will be overwritten.'
    )
    return msg

def deposition_bound_state_preparation_step_two_text():
    msg = (
        '2. Press the button below to copy all the mmcif files into the "group deposition directory" (see settings tab).\n'
        'All mmcif files will be bundled into a single, bzipped tar archive which can be uploaded into via the PDB group deposition website.'
    )
    return msg

def pdb_group_deposition_instruction_one():
    msg = (
        '3. Go to the group deposition website, create a session and upload the ligand-bound.tar.bz2 file from the group deposition directory.'
    )
    return msg

def pdb_group_deposition_link():
    lnk = (
        "<a href=\"https://deposit-group-1.rcsb.rutgers.edu/groupdeposit\">'https://deposit-group-1.rcsb.rutgers.edu/groupdeposit'</a>"
    )
    return lnk

def pdb_group_deposition_instruction_two():
    msg = (
        '     user: grouptester\n'
        '     password: !2016rcsbpdb '
    )
    return msg

def deposition_ground_state_prerequisites():
    msg = (
        '1. Convert all apo MTZ files to mmcif format.\n'
        '     Swtich to the pandda tab, select the PanDDA directory which contains the analysis you want to deposit.\n'
        '     Then select and run "apo -> mmcif" from the Hit Identification action box.\n'
        '     Note: you only need to do this once.'
    )
    return msg


def deposition_ground_state_preparation_step_one_text():
    msg = (
        '1. Select ground-state PDB file.\n'
        '     Note: the file is usually in the reference directory.'
    )
    return msg

def deposition_ground_state_log_info():
    msg = (
            '     IMPORTANT:\n'
            '     Make sure that a file or symlink to an AIMLESS logile with the same name root as your ground-state PDB file\n'
            '     is present in the same folder. This is currently not generated automatically! '
    )
    return msg

def deposition_ground_state_preparation_step_two_text():
    msg = (
        '2. Select ground-state MTZ file.\n'
        '     Notes: - this is the MTZ file of the ground-state mean map\n'
        '            - the file should have been created after the pandda pre-run\n'
        '            - hence, you will usually find it in the reference directory'
    )
    return msg

def deposition_ground_state_preparation_step_three_text():
    msg = (
        '3. Please check the settings tab that you have selected the correct pandda directory.\n'
        '   (Note: we will take all the apo mmcif files from this directory)\n'
        '   Current PanDDA directory:'
    )
    return msg

def deposition_ground_state_preparation_step_four_text():
    msg = (
        '4. Add the ground-state entry to the database.'
    )
    return msg

def deposition_ground_state_preparation_step_five_text():
    msg = (
            '5. Enter meta-data for ground-state model:\n'
            '     - Open "Deposition -> Edit information"\n'
            '     - Fill out form or load .deposit file\n'
            '     - Press "Save to Database"\n'
            '     - Press "OK"'
    )
    return msg

def deposition_ground_state_preparation_step_six_text():
    msg = (
            '6. Prepare the ground-state mmcif file.\n'
            '     Note: the mmcif files are saved into the selected pandda directory'
    )
    return msg

def deposition_ground_state_preparation_step_seven_text():
    msg = (
        '7. Press the button below to copy the structire and structure factor mmcif files into the "group deposition directory" (see settings tab).\n'
        'Both mmcif files will be bundled into a single, bzipped tar archive which can be uploaded into via the PDB group deposition website.'
    )
    return msg

def deposition_ground_state_preparation_step_eight_text():
    msg = (
        '8. Go to the group deposition website, create a session and upload the ligand-bound.tar.bz2 file from the group deposition directory.'
    )
    return msg

def after_deposition_step_one_text():
    msg = (
        'After you have successfully submitted the ligand-bound structures via the PDB group deposition interface, you will immediately get\n'
        'an email with the PDB codes. There will be a single line for each PDB submission. Highlight and copy the text! Then go to the "Deposition" menu\n'
        'and select "Update DB with PDB codes". A pop-up window will appear, paste the text into the window and press "Update Database".'
    )
    return msg

def deposition_bounnd_state_preparation_ignore_event_map():
    msg = (
        'do NOT include PanDDA event maps (ONLY USE IN CASE IF DATA WERE NOT ANALYSED WITH PANDDA!)'
    )
    return msg

def second_cif_file_info(second_cif_file):
    if os.path.isfile(str(second_cif_file)):
        msg = (
            'You have selected the following restraints file\n'
            'for merging into the the ligand CIF files:\n'
            +second_cif_file+
            '\n'
            'Press NO in case you do not want to continue\n'
            'In case you want to use another file, open the preferences menu and set it there\n'
            'Please check the XCE logfile in case of unexpected behaviour'
        )
    elif second_cif_file.replace(' ','') == '' or second_cif_file.lower() == 'none':
        msg = (
            'No restraints file was selected!'
        )
    else:
        msg = (
            'The selected restraints file does not exist:\n'
            +second_cif_file+
            '\n'
            'Open the preferences menu and set it there\n'
        )
    return msg