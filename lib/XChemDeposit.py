# last edited: 16/05/2017, 15:00

from XChemUtils import parse
from XChemUtils import smilestools
from XChemUtils import mtztools
from XChemUtils import pdbtools
import XChemMain
import XChemToolTips
import XChemDB
import XChemLog
import sys
import fileinput
import os
import glob

from PyQt4 import QtGui, QtCore, QtWebKit

from lxml import etree

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')


def create_SF_mmcif(outDir, mtzList):
    print('hallo')


def get_protein_sequence(database, xtalID):
    print('hallo')


def check_depositDict(depositDict):
    # check depositDict
    for entry in depositDict:
        if 'middle_name' in depositDict[entry]:
            continue
        elif 'State_or_Province' in depositDict[entry]:
            continue
        elif depositDict[entry] == '':
            print('ERROR')


def update_title(depositDict):

    print('hallo')


def create_data_template_text():

    data_template_text = data_template(depositDict, sequence)


def create_Model_mmcif(outDir, pdbList):
    print('hallo')

# def update_file_locations_of_apo_structuresin_DB(database,projectDir,xce_logfile):
#    Logfile=XChemLog.updateLog(xce_logfile)
#    Logfile.insert('updating file information for apo structures')
#    db=XChemDB.data_source(database)
#    apo=db.execute_statement("select CrystalName from depositTable where StructureType is 'apo';")
#    for item in apo:
#        xtal=str(item[0])
#        db_dict={}
#        db_dict['label']=xtal+'-apo'
#        db_dict['description']='apo structure for pandda.analyse'
#        if os.path.isfile(os.path.join(projectDir,xtal,'dimple.pdb')):
#            db_dict['PDB_file']=os.path.realpath(os.path.join(projectDir,xtal,'dimple.pdb'))
#            if os.path.isfile(os.path.join(projectDir,xtal,'dimple.mtz')):
#                db_dict['MTZ_file']=os.path.realpath(os.path.join(projectDir,xtal,'dimple.mtz'))
#                Logfile.insert('updating depositTable for apo structure '+xtal)
#                db.update_depositTable(xtal,'apo',db_dict)


class templates:

    def data_template_cif(self, depositDict):

        taxonomy_dict = XChemMain.NCBI_taxonomy_ID()
        for key in taxonomy_dict:
            if taxonomy_dict[key] == depositDict['Source_organism_scientific_name']:
                pdbx_gene_src_ncbi_taxonomy_id = key
            if taxonomy_dict[key] == depositDict['Expression_system_scientific_name']:
                pdbx_host_org_ncbi_taxonomy_id = key

#        structure_author_name=''
#        for name in depositDict['structure_author_name'].split(';'):
#            structure_author_name+='<structure_author_name=  %s>\n' %name

        audit_author_name = ''
        # one name must be within quotation, last name and first initial must be separated by comma and space
        for name in depositDict['structure_author_name'].split(';'):
            if name.replace(' ', '') == '':
                continue
            if name[name.find(',')+1:name.find(',')+2] != ' ':
                name = name.replace(',', ', ')
            audit_author_name += "'{0!s}'\n".format(name)

        primary_citation_author_name = ''
        # one name must be within quotation, last name and first initial must be separated by comma and space
        for name in depositDict['primary_citation_author_name'].split(';'):
            if name.replace(' ', '') == '':
                continue
            if name[name.find(',')+1:name.find(',')+2] != ' ':
                name = name.replace(',', ', ')
            primary_citation_author_name += "primary '{0!s}'\n".format(name)

        molecule_one_letter_sequence = ';'
        counter = 1
        for aa in depositDict['molecule_one_letter_sequence']:
            if counter < 70:
                molecule_one_letter_sequence += aa
            if counter == 70:
                molecule_one_letter_sequence += '\n'+aa
                counter = 0
            counter += 1

        if depositDict['molecule_name_two'].replace(' ', '') == '' or depositDict['molecule_name_two'].replace(' ', '').lower() == 'none':
            try:
                entity = (
                    'loop_\n'
                    '_entity.id\n'
                    '_entity.type\n'
                    '_entity.src_method\n'
                    '_entity.pdbx_description\n'
                    '_entity.pdbx_mutation\n'
                    '1 polymer     man "%s" %s\n' % (depositDict['Source_organism_gene'], depositDict['fragment_name_one_specific_mutation']) +
                    '#\n'
                    'loop_\n'
                    '_entity_poly.entity_id\n'
                    '_entity_poly.type\n'
                    '_entity_poly.pdbx_seq_one_letter_code\n'
                    '_entity_poly.pdbx_strand_id\n'
                    '_entity_poly.pdbx_seq_db_id\n'
                    '_entity_poly.pdbx_seq_db_name\n'
                    '1 "polypeptide(L)"\n'
                    + molecule_one_letter_sequence + '\n'
                    ';\n'
                    '%s %s UNP\n' % (depositDict['protein_chains'], depositDict['molecule_one_letter_sequence_uniprot_id']) +
                    '#\n'
                    'loop_\n'
                    '_entity_src_gen.entity_id\n'
                    '_entity_src_gen.gene_src_strain\n'
                    '_entity_src_gen.pdbx_gene_src_scientific_name\n'
                    '_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id\n'
                    '_entity_src_gen.pdbx_host_org_scientific_name\n'
                    '_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id\n'
                    '1 ? "%s" %s  "%s" %s\n' % (depositDict['Source_organism_scientific_name'], pdbx_gene_src_ncbi_taxonomy_id, depositDict['Expression_system_scientific_name'], pdbx_host_org_ncbi_taxonomy_id) +
                    '#\n'
                )
            except UnboundLocalError:
                print('Error: Something went wrong! Please check if you have saved the .deposit file to the database: Menu Deposition -> Edit Information')
        else:
            molecule_two_letter_sequence = ';'
            counter = 1
            for aa in depositDict['molecule_two_letter_sequence']:
                if counter < 70:
                    molecule_two_letter_sequence += aa
                if counter == 70:
                    molecule_two_letter_sequence += '\n'+aa
                    counter = 0
                counter += 1

            entity = (
                'loop_\n'
                '_entity.id\n'
                '_entity.type\n'
                '_entity.src_method\n'
                '_entity.pdbx_description\n'
                '_entity.pdbx_mutation\n'
                '1 polymer     man "%s" %s\n' % (depositDict['Source_organism_gene'], depositDict['fragment_name_one_specific_mutation']) +
                '2 polymer     man "%s" %s\n' % (depositDict['Source_organism_gene_two'], depositDict['fragment_name_two_specific_mutation']) +
                '#\n'
                'loop_\n'
                '_entity_poly.entity_id\n'
                '_entity_poly.type\n'
                '_entity_poly.pdbx_seq_one_letter_code\n'
                '_entity_poly.pdbx_strand_id\n'
                '_entity_poly.pdbx_seq_db_id\n'
                '_entity_poly.pdbx_seq_db_name\n'
                '1 "polypeptide(L)"\n'
                + molecule_one_letter_sequence + '\n'
                ';\n'
                '%s %s UNP\n' % (depositDict['molecule_chain_one'], depositDict['molecule_one_letter_sequence_uniprot_id']) +
                #                    ';\n'
                '2 "polypeptide(L)"\n'
                + molecule_two_letter_sequence + '\n'
                ';\n'
                '%s %s UNP\n' % (depositDict['molecule_chain_two'], depositDict['molecule_two_letter_sequence_uniprot_id']) +
                '#\n'
                'loop_\n'
                '_entity_src_gen.entity_id\n'
                '_entity_src_gen.gene_src_strain\n'
                '_entity_src_gen.pdbx_gene_src_scientific_name\n'
                '_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id\n'
                '_entity_src_gen.pdbx_host_org_scientific_name\n'
                '_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id\n'
                '1 ? "%s" %s  "%s" %s\n' % (depositDict['Source_organism_scientific_name'], pdbx_gene_src_ncbi_taxonomy_id, depositDict['Expression_system_scientific_name'], pdbx_host_org_ncbi_taxonomy_id) +
                '2 ? "%s" %s  "%s" %s\n' % (depositDict['Source_organism_scientific_name'], pdbx_gene_src_ncbi_taxonomy_id, depositDict['Expression_system_scientific_name'], pdbx_host_org_ncbi_taxonomy_id) +
                '#\n'
            )

        data_template_cif = (
            'data_UNNAMED\n'
            '#\n'
            '_pdbx_database_status.entry_id                       UNNAMED\n'
            "_pdbx_database_status.dep_release_code_coordinates   '%s'\n" % depositDict['Release_status_for_coordinates'] +
            "_pdbx_database_status.dep_release_code_sequence      '{0!s}'\n".format(depositDict['Release_status_for_sequence']) +
            '#\n'
            '_pdbx_deposit_group.group_id	   UNNAMED\n'
            '_pdbx_deposit_group.group_description  "%s"\n' % depositDict['group_description'] +
            '_pdbx_deposit_group.group_title        "{0!s}"\n'.format(depositDict['group_title']) +
            '_pdbx_deposit_group.group_type         "{0!s}"\n'.format(depositDict['group_type']) +
            '#\n'
            '_exptl_crystal_grow.crystal_id      1\n'
            "_exptl_crystal_grow.method          '%s'\n" % depositDict['crystallization_method'] +
            '_exptl_crystal_grow.pH              {0!s}\n'.format(depositDict['crystallization_pH']) +
            '_exptl_crystal_grow.temp            {0!s}\n'.format(depositDict['crystallization_temperature']) +
            '_exptl_crystal_grow.pdbx_details    "{0!s}"\n'.format(depositDict['crystallization_details']) +
            '#\n'
            '_diffrn.id                     1\n'
            '_diffrn.ambient_temp           %s\n' % depositDict['data_collection_temperature'] +
            '_diffrn.crystal_id             1\n'
            '#\n'
            '_diffrn_source.diffrn_id                       1\n'
            '_diffrn_source.source                          %s\n' % depositDict['radiation_source'] +
            '_diffrn_source.type                            "{0!s}"\n'.format(depositDict['radiation_source_type']) +
            '_diffrn_source.pdbx_wavelength_list            {0!s}\n'.format(depositDict['radiation_wavelengths']) +
            '#\n'
            '_diffrn_detector.detector               %s\n' % depositDict['radiation_detector'] +
            "_diffrn_detector.type                   '{0!s}'\n".format(depositDict['radiation_detector_type']) +
            '_diffrn_detector.pdbx_collection_date   {0!s}\n'.format(depositDict['data_collection_date']) +
            '_diffrn_detector.diffrn_id              1\n'
            '#\n'
            '_diffrn_radiation.diffrn_id                        1\n'
            '_diffrn_radiation.wavelength_id                    1\n'
            "_diffrn_radiation.pdbx_diffrn_protocol             'SINGLE WAVELENGTH'\n"
            '#\n'
            '_diffrn_radiation_wavelength.id           1\n'
            '_diffrn_radiation_wavelength.wavelength   %s\n' % depositDict['radiation_wavelengths'] +
            '#\n'
            '#\n'
            + entity +
            #            'loop_\n'
            #            '_entity.id\n'
            #            '_entity.type\n'
            #            '_entity.src_method\n'
            #            '_entity.pdbx_description\n'
            #            '_entity.pdbx_mutation\n'
            #            '1 polymer     man "%s" %s\n'                                                          %(depositDict['Source_organism_gene'],depositDict['fragment_name_one_specific_mutation'])+
            #            '#\n'
            #            'loop_\n'
            #            '_entity_poly.entity_id\n'
            #            '_entity_poly.type\n'
            #            '_entity_poly.pdbx_seq_one_letter_code\n'
            #            '_entity_poly.pdbx_strand_id\n'
            #            '_entity_poly.pdbx_seq_db_id\n'
            #            '_entity_poly.pdbx_seq_db_name\n'
            #            '1 "polypeptide(L)"\n'
            #            +molecule_one_letter_sequence+'\n'
            #            ';\n'
            #            '%s %s UNP\n'                                        %(depositDict['protein_chains'],depositDict['molecule_one_letter_sequence_uniprot_id'])+
            #            '#\n'
            #            'loop_\n'
            #            '_entity_src_gen.entity_id\n'
            #            '_entity_src_gen.gene_src_strain\n'
            #            '_entity_src_gen.pdbx_gene_src_scientific_name\n'
            #            '_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id\n'
            #            '_entity_src_gen.pdbx_host_org_scientific_name\n'
            #            '_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id\n'
            #            '1 ? "%s" %s  "%s" %s\n'                %(depositDict['Source_organism_scientific_name'],pdbx_gene_src_ncbi_taxonomy_id,depositDict['Expression_system_scientific_name'],pdbx_host_org_ncbi_taxonomy_id)+
            #            '#\n'
            'loop_\n'
            '_pdbx_contact_author.id                  \n'
            "_pdbx_contact_author.address_1           \n"
            '_pdbx_contact_author.address_2           \n'
            '_pdbx_contact_author.city                \n'
            "_pdbx_contact_author.state_province      \n"
            '_pdbx_contact_author.postal_code         \n'
            '_pdbx_contact_author.email               \n'
            '_pdbx_contact_author.name_first          \n'
            '_pdbx_contact_author.name_last           \n'
            '_pdbx_contact_author.country             \n'
            '_pdbx_contact_author.phone               \n'
            '_pdbx_contact_author.role                \n'
            '_pdbx_contact_author.organization_type   \n'
            '_pdbx_contact_author.identifier_ORCID    \n'
            "1 '%s' '%s' '%s' '?' '%s' %s %s '%s' '%s' '%s' '%s' %s %s\n" % (depositDict['contact_author_PI_address'], depositDict['contact_author_PI_organization_name'], depositDict['contact_author_PI_city'], depositDict['contact_author_PI_Zip_Code'], depositDict['contact_author_PI_email'], depositDict['contact_author_PI_first_name'], depositDict['contact_author_PI_last_name'], depositDict['contact_author_PI_Country'], depositDict['contact_author_PI_phone_number'], depositDict['contact_author_PI_role'], depositDict['contact_author_PI_organization_type'], depositDict['contact_author_PI_ORCID']) +
            "2 '{0!s}' '{1!s}' '{2!s}' '?' '{3!s}' '{4!s}' {5!s} {6!s} '{7!s}' '{8!s}' '{9!s}' '{10!s}' {11!s}\n".format(depositDict['contact_author_address'], depositDict['contact_author_organization_name'], depositDict['contact_author_city'], depositDict['contact_author_Zip_Code'].replace(' ', ''), depositDict['contact_author_email'], depositDict['contact_author_first_name'], depositDict['contact_author_last_name'], depositDict['contact_author_Country'], depositDict['contact_author_phone_number'], depositDict['contact_author_role'], depositDict['contact_author_organization_type'], depositDict['contact_author_ORCID']) +
            '#\n'
            'loop_\n'
            '_audit_author.name\n'
            + audit_author_name +
            #            "'%s, %s.'\n" %(depositDict['contact_author_last_name'],depositDict['contact_author_first_name'][0])+
            '#\n'
            '_citation.id                        primary\n'
            "_citation.title                     '%s'\n" % depositDict['group_title'] +
            "_citation.journal_abbrev            'To Be Published'\n"
            '#\n'
            'loop_\n'
            '_citation_author.citation_id\n'
            '_citation_author.name\n'
            + primary_citation_author_name +
            #        "primary 'Krojer, T.'\n"
            #        "primary 'Von Delft, F.'\n"
            '#\n'
            '_struct.entry_id                     UNNAMED\n'
            '_struct.title\n'
            ';%s\n' % depositDict['title'] +
            ';\n'
            '#\n'
            '_struct_keywords.entry_id        UNNAMED\n'
            '_struct_keywords.text            "%s"\n' % depositDict['structure_keywords'] +
            '#\n'
            '_pdbx_struct_assembly_depositor_info.id                   1\n'
            "_pdbx_struct_assembly_depositor_info.method_details       PISA\n"
            '_pdbx_struct_assembly_depositor_info.oligomeric_count     %s\n' % depositDict['biological_assembly_chain_number'] +
            '#\n'
        )

        return data_template_cif


class update_depositTable(QtCore.QThread):
    def __init__(self, deposit_dict, database, xce_logfile):
        QtCore.QThread.__init__(self)
        self.deposit_dict = deposit_dict
        self.database = database
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.db = XChemDB.data_source(database)

    def run(self):
        self.Logfile.insert(
            'all entries in the depositTable will be updated with the following values:')
        for key in self.deposit_dict:
            self.Logfile.insert(key+': '+self.deposit_dict[key])
        dbEntries = self.db.execute_statement(
            "select CrystalName,StructureType from depositTable;")
        for item in dbEntries:
            xtal = str(item[0])
            type = str(item[1])
            # need to do this because individual fields might need updating for some xtals
            db_dict = self.deposit_dict

            # try to get information about the diffraction experiment
            try:
                diffractionExperiment = self.db.execute_statement(
                    "select DataCollectionBeamline,DataCollectionDate from mainTable where CrystalName is '{0!s}'".format(xtal))
                beamline = str(diffractionExperiment[0][0])
                date = str(diffractionExperiment[0][1])
            except (UnboundLocalError, IndexError):
                self.Logfile.warning(
                    '%s: cannot find details about diffraction experiment in mainTable' % xtal)
                beamline = db_dict['radiation_source']
                date = db_dict['data_collection_date']
                self.Logfile.warning(
                    '%s: using values provided in depositTable for beamline and data collection date' % xtal)
            if beamline.lower() != 'none':
                db_dict = self.tweak_deposit_dict(xtal, db_dict)
            if date.lower() != 'none':
                db_dict['data_collection_date'] = date.split()[0]

            self.Logfile.insert('updating depositTable for '+xtal+' @ '+type)
            self.db.update_depositTable(xtal, type, db_dict)
        self.Logfile.insert('Note: use DBbrowser to edit individual entries')

    def tweak_deposit_dict(self, xtal, db_dict):
        dls_beamlines = ['i02', 'i03', 'i04', 'i04-1', 'i23', 'i24']
        dls_beamline_dict = {'i02':      ['DIAMOND BEAMLINE I02',    'DECTRIS PILATUS 6M'],
                             'i03':      ['DIAMOND BEAMLINE I03',    'DECTRIS PILATUS 6M'],
                             'i04':      ['DIAMOND BEAMLINE I04',    'DECTRIS PILATUS 6M'],
                             'i04-1':    ['DIAMOND BEAMLINE I04-1',  'DECTRIS PILATUS 6M'],
                             'i23':      ['DIAMOND BEAMLINE I23',    'DECTRIS PILATUS 12M'],
                             'i24':      ['DIAMOND BEAMLINE I24',    'DECTRIS PILATUS 6M'], }

        if db_dict['radiation_source_type'] in dls_beamlines:
            db_dict['radiation_source_type'] = dls_beamline_dict[db_dict['radiation_source_type']][0]
            db_dict['radiation_detector_type'] = dls_beamline_dict[db_dict['radiation_source_type']][1]
            db_dict['radiation_detector'] = 'PIXEL'
            db_dict['radiation_source'] = 'SYNCHROTRON'

        return db_dict


class prepare_mmcif_files_for_deposition(QtCore.QThread):

    def __init__(self, database, xce_logfile, overwrite_existing_mmcif, projectDir, ground_state, ignore_event_map):
        QtCore.QThread.__init__(self)
        self.database = database
        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.db = XChemDB.data_source(database)
        self.overwrite_existing_mmcif = overwrite_existing_mmcif
        self.projectDir = projectDir
        self.data_template_dict = {}

        self.errorList = []
        self.eventList = []
        self.db_dict = None
        self.data_template_dict = None
        self.zenodo_dict = None
        self.pdb = None
        self.mtz = None
        self.logDir = None

        self.ground_state = False
        self.ground_state_pdb = ''
#        self.ground_state_mean_mtz = ''
        self.panddaDir = ''
        self.ignore_event_map = ignore_event_map
        if ground_state:
            self.ground_state = True
            self.ground_state_pdb = ground_state[0]
            self.ground_state_mtz = ground_state[1]
            self.panddaDir = ground_state[2]
            self.logDir = self.projectDir
            self.projectDir = self.panddaDir
            self.pdb = pdbtools(self.ground_state_pdb)
            self.mtz = mtztools(self.ground_state_mtz)
#            self.mtz = mtztools(self.ground_state_mean_mtz)

    def run(self):

        self.Logfile.insert(
            '======= preparing mmcif files for wwPDB deposition =======')
        self.Logfile.insert('checking DB for structures to deposit...')
        if self.ground_state:
            toDeposit = self.db.execute_statement(
                "select CrystalName from depositTable where DimplePANDDApath = '%s';" % self.panddaDir)
        else:
            toDeposit = self.db.execute_statement(
                "select CrystalName from mainTable where RefinementOutcome like '5%';")
        self.Logfile.insert('found ' + str(len(toDeposit)) +
                            ' samples ready for deposition')

        progress_step = 1
        if len(toDeposit) != 0:
            progress_step = 100/float(len(toDeposit))
        else:
            progress_step = 1
        progress = 0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        for item in sorted(toDeposit):
            xtal = str(item[0])
            if self.ground_state:
                # projectDir == referenceDir in this case
                os.chdir(self.projectDir)
            else:
                os.chdir(os.path.join(self.projectDir, xtal))
            self.Logfile.insert(
                '%s: ----- preparing files for deposition -----' % xtal)

            if self.ground_state:
                if not self.data_template_dict_exists(xtal):
                    continue

#                if not self.zenodo_dict_exists(xtal):
#                    continue

                if not self.save_data_template_dict(xtal):
                    continue

                if not self.create_model_mmcif(xtal):
                    continue

                self.add_funding_information(xtal)

#                if not self.create_sf_mmcif(xtal):
#                    continue

                if not self.apo_mmcif_exists():
                    continue

                if not self.add_apo_sf_mmcif_to_ground_state_mmcif():
                    continue

                if not self.add_data_increment_to_apo_mmcif():
                    continue

            else:
                if not self.mmcif_files_can_be_replaced(xtal):
                    continue

                if not self.data_template_dict_exists(xtal):
                    continue

                if not self.db_dict_exists(xtal):
                    continue

#                if not self.zenodo_dict_exists(xtal):
#                    continue

                if not self.refine_bound_exists(xtal):
                    continue

                if not self.refine_mtz_exists(xtal):
                    continue

                if not self.mtzFree_exists(xtal):
                    continue

                if not self.aimless_logfile_exists(xtal):
                    continue

                if not self.ligand_in_pdb_file(xtal):
                    continue

                if not self.eventMTZ_exists((xtal)):
                    continue

                if not self.find_matching_event_map(xtal):
                    continue

                if not self.save_data_template_dict(xtal):
                    continue

                if not self.create_model_mmcif(xtal):
                    continue

                self.add_funding_information(xtal)

                if not self.add_ligand_cif_to_model_mmcif(xtal):
                    continue

                if not self.create_sf_mmcif(xtal):
                    continue

                if not self.event_maps_exist_in_sf_mmcif(xtal):
                    continue

                self.make_table_one(xtal)

        self.print_errorlist()
        self.Logfile.insert(
            '======= finished preparing mmcif files for wwPDB deposition =======')

    def ground_state_mmcif_exists(self):
        mmcifStatus = False
        for dirs in glob.glob(os.path.join(self.panddaDir, 'processed_datasets', '*')):
            xtal = dirs[dirs.rfind('/') + 1:]
            if os.path.isfile(os.path.join(dirs, xtal + '_sf.mmcif')):
                self.Logfile.insert(
                    '%s: found mmcif file for apo structure' % xtal)
                mmcifStatus = True

        return mmcifStatus

    def data_template_dict_exists(self, xtal):
        dictStatus = False
        self.data_template_dict = None
        self.Logfile.insert(
            '%s: reading information from depositTable for sample' % xtal)
        self.data_template_dict = self.db.get_deposit_dict_for_sample(xtal)
        if self.data_template_dict == {}:
            self.Logfile.error(
                '%s: cannot find data_template dictionary in depositTable; moving to next dataset...' % xtal)
            self.add_to_errorList(xtal)
        else:
            self.Logfile.insert(
                '%s: found data_template dictionary in depositTable' % xtal)
            dictStatus = True
        return dictStatus

    def update_beamline_info_data_template_dict(self, xtal):
        dls_beamlines = ['i02', 'i03', 'i04', 'i04-1', 'i23', 'i24']
        dls_beamline_dict = {'i02':      ['DIAMOND BEAMLINE I02',    'DECTRIS PILATUS 6M'],
                             'i03':      ['DIAMOND BEAMLINE I03',    'DECTRIS PILATUS 6M'],
                             'i04':      ['DIAMOND BEAMLINE I04',    'DECTRIS PILATUS 6M'],
                             'i04-1':    ['DIAMOND BEAMLINE I04-1',  'DECTRIS PILATUS 6M'],
                             'i23':      ['DIAMOND BEAMLINE I23',    'DECTRIS PILATUS 12M'],
                             'i24':      ['DIAMOND BEAMLINE I24',    'DECTRIS PILATUS 6M'], }

        if self.db_dict['DataCollectionBeamline'] in dls_beamlines:
            self.data_template_dict['radiation_source_type'] = dls_beamline_dict[self.db_dict['DataCollectionBeamline']][0]
            self.data_template_dict['radiation_detector_type'] = dls_beamline_dict[self.db_dict['DataCollectionBeamline']][1]
            self.data_template_dict['radiation_detector'] = 'PIXEL'
            self.data_template_dict['radiation_source'] = 'SYNCHROTRON'
            self.Logfile.insert(('%s: setting data collection beamline to %s' % (
                xtal, self.data_template_dict['radiation_source_type'])))

    def db_dict_exists(self, xtal):
        dictStatus = False
        self.db_dict = None
        self.Logfile.insert(
            '%s: reading information from mainTable for sample' % xtal)
        self.db_dict = self.db.get_db_dict_for_sample(xtal)
        if self.db_dict == {}:
            self.Logfile.error(
                '%s: cannot find db_dict dictionary in mainTable; moving to next dataset...' % xtal)
            self.add_to_errorList(xtal)
        else:
            self.Logfile.insert(
                '%s: found db_dict dictionary in mainTable' % xtal)
            self.update_beamline_info_data_template_dict(xtal)
            dictStatus = True
        return dictStatus

    def zenodo_dict_exists(self, xtal):
        dictStatus = False
        self.zenodo_dict = None
        if self.ground_state:
            self.zenodo_dict = self.db.get_zenodo_dict_for_pandda_analysis(
                self.panddaDir)
        else:
            self.Logfile.insert('%s: reading information from zenodoTable for pandda run: %s' % (
                xtal, self.db_dict['DimplePANDDApath']))
            self.zenodo_dict = self.db.get_zenodo_dict_for_pandda_analysis(
                self.db_dict['DimplePANDDApath'])
        if self.zenodo_dict == {}:
            dictStatus = True
            self.zenodo_dict['ZenodoDOI'] = ''
            self.Logfile.warning(
                '%s: cannot find information about zenodo deposition in zenodoTable!!!' % xtal)
#            self.Logfile.error('%s: cannot find information about zenodo deposition in zenodoTable; moving to next dataset...' %xtal)
#            self.add_to_errorList(xtal)
        else:
            self.Logfile.insert(
                '%s: found zenodo_dict dictionary in zenodoTable' % xtal)
            dictStatus = True
        return dictStatus

    def mmcif_files_can_be_replaced(self, xtal):
        status = True
        if self.overwrite_existing_mmcif:
            self.Logfile.insert(
                '%s: removing existing mmcif files as chosen by user' % xtal)
            self.db.execute_statement(
                "update depositTable set mmCIF_model_file='',mmCIF_SF_file='' where CrystalName is '{0!s}'".format(xtal))
            for mmcif in glob.glob('*.mmcif'):
                self.Logfile.warning('%s: removing %s' % (xtal, mmcif))
                os.system('/bin/rm ' + mmcif)
        else:
            for mmcif in glob.glob('*.mmcif'):
                self.Logfile.warning(
                    '%s: %s exists; skipping...' % (xtal, mmcif))
                status = False
        return status

    def refine_bound_exists(self, xtal):
        self.pdb = None
        self.Logfile.insert(
            '%s: checking if refine.split.bound-state.pdb exists' % xtal)
        fileStatus = False
        if os.path.isfile('refine.split.bound-state.pdb'):
            self.Logfile.insert(
                '%s: found refine.split.bound-state.pdb' % xtal)
            self.pdb = pdbtools('refine.split.bound-state.pdb')
            fileStatus = True
        else:
            self.Logfile.error(
                '%s: cannot find refine.split.bound-state.pdb; moving to next dataset...' % xtal)
            self.add_to_errorList(xtal)
        return fileStatus

    def refine_mtz_exists(self, xtal):
        self.mtz = None
        self.Logfile.insert('%s: checking if refine.mtz exists' % xtal)
        fileStatus = False
        if os.path.isfile('refine.mtz'):
            self.Logfile.insert('%s: found refine.mtz' % xtal)
            self.mtz = mtztools('refine.mtz')
            fileStatus = True
        else:
            self.Logfile.error(
                '%s: cannot find refine.mtz; moving to next dataset...' % xtal)
            self.add_to_errorList(xtal)
        return fileStatus

    def run_aimless_merge_only(self, xtal, unmerged, APpath):
        cmd = (
            'aimless hklin %s hklout mergeonly.mtz << eof > %s\n' % (unmerged, unmerged.replace('.mtz', '.log')) +
            ' onlymerge\n'
            'eof'
        )
        self.Logfile.insert('%s: running AIMLESS in onlymerge mode...' % xtal)
        os.system(cmd)
        os.chdir(os.path.join(self.projectDir, xtal))
        os.system('/bin/rm %s.log' % xtal)
        os.system('ln -s %s/%s %s.log' %
                  (APpath, unmerged.replace('.mtz', '.log'), xtal))
        self.Logfile.insert('%s: finished running AIMLESS' % xtal)
        if os.path.isfile(xtal+'.log'):
            self.Logfile.insert(
                '%s: AIMLESS logfile successfully created' % xtal)
        else:
            self.Logfile.error('%s: cannot find AIMLESS logfile...' % xtal)

    def prepare_aimless_log(self, xtal):
        parse().make_pseudo_aimless_log_from_json(xtal + '.log')

    def aimless_logfile_exists(self, xtal):
        self.Logfile.insert(
            '%s: checking if aimless logfile, i.e. %s.log, exists' % (xtal, xtal))
        fileStatus = False
        if os.path.isfile('%s.log' % xtal):
            self.Logfile.insert('%s: found %s.log' % (xtal, xtal))
            for n, line in enumerate(open('%s.log' % xtal)):
                if 'AIMLESS' in line:
                    fileStatus = True
                    break
#                if n < 10:
#                    print line
            if not fileStatus:
                if os.path.realpath('%s.log' % xtal).endswith('.table1'):
                    self.Logfile.warning(
                        '{0!s}: {1!s}.log seems to be a staraniso .table1 file'.format(xtal, xtal))
                    fileStatus = True
            if not fileStatus:
                self.Logfile.warning(
                    '%s: this does not seem to be an AIMLESS logfile' % xtal)
                Filepath = os.path.relpath(os.path.realpath('%s.log' % xtal))
                APpath = Filepath[:Filepath.rfind('/')]
                self.Logfile.insert(
                    '%s: relative path to logfile %s' % (xtal, APpath))
                self.Logfile.insert(
                    '%s: file path to logfile %s' % (xtal, Filepath))
                if os.path.isdir(APpath):
                    os.chdir(APpath)
                    for unmerged in glob.glob('*_scaled_unmerged.mtz'):
                        self.Logfile.insert('%s: found %s in %s' %
                                            (xtal, unmerged, APpath))
                        self.run_aimless_merge_only(xtal, unmerged, APpath)
                        fileStatus = True
                        break
                else:
                    self.Logfile.error(
                        '%s: %s is not a directory' % (xtal, APpath))
                    self.add_to_errorList(xtal)
                if not fileStatus:
                    os.chdir(os.path.join(self.projectDir, xtal))
                    self.Logfile.insert(
                        '%s: trying to prepare a pseudo-aimless file from json file...' % xtal)
                    self.prepare_aimless_log(xtal)
                    if os.path.isfile('aimless_dials.log'):
                        self.Logfile.insert(
                            '%s: found aimless_dials.log' % xtal)
                        fileStatus = True
                if not fileStatus:
                    self.Logfile.error(
                        '%s: cannot find a suitable AIMLESS logfile' % xtal)
                    self.add_to_errorList(xtal)

#            fileStatus = True
        else:
            self.Logfile.error(
                '%s: cannot find %s.log; moving to next dataset...' % (xtal, xtal))
            self.add_to_errorList(xtal)
        return fileStatus

    def mtzFree_exists(self, xtal):
        self.Logfile.insert(
            '%s: checking if %s.free.mtz exists' % (xtal, xtal))
        fileStatus = False
        if os.path.isfile('%s.free.mtz' % xtal):
            self.Logfile.insert('%s: found %s.free.mtz' % (xtal, xtal))
            fileStatus = True
        else:
            self.Logfile.error(
                '%s: cannot find %s.free.mtz; moving to next dataset...' % (xtal, xtal))
            self.add_to_errorList(xtal)
        return fileStatus

    def ligand_in_pdb_file(self, xtal):
        self.Logfile.insert(
            '%s: checking if refine.split.bound-state.pdb contains ligands of type LIG' % xtal)
        ligandStatus = False
        ligList = pdbtools(
            'refine.split.bound-state.pdb').get_residues_with_resname('LIG')
        if ligList is []:
            self.Logfile.error(
                '%s: refine.split.bound-state.pdb does not contain any modelled ligands of type LIG' % xtal)
            self.add_to_errorList(xtal)
        else:
            self.Logfile.insert(xtal + ': found ' +
                                str(len(ligList)) + ' ligands of type LIG')
            ligandStatus = True
        return ligandStatus

    def eventMTZ_exists(self, xtal):
        self.Logfile.insert('%s: checking if mtz of event maps exists' % xtal)
        eventMTZlist = []
        eventMTZexists = False
        if os.path.isfile('no_pandda_analysis_performed') or self.ignore_event_map:
            self.Logfile.warning(
                '%s: found empty file named "no_pandda_analysis_performed" which suggests we will ignore event maps for this sample' % xtal)
            eventMTZexists = True
        elif self.ignore_event_map:
            self.Logfile.warning(
                '%s: user selected to not include event map in SF mmcif file' % xtal)
            eventMTZexists = True
        else:
            for mtz in glob.glob('*event*.native*.mtz'):
                eventMTZlist.append(mtz[mtz.rfind('/')+1:])
            if eventMTZlist is []:
                self.Logfile.error(
                    '%s: MTZ files of event maps do not exists! Go to PANDDA tab and run "Event Map -> SF"' % xtal)
                self.add_to_errorList(xtal)
            else:
                self.Logfile.insert(
                    xtal + ': found ' + str(len(eventMTZlist)) + ' MTZ files of event maps')
                eventMTZexists = True
        return eventMTZexists

    def find_matching_event_map(self, xtal):
        self.eventList = []
        self.Logfile.insert(
            '%s: trying to find fitting event maps for modelled ligands' % xtal)
        ligList = self.pdb.save_residues_with_resname(
            os.path.join(self.projectDir, xtal), 'LIG')
        foundMatchingMap = None
        if os.path.isfile('no_pandda_analysis_performed') or self.ignore_event_map:
            self.Logfile.warning(
                '%s: found empty file named "no_pandda_analysis_performed" which suggests we will ignore event maps for this sample' % xtal)
            foundMatchingMap = True
            ligList = []

        self.Logfile.insert(
            '%s: looking for event maps for the following ligands -> %s' % (xtal, str(ligList)))
        for lig in sorted(ligList):
            ligID = lig.replace('.pdb', '')
#            if os.path.isfile('no_pandda_analysis_performed'):
#                self.Logfile.warning('%s: no pandda analysis performed; skipping this step...' %xtal)
#                foundMatchingMap = True
#                break
            ligCC = []

            # for newer BUSTER export
#            print '>>>> %s-event_*.native_%s.mtz' %(xtal,lig.replace('.pdb',''))
            for mtz in glob.glob(('%s-event_*.native_%s.mtz' % (xtal, lig.replace('.pdb', '')))):
                self.Logfile.insert(xtal + ': found ' + mtz)
                foundMatchingMap = True
                self.eventList.append(mtz)
                break

            if foundMatchingMap:
                continue

            for mtz in sorted(glob.glob('*event*.native.mtz')):
                self.get_lig_cc(xtal, mtz, lig)
                cc = self.check_lig_cc(mtz.replace('.mtz', '_CC'+ligID+'.log'))
                self.Logfile.insert('%s: %s -> CC = %s for %s' %
                                    (xtal, ligID, cc, mtz))
                try:
                    ligCC.append([mtz, float(cc)])
                except ValueError:
                    ligCC.append([mtz, 0.00])

            for mtz in sorted(glob.glob('*event*.native*P1.mtz')):
                self.get_lig_cc(xtal, mtz, lig)
                cc = self.check_lig_cc(mtz.replace('.mtz', '_CC'+ligID+'.log'))
                self.Logfile.insert('%s: %s -> CC = %s for %s' %
                                    (xtal, ligID, cc, mtz))
                try:
                    ligCC.append([mtz, float(cc)])
                except ValueError:
                    ligCC.append([mtz, 0.00])
            try:
                for cm in ligCC:
                    self.Logfile.insert('%s: cc = %s - %s' %
                                        (xtal, cm[1], cm[0]))
#                highestCC = max(ligCC, key=lambda x: x[0])[1]
                highestCCeventmap = max(ligCC, key=lambda x: float(x[1]))[0]
            except ValueError:
                highestCC = 0.00
                highestCCeventmap = None

#            if highestCC == 0.00 or ligCC is []:
            if highestCCeventmap is None or ligCC is []:
                self.Logfile.error(
                    '%s: best CC of ligand %s for any event map is 0!' % (xtal, lig))
                self.add_to_errorList(xtal)
                foundMatchingMap = False
            else:
                self.Logfile.insert('%s: selected event map for ligand %s is %s' % (
                    xtal, lig, highestCCeventmap))
                if os.path.isfile(highestCCeventmap.replace('.mtz', '_'+ligID+'.mtz')):
                    self.Logfile.warning('%s: symlink exists %s' % (
                        xtal, highestCCeventmap.replace('.mtz', '_'+ligID+'.mtz')))
                else:
                    self.Logfile.insert('%s: making symlink %s' % (
                        xtal, highestCCeventmap.replace('.mtz', '_'+ligID+'.mtz')))
                    os.system('ln -s %s %s' % (highestCCeventmap,
                                               highestCCeventmap.replace('.mtz', '_'+ligID+'.mtz')))
                if highestCCeventmap not in self.eventList:
                    self.eventList.append(highestCCeventmap)
                if foundMatchingMap is None:
                    foundMatchingMap = True
        return foundMatchingMap

    def get_lig_cc(self, xtal, mtz, lig):
        ligID = lig.replace('.pdb', '')
        self.Logfile.insert('%s: calculating CC for %s in %s' %
                            (xtal, lig, mtz))
        if os.path.isfile(mtz.replace('.mtz', '_CC'+ligID+'.log')):
            self.Logfile.warning('logfile of CC analysis exists; skipping...')
            return
        cmd = ('module load phenix\n'
               'phenix.get_cc_mtz_pdb %s %s > %s' % (mtz, lig, mtz.replace('.mtz', '_CC'+ligID+'.log')))
        os.system(cmd)

    def check_lig_cc(self, log):
        cc = 'n/a'
        if os.path.isfile(log):
            for line in open(log):
                if line.startswith('local'):
                    cc = line.split()[len(line.split()) - 1]
        else:
            self.Logfile.error('logfile does not exist: %s' % log)
        return cc

    def add_to_errorList(self, xtal):
        if xtal.replace(' ', '') == '':
            self.Logfile.warning(
                'trying to add xtal to error list, but xtal string is empty')
        else:
            if xtal not in self.errorList:
                self.errorList.append(xtal)

    def print_errorlist(self):
        if not self.errorList:
            self.Logfile.insert('XCE did not detect any problems during mmcif file preparation. '
                                'It is however recommended to check the logfile.')
        else:
            self.Logfile.warning('The following samples had problems during mmcif creation. '
                                 'Please check the logfile for details!')
            for xtal in self.errorList:
                self.Logfile.error(xtal)

    def save_data_template_dict(self, xtal):
        # check if file exists
        noError = True
        self.Logfile.insert('%s: preparing data_template.cif file' % xtal)
        if self.overwrite_existing_mmcif:
            self.data_template_dict['radiation_wavelengths'] = self.mtz.get_wavelength(
            )
            if str(self.data_template_dict['radiation_wavelengths']).startswith('0.0'):
                self.Logfile.error('%s: this does not seem to be the true experimental wavelength: %s' % (
                    xtal, str(self.data_template_dict['radiation_wavelengths'])))
                self.Logfile.insert(
                    '%s: trying to find it from %s.free.mtz...' % (xtal, xtal))
                if os.path.isfile('%s.free.mtz' % xtal):
                    self.data_template_dict['radiation_wavelengths'] = mtztools(
                        xtal+'.free.mtz').get_wavelength()
                    self.Logfile.warning('%s: found the following wavelength -> %s' % (
                        xtal, str(self.data_template_dict['radiation_wavelengths'])))
            self.Logfile.insert('%s: experimental wavelength according to %s is %s' % (
                xtal, self.mtz, self.data_template_dict['radiation_wavelengths']))
            if self.ground_state:
                os.chdir(self.projectDir)
#                self.data_template_dict['radiation_wavelengths'] = '1.000'
                self.data_template_dict['group_type'] = 'ground state'
                self.data_template_dict['group_title'] = 'PanDDA analysis group deposition of ground-state model'
                self.data_template_dict['group_description'] = self.data_template_dict['group_description'].replace('$ProteinName',
                                                                                                                    self.data_template_dict[
                                                                                                                        'Source_organism_gene'])
                self.data_template_dict['title'] = self.data_template_dict['structure_title_apo'].replace('$ProteinName',
                                                                                                          self.data_template_dict[
                                                                                                              'Source_organism_gene'])
            else:
                os.chdir(os.path.join(self.projectDir, xtal))
                # edit wavelength
#                self.data_template_dict['radiation_wavelengths'] = self.mtz.get_wavelength()

                title = self.data_template_dict['structure_title'].replace('$ProteinName', self.data_template_dict[
                    'Source_organism_gene']).replace('$CompoundName', self.db_dict['CompoundCode']).replace('($SampleID)', '('+xtal+')')

                self.data_template_dict['group_type'] = 'changed state'

                # edit title
                self.data_template_dict['group_title'] = self.data_template_dict['group_deposition_title'].replace('$ProteinName',
                                                                                                                   self.data_template_dict[
                                                                                                                       'Source_organism_gene']).replace(
                    '$CompoundName', self.db_dict['CompoundCode'])

#            self.data_template_dict[
#                    'group_title'] = 'PanDDA analysis group deposition of models with modelled events (e.g. bound ligands)'
                self.data_template_dict['group_description'] = self.data_template_dict['group_description'].replace('$ProteinName',
                                                                                                                    self.data_template_dict[
                                                                                                                        'Source_organism_gene'])
                self.data_template_dict['title'] = self.data_template_dict['group_title'] + ' -- ' + title

                if ('$ProteinName' or '$CompoundName') in self.data_template_dict['title']:
                    self.Logfile.error(
                        '%s: data_template - title not correctly formatted')
                    self.add_to_errorList(xtal)
                    noError = False

            # mutations
            mutations = self.data_template_dict['fragment_name_one_specific_mutation']
            if mutations.lower().replace(' ', '').replace('none', '').replace('null', '') == '':
                self.data_template_dict['fragment_name_one_specific_mutation'] = '?'
            else:
                self.data_template_dict['fragment_name_one_specific_mutation'] = '"' + \
                    mutations.replace(' ', '') + '"'

            # get protein chains
            self.data_template_dict['protein_chains'] = ''
            chains = self.pdb.GetProteinChains()
            for item in chains:
                self.data_template_dict['protein_chains'] += item + ','
            self.data_template_dict['protein_chains'] = self.data_template_dict['protein_chains'][:-1]

            data_template = templates().data_template_cif(self.data_template_dict)
#            site_details = self.make_site_description(xtal)
#            data_template += site_details
            if self.ground_state:
                f = open(os.path.join(self.projectDir,
                                      'data_template.cif'), 'w')
            else:
                f = open(os.path.join(self.projectDir,
                                      xtal, 'data_template.cif'), 'w')

            f.write(data_template)
            f.close()

        return noError

    def create_model_mmcif(self, xtal):
        fileStatus = False
        if self.ground_state:
            os.chdir(os.path.join(self.projectDir))
        else:
            os.chdir(os.path.join(self.projectDir, xtal))
        refSoft = self.pdb.get_refinement_program()

        if os.path.isdir('/dls'):
            #            if os.path.isfile(self.db_dict['RefinementMMCIFmodel_latest']):
            #                pdb_extract_init = 'source /dls/science/groups/i04-1/software/pdb-extract-v3.26/pdb-extract-v3.26-prod-src/setup.sh\n'
            #                pdb_extract_init += '/dls/science/groups/i04-1/software/pdb-extract-v3.26/pdb-extract-v3.26-prod-src/bin/pdb_extract'
            #            else:
            pdb_extract_init = 'source /dls/science/groups/i04-1/software/pdb-extract-prod/setup.sh\n'
            pdb_extract_init += '/dls/science/groups/i04-1/software/pdb-extract-prod/bin/pdb_extract'
        else:
            if os.path.isfile(os.path.join(os.getenv('XChemExplorer_DIR'), 'pdb_extract/pdb-extract-prod/bin/pdb_extract')):
                pdb_extract_init = 'source ' + \
                    os.path.join(os.getenv('XChemExplorer_DIR'),
                                 'pdb_extract/pdb-extract-prod/setup.sh') + '\n'
                pdb_extract_init += os.path.join(
                    os.getenv('XChemExplorer_DIR'), 'pdb_extract/pdb-extract-prod/bin/pdb_extract')
            else:
                pdb_extract_init = 'pdb_extract'

        if self.ground_state:
            refXtal = self.ground_state_pdb.split(
                '/')[len(self.ground_state_pdb.split('/')) - 2]
            self.Logfile.insert(
                'ground_state deposition reference dataset: {0!s}'.format(refXtal))
            aimless = os.path.join(self.logDir, refXtal, refXtal+'.log')
        else:
            aimless = '%s.log' % xtal

        isAimlessFile = False
        for n, line in enumerate(open(aimless)):
            if 'AIMLESS' in line:
                isAimlessFile = True
                break
        if not isAimlessFile:
            self.Logfile.warning(
                'processing log file does not seem to be an aimless file: {0!s}'.format(aimless))
            self.Logfile.insert('realpath of processing logfile: {0!s}'.format(
                os.path.realpath(aimless)))
            if os.path.realpath(aimless).endswith('.table1'):
                self.Logfile.warning(
                    '{0!s}: {1!s}.log seems to be a staraniso .table1 file'.format(xtal, xtal))
                isAimlessFile = True
            self.Logfile.warning(
                'it does not seem to originate from staraniso either')
        if not isAimlessFile:
            if os.path.isfile('aimless_dials.log'):
                aimless = 'aimless_dials.log'
            else:
                parse().make_pseudo_aimless_log_from_json(aimless)
                aimless = 'aimless_dials.log'

#        aimless = None
#        for n,line in enumerate(open('%s.log' %xtal)):
#            if 'AIMLESS' in line:
#                aimless = '%s.log' %xtal
#                break
#        if not aimless:
#            if os.path.isfile('aimless_dials.log'):
#                aimless = 'aimless_dials.log'

        if self.ground_state:
            #            refXtal = self.ground_state_pdb.split('/')[len(self.ground_state_pdb.split('/')) - 2]
            #            aimless = os.path.join(self.logDir,refXtal,refXtal+'.log')
            self.Logfile.insert('aimless.log file: ' + aimless)
            Cmd = (pdb_extract_init +
                   ' -r {0!s}'.format(refSoft) +
                   ' -iPDB {0!s}'.format(self.ground_state_pdb) +
                   ' -e MR'
                   ' -s AIMLESS'
                   ' -iLOG {0!s}'.format(aimless) +
                   ' -iENT data_template.cif'
                   ' -o {0!s}.mmcif > {1!s}.mmcif.log'.format(xtal, xtal))
        else:
            #            if os.path.isfile(self.db_dict['RefinementMMCIFmodel_latest']):
            #                self.Logfile.insert('%s: found MMCIF file; will use instead of PDB: %s' %(xtal,self.db_dict['RefinementMMCIFmodel_latest']))
            #                Cmd = (pdb_extract_init +
            #                       ' -r {0!s}'.format(refSoft) +
            #                       ' -iCIF {0!s}'.format(self.db_dict['RefinementMMCIFmodel_latest']) +
            #                       ' -e MR'
            #                       ' -s AIMLESS'
            #                       ' -iLOG {0!s}.log'.format(xtal) +
            #                       ' -iENT data_template.cif'
            #                       ' -o {0!s}.mmcif > {1!s}.mmcif.log'.format(xtal, xtal))
            #            else:
            Cmd = (pdb_extract_init +
                   ' -r {0!s}'.format(refSoft) +
                   ' -iPDB {0!s}'.format('refine.split.bound-state.pdb') +
                   ' -e MR'
                   ' -s AIMLESS'
                   ' -iLOG {0!s}'.format(aimless) +
                   ' -iENT data_template.cif'
                   ' -o {0!s}.mmcif > {1!s}.mmcif.log'.format(xtal, xtal))

        self.Logfile.insert(xtal + ': running pdb_extract: ' + Cmd)
        os.system(Cmd)

        self.update_model_mmcif_header(xtal)

        if os.path.isfile(xtal+'.mmcif') and os.path.getsize(xtal+'.mmcif') > 20000:
            self.Logfile.insert(
                '%s: model mmcif file successfully created' % xtal)
            if self.ground_state:
                self.db.execute_statement(
                    "update depositTable set mmCIF_model_file='{0!s}.mmcif' where CrystalName is '{1!s}' and DimplePANDDApath is '{2!s}'".format(xtal,
                                                                                                                                                 xtal, self.panddaDir))
            else:
                self.db.execute_statement(
                    "update depositTable set mmCIF_model_file='{0!s}.mmcif' where CrystalName is '{1!s}'".format(xtal, xtal))
            fileStatus = True
        else:
            self.Logfile.error(
                '%s: model mmcif file was not created successfully')
            self.add_to_errorList(xtal)

        return fileStatus

    def update_model_mmcif_header(self, xtal):
        self.Logfile.insert('%s: updating header of model mmcif file' % xtal)
        foundSoftwareBlock = False
        amendSoftwareBlock = False
        softwareEntry = []
        for i, line in enumerate(fileinput.input(xtal + '.mmcif', inplace=1)):
            #            if i == 4: sys.stdout.write('\n')  # write a blank line after the 5th line
            if '_software.pdbx_ordinal' in line:
                foundSoftwareBlock = True
            if foundSoftwareBlock:
                if not line.startswith('_'):
                    try:
                        softwareEntry.append(int(line.split()[0]))
                    except (ValueError, IndexError):
                        pass
                if '#' in line:
                    amendSoftwareBlock = True
                    foundSoftwareBlock = False
            if '_refine.pdbx_ls_cross_valid_method' in line:
                sys.stdout.write(
                    '_refine.pdbx_ls_cross_valid_method               THROUGHOUT \n')

            elif '_refine.pdbx_starting_model' in line:
                sys.stdout.write('_refine.pdbx_starting_model                      {0!s} \n'.format(
                    self.data_template_dict['pdbx_starting_model']))

            elif '_refine.pdbx_method_to_determine_struct' in line:
                sys.stdout.write(
                    "_refine.pdbx_method_to_determine_struct          'FOURIER SYNTHESIS'\n")
            elif '_struct.title      ---' in line:
                #                self.Logfile.warning('structure title was not created by pdb_extract; there might be an issue with REFMAC 5.8.0238')
                #                self.Logfile.insert('trying to get title from data_template.cif file')
                Title = ''
                foundTitle = False
                for li in open('data_template.cif'):
                    if li.startswith('_struct.title'):
                        foundTitle = True
                    if foundTitle:
                        if li.replace(' ', '').replace('\n', '').replace('\r', '') == ';':
                            Title += li
                            break
                        Title += li
                sys.stdout.write(Title)
            elif amendSoftwareBlock:
                cifItem = (
                    "{0!s} {1!s} ? ? program ? ? 'data reduction' ? ?\n".format(str(max(softwareEntry) + 1),
                                                                                self.data_template_dict[
                        'data_integration_software']) +
                    '{0!s} {1!s} ? ? program ? ? phasing ? ?\n'.format(str(max(softwareEntry) + 2),
                                                                       self.data_template_dict[
                        'phasing_software'])
                    #                    '#\n'
                    #                    "loop_\n"
                    #                    "_pdbx_related_exp_data_set.ordinal\n"
                    #                    "_pdbx_related_exp_data_set.data_reference\n"
                    #                    "_pdbx_related_exp_data_set.metadata_reference\n"
                    #                    "_pdbx_related_exp_data_set.data_set_type\n"
                    #                    "_pdbx_related_exp_data_set.details\n"
                    #                    " 1  "
                    #                    " 'doi:%s' "  %self.zenodo_dict['ZenodoDOI']+
                    #                    " 'doi:%s' "  %self.zenodo_dict['ZenodoDOI']+
                    #                    " 'other data'  "                 # 'other data' is only thing wwPDB accepts at the moment
                    #                    " 'Complete PanDDA analysis'\n"
                    #                    "\n"
                )

                sys.stdout.write(cifItem)
                amendSoftwareBlock = False

            else:
                sys.stdout.write(line)


#            elif '_reflns.d_resolution_low' in line and len(line.split()) == 2:
#                sys.stdout.write('_reflns.d_resolution_low             {0!s}\n'.format(min(low_reso_list)))
#
#            elif '_refine.ls_d_res_low' in line and len(line.split()) == 2:
#                sys.stdout.write('_refine.ls_d_res_low                             {0!s}\n'.format(min(low_reso_list)))

    def add_funding_information(self, xtal):
        pdbx_funding_ordinal_one = self.data_template_dict['pdbx_funding_ordinal_one']
        if pdbx_funding_ordinal_one.lower().replace(' ', '').replace('none', '').replace('null', '') == '':
            self.data_template_dict['pdbx_funding_ordinal_one'] = ''
            self.data_template_dict['pdbx_funding_organization_one'] = ''
            self.data_template_dict['pdbx_grant_number_one'] = ''
            self.data_template_dict['pdbx_grant_country_one'] = ''
            funding_one = ''
        else:
            funding_one = "%s '%s' '%s' '%s'\n" % (self.data_template_dict['pdbx_funding_ordinal_one'], self.data_template_dict[
                                                   'pdbx_funding_organization_one'], self.data_template_dict['pdbx_grant_number_one'], self.data_template_dict['pdbx_grant_country_one'])

        pdbx_funding_ordinal_two = self.data_template_dict['pdbx_funding_ordinal_two']
        if pdbx_funding_ordinal_two.lower().replace(' ', '').replace('none', '').replace('null', '') == '':
            self.data_template_dict['pdbx_funding_ordinal_two'] = ''
            self.data_template_dict['pdbx_funding_organization_two'] = ''
            self.data_template_dict['pdbx_grant_number_two'] = ''
            self.data_template_dict['pdbx_grant_country_two'] = ''
            funding_two = ''
        else:
            funding_two = "%s '%s' '%s' '%s'\n" % (self.data_template_dict['pdbx_funding_ordinal_two'], self.data_template_dict[
                                                   'pdbx_funding_organization_two'], self.data_template_dict['pdbx_grant_number_two'], self.data_template_dict['pdbx_grant_country_two'])

        pdbx_funding_ordinal_three = self.data_template_dict['pdbx_funding_ordinal_three']
        if pdbx_funding_ordinal_three.lower().replace(' ', '').replace('none', '').replace('null', '') == '':
            self.data_template_dict['pdbx_funding_ordinal_three'] = ''
            self.data_template_dict['pdbx_funding_organization_three'] = ''
            self.data_template_dict['pdbx_grant_number_three'] = ''
            self.data_template_dict['pdbx_grant_country_three'] = ''
            funding_three = ''
        else:
            funding_three = "%s '%s' '%s' '%s'\n" % (self.data_template_dict['pdbx_funding_ordinal_three'], self.data_template_dict[
                                                     'pdbx_funding_organization_three'], self.data_template_dict['pdbx_grant_number_three'], self.data_template_dict['pdbx_grant_country_three'])

        funding_info = (
            '#\n'
            'loop_\n'
            '_pdbx_audit_support.ordinal                   \n'
            '_pdbx_audit_support.funding_organization      \n'
            '_pdbx_audit_support.grant_number              \n'
            '_pdbx_audit_support.country                   \n'
            + funding_one + funding_two + funding_three +
            '#\n'
        )

        f = open(xtal + '.mmcif', 'a')
        f.write(funding_info)
        f.close()

    def add_ligand_cif_to_model_mmcif(self, xtal):
        filestatus = False
        self.Logfile.insert('%s: looking for ligand restraints file...' % xtal)
        os.chdir(os.path.join(self.projectDir, xtal))
        if os.path.isfile(self.db_dict['RefinementMMCIFmodel_latest']):
            self.Logfile.insert('%s: found %s; assuming that ligand cif dictionary is already included...' % (
                xtal, self.db_dict['RefinementMMCIFmodel_latest']))
            filestatus = True
        else:
            if os.path.isfile(self.db_dict['CompoundCode']+'.cif'):
                self.Logfile.insert('%s: found ligand restraints file -> %s' %
                                    (xtal, self.db_dict['CompoundCode']+'.cif'))
                self.Logfile.insert(
                    '%s: adding ligand restraints file to model mmcif' % xtal)
                cif = ''
                for line in open(self.db_dict['CompoundCode']+'.cif'):
                    cif += line
                f = open(xtal+'.mmcif', 'a')
                f.write(cif)
                f.close()
                filestatus = True
            else:
                self.Logfile.warning('%s: could not find %s' % (
                    xtal, self.db_dict['CompoundCode']+'.cif'))
        return filestatus

    def make_table_one(self, xtal):
        os.chdir(os.path.join(self.projectDir, xtal))
        if os.path.isfile(xtal + '.mmcif') and os.path.getsize(xtal + '.mmcif') > 20000:
            self.Logfile.insert('making table_1 for %s.mmcif' % xtal)
            if os.path.isdir('/dls'):
                extract_table_init = 'source /dls/science/groups/i04-1/software/pdb-extract-prod/setup.sh\n'
                extract_table_init += '/dls/science/groups/i04-1/software/pdb-extract-prod/bin/extract_table'
            else:
                extract_table_init = 'source ' + os.path.join(os.getenv('XChemExplorer_DIR'),
                                                              'pdb_extract/pdb-extract-prod/setup.sh') + '\n'
                extract_table_init += os.path.join(os.getenv('XChemExplorer_DIR'),
                                                   'pdb_extract/pdb-extract-prod/bin/extract_table')

            Cmd = extract_table_init + ' ' + xtal + '.mmcif'

            self.Logfile.insert(xtal + ': running sf_convert: ' + Cmd)
            os.system(Cmd)

            if os.path.isfile('cryst-table-1.out'):
                os.system('/bin/mv cryst-table-1.out %s-table-1.txt' % xtal)
                self.Logfile.insert(
                    '%s: table_1 successfully created; updating database...' % xtal)
                self.db.execute_statement(
                    "update mainTable set table_one='{0!s}-table-1.txt' where CrystalName is '{1!s}'".format(xtal, xtal))
            else:
                self.Logfile.warning('%s: could not create table_1' % xtal)

    def create_sf_mmcif(self, xtal):
        fileStatus = False
        if self.ground_state:
            os.chdir(self.projectDir)
        else:
            os.chdir(os.path.join(self.projectDir, xtal))

        if os.path.isfile('no_pandda_analysis_performed') or self.ignore_event_map:
            mtzin = 'refine.mtz '
        else:
            mtzin = 'refine.mtz ' + xtal + '.free.mtz '
            for event in self.eventList:
                mtzin += event + ' '

#        if self.ground_state:
#            mtzin = self.ground_state_mean_mtz

        if os.path.isdir('/dls'):
            pdb_extract_init = 'source /dls/science/groups/i04-1/software/pdb-extract-prod/setup.sh\n'
            pdb_extract_init += '/dls/science/groups/i04-1/software/pdb-extract-prod/bin/sf_convert'
        else:
            if os.path.isfile(os.path.join(os.getenv('XChemExplorer_DIR'), 'pdb_extract/pdb-extract-prod/bin/sf_convert')):
                pdb_extract_init = 'source ' + os.path.join(os.getenv('XChemExplorer_DIR'),
                                                            'pdb_extract/pdb-extract-prod/setup.sh') + '\n'
                pdb_extract_init += os.path.join(
                    os.getenv('XChemExplorer_DIR'), 'pdb_extract/pdb-extract-prod/bin/sf_convert')
            else:
                pdb_extract_init = 'sf_convert'

        Cmd = (pdb_extract_init +
               ' -o mmcif'
               ' -sf %s' % mtzin +
               ' -out {0!s}_sf.mmcif  > {1!s}.sf_mmcif.log'.format(xtal, xtal))

        self.Logfile.insert(xtal + ': running sf_convert: ' + Cmd)
        os.system(Cmd)
        os.system(
            '/bin/rm sf_format_guess.text mtzdmp.log SF_4_validate.cif sf_information.cif')

        self.update_sf_mmcif_file(xtal)

        if os.path.isfile(xtal+'_sf.mmcif') and os.path.getsize(xtal+'_sf.mmcif') > 20000:
            self.Logfile.insert(
                '%s: SF mmcif file successfully created' % xtal)
            if self.ground_state:
                self.db.execute_statement(
                    "update depositTable set mmCIF_SF_file='{0!s}_sf.mmcif' where CrystalName is '{1!s}' and DimplePANDDApath is '{2!s}'".format(xtal, xtal, self.panddaDir))
            else:
                self.db.execute_statement(
                    "update depositTable set mmCIF_SF_file='{0!s}_sf.mmcif' where CrystalName is '{1!s}'".format(xtal, xtal))
            fileStatus = True
        else:
            self.Logfile.error(
                '%s: SF mmcif file was not created successfully')
            self.add_to_errorList(xtal)

        return fileStatus

    def apo_mmcif_exists(self):
        fileStatus = False
        self.Logfile.insert('checking if mmcif files of apo structures exist')
        counter = 0
        for mmcif in glob.glob(os.path.join(self.panddaDir, 'processed_datasets', '*', '*.mmcif')):
            if os.path.isfile(mmcif):
                counter += 1
        if counter < 40:
            self.Logfile.error('found only %s apo mmcif files' % str(counter))
            self.Logfile.warning(
                'you may need to run "PanDDA tab"/"apo -> mmcif"')
        else:
            self.Logfile.insert(
                'found %s apo mmcif files; seems OK!' % str(counter))
            fileStatus = True
        return fileStatus

    def add_apo_sf_mmcif_to_ground_state_mmcif(self):
        os.chdir(self.projectDir)
        self.Logfile.insert(
            'checking pandda directory for apo mmcif files: '+self.panddaDir)
        f = open('ground_state_sf.mmcif', 'w')

        refXtal = self.ground_state_pdb.split(
            '/')[len(self.ground_state_pdb.split('/'))-2]
        # make sure that mmcof belonging to ref PDB file is first
        xtalList = [refXtal]
        for dirs in glob.glob(os.path.join(self.panddaDir, 'processed_datasets', '*')):
            xtal = dirs[dirs.rfind('/') + 1:]
            if xtal not in xtalList:
                xtalList.append(xtal)

        counter = 1
        for xtal in xtalList:
            # this is needed in case single files are in processed_datasets
            if not os.path.isdir(os.path.join(self.panddaDir, 'processed_datasets', xtal)):
                continue
            else:
                dirs = os.path.join(self.panddaDir, 'processed_datasets', xtal)
            self.Logfile.insert(
                '%s: reading saoked compound information from database' % xtal)
            xtalDict = self.db.get_db_dict_for_sample(xtal)
            if xtalDict['CompoundSMILES'].lower().replace(' ', '') == '':
                smiles = 'none'
            elif 'none' in xtalDict['CompoundSMILES'].lower().replace(' ', ''):
                smiles = 'none'
            elif 'null' in xtalDict['CompoundSMILES'].lower().replace(' ', ''):
                smiles = 'none'
            else:
                smiles = xtalDict['CompoundSMILES'].replace(' ', '')
            self.Logfile.insert('%s: compound SMILES -> %s' % (xtal, smiles))
            if os.path.isfile(os.path.join(dirs, xtal+'_sf.mmcif')):
                self.Logfile.insert(
                    'adding %s_sf.mmcif to ground-state_sf.mmcif' % xtal)
                for line in open(os.path.join(dirs, xtal+'_sf.mmcif')):
                    if line.startswith('_cell.angle_gamma'):
                        newLine = line
                        newLine += '#\n'
                        newLine += '_diffrn.id                  1\n'
                        newLine += '_diffrn.details             "diffraction data from crystal %s; soaked compound: %s"\n' % (
                            str(counter), smiles.replace('\n', '').replace('\r', ''))
                        f.write(newLine)
                        counter += 1
                    else:
                        f.write(line)
        f.close()
        self.Logfile.insert(
            'added %s apo mmcif files to ground-state mmcif' % str(counter))
        return True

    def add_data_increment_to_apo_mmcif(self):
        self.Logfile.insert(
            'inrementing data_rxxxxsf in ground-state_sf.mmcif')
        x = ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
             'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
        a = 0
        b = 0
        c = 0

        foundFirstLine = False
        foundCulprit = False
        datasetCounter = 0
        if os.path.isfile(os.path.join(self.panddaDir, 'ground_state_sf.mmcif')):
            f = open('ground_state_sf_tmp.mmcif', 'w')
            for n, line in enumerate(open(os.path.join(self.panddaDir, 'ground_state_sf.mmcif'))):
                if line.startswith('data_rxxxxsf') and not foundFirstLine:
                    foundFirstLine = True
                    a += 1
                    f.write(line)
                elif line.startswith('data_rxxxxsf') and foundFirstLine:
                    if a == len(x):
                        a = 1
                        b += 1
                    if b == len(x):
                        a = 1
                        b = 1
                        c += 1
                    newLine = line.replace('xsf', 's%ssf' %
                                           str(x[a]+x[b]+x[c]))
                    datasetCounter += 1
                    f.write(newLine)
                    a += 1
                    self.Logfile.insert('new dataset block: %s -> %s' % (
                        str(datasetCounter), newLine.replace('\n', '').replace('\r', '')))
#                    if datasetCounter % 50 == 0:
#                        self.Logfile.insert('%s data_rxxxxsf records edited...' %str(datasetCounter))
                else:
                    f.write(line)
            f.close()
        os.chdir(self.panddaDir)
        os.system('/bin/mv ground_state_sf_tmp.mmcif ground_state_sf.mmcif')

        if os.path.isfile('ground_state_sf.mmcif') and os.path.getsize('ground_state_sf.mmcif') > 20000:
            self.Logfile.insert(
                'ground_state: SF mmcif file successfully created')
            self.db.execute_statement(
                "update depositTable set mmCIF_SF_file='ground_state_sf.mmcif' where CrystalName is 'ground_state' and DimplePANDDApath is '{0!s}'".format(self.panddaDir))
        else:
            self.Logfile.error(
                '%s: SF mmcif file was not created successfully')
            self.add_to_errorList(xtal)

        return True

    def event_maps_exist_in_sf_mmcif(self, xtal):
        fileOK = False
        # set to -2 since first two data blocks are initial/final.mtz and data.mtz
        n_eventMTZ_found = -2
        if os.path.isfile('no_pandda_analysis_performed') or self.ignore_event_map:
            self.Logfile.warning(
                '%s: no pandda analysis performed; skipping this step...' % xtal)
            fileOK = True
        else:
            for line in open(xtal + '_sf.mmcif'):
                if line.startswith('_refln.crystal_id'):
                    n_eventMTZ_found += 1
            if n_eventMTZ_found == len(self.eventList):
                fileOK = True
                self.Logfile.insert('%s: %s_sf.mmcif should contains %s of %s event maps' % (
                    xtal, xtal, n_eventMTZ_found, len(self.eventList)))
            else:
                self.Logfile.error('%s: %s_sf.mmcif should contains only %s of %s event maps' % (
                    xtal, xtal, n_eventMTZ_found, len(self.eventList)))
                self.add_to_errorList(xtal)
        return fileOK

    def update_sf_mmcif_file(self, xtal):
        self.Logfile.insert('%s: updating %s_sf.mmcif' % (xtal, xtal))

        if self.ground_state:
            bound = ["data for PanDDA ground-state-mean-map"]
        else:
            bound = ["data from final refinement with ligand, final.mtz",
                     "data from original reflections, data.mtz",
                     "data for ligand evidence map (PanDDA event map), event_map_$.mtz"]

        block = -1

        self.Logfile.insert('%s: reading wavelength from mtz file; lambda = %s' % (
            xtal, str(self.data_template_dict['radiation_wavelengths'])))

        if os.path.isfile('no_pandda_analysis_performed') or self.ignore_event_map:
            self.Logfile.warning(
                '%s: apparently not a pandda deposition; will skip this step...' % xtal)
            return None

        for i, line in enumerate(fileinput.input(xtal + '_sf.mmcif', inplace=1)):

            if line.startswith('_cell.length_a'):
                block += 1

            if line.startswith('_cell.angle_gamma'):
                if block >= 2:
                    n = 2
                else:
                    n = block
                sys.stdout.write(line)
                newLines = ('#\n'
                            '_diffrn.id                  1\n'
                            '_diffrn.details             "%s"\n' % bound[n]).replace('$', str(block-1))
                sys.stdout.write(newLines)
            elif line.startswith('_diffrn_radiation_wavelength.wavelength'):
                sys.stdout.write('_diffrn_radiation_wavelength.wavelength   {0!s}\n'.format(
                    str(self.data_template_dict['radiation_wavelengths'])))
            else:
                sys.stdout.write(line)


class prepare_for_group_deposition_upload(QtCore.QThread):

    def __init__(self, database, xce_logfile, depositDir, projectDir, type):
        QtCore.QThread.__init__(self)
        self.database = database
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.db = XChemDB.data_source(database)
        self.depositDir = depositDir
        self.projectDir = projectDir
        self.type = type

    def run(self):

        TextIndex = ''
        os.chdir(self.depositDir)

        # ligand bound structures
        if self.type == 'ligand_bound':
            self.Logfile.insert(
                'checking depositionTable for mmcif files of ligand-bound structures')
            depositList = self.db.execute_statement(
                "select CrystalName from mainTable where RefinementOutcome like '5%';")
            xtalString = '('
            for item in depositList:
                xtal = str(item[0])
                self.Logfile.insert(
                    '%s: adding mmcif files to final tar.bz2 file' % xtal)
                xtalString += "CrystalName = '"+xtal+"' or "
            xtalString = xtalString[:-4] + ')'
            toDeposit = self.db.execute_statement(
                "select CrystalName,mmCIF_model_file,mmCIF_SF_file,DimplePANDDApath from depositTable where StructureType is 'ligand_bound' and %s;" % xtalString)
        elif self.type == 'ground_state':
            self.Logfile.insert(
                'checking depositionTable for mmcif files of ground-state structures')
            toDeposit = self.db.execute_statement(
                "select CrystalName,mmCIF_model_file,mmCIF_SF_file,DimplePANDDApath from depositTable where StructureType is 'ground_state';")
        else:
            return

        for n, item in enumerate(sorted(toDeposit)):
            xtal = str(item[0])
            if self.type == 'ligand_bound':
                mmcif = os.path.join(self.projectDir, xtal, str(item[1]))
                mmcif_sf = os.path.join(self.projectDir, xtal, str(item[2]))
            elif self.type == 'ground_state':
                mmcif = os.path.join(str(item[3]), str(item[1]))
                mmcif_sf = os.path.join(str(item[3]), str(item[2]))
            else:
                continue
            self.Logfile.insert('%s: %s/ %s' % (xtal, mmcif, mmcif_sf))
            if os.path.isfile(mmcif) and os.path.isfile(mmcif_sf):
                self.Logfile.insert(
                    'copying {0!s} to {1!s}'.format(mmcif, self.depositDir))
                os.system('/bin/cp {0!s} .'.format(mmcif))
                if self.type == 'ground_state':
                    os.system(
                        '/bin/mv ground_state.mmcif ground_state_{0!s}.mmcif'.format(str(n)))
                    mmcif = mmcif.replace(
                        'ground_state.mmcif', 'ground_state_{0!s}.mmcif'.format(str(n)))
                self.Logfile.insert('copying {0!s} to {1!s}'.format(
                    mmcif_sf, self.depositDir))
                os.system('/bin/cp {0!s} .'.format(mmcif_sf))
                if self.type == 'ground_state':
                    os.system(
                        '/bin/mv ground_state_sf.mmcif ground_state_{0!s}_sf.mmcif'.format(str(n)))
                    mmcif_sf = mmcif_sf.replace(
                        'ground_state_sf.mmcif', 'ground_state_{0!s}_sf.mmcif'.format(str(n)))
            else:
                self.Logfile.error('cannot find mmcif file for '+xtal)

            text = ('label: {0!s}-{1!s}\n'.format(xtal, self.type) +
                    'description: {0!s} structure of {1!s}\n'.format(self.type, xtal) +
                    'model: {0!s}\n'.format(mmcif[mmcif.rfind('/')+1:]) +
                    'sf: {0!s}\n\n'.format(mmcif_sf[mmcif_sf.rfind('/')+1:]))
            TextIndex += text

        f = open('index.txt', 'w')
        f.write(TextIndex)
        f.close()

        # checking of tar.bz2 files exisit
        fileList = []
        for i in sorted(glob.glob('%s_structures.tar.bz2.*' % self.type)):
            fileList.append(int(i[i.rfind('.')+1:]))

        if os.path.isfile('%s_structures.tar.bz2' % self.type):
            if fileList == []:
                self.Logfile.warning(
                    'moving existing %s_structures.tar.bz2 to %s_structures.tar.bz2.1' % (self.type, self.type))
                os.system(
                    '/bin/mv %s_structures.tar.bz2 %s_structures.tar.bz2.1' % (self.type, self.type))
            else:
                self.Logfile.warning('moving existing %s_structures.tar.bz2 %s_structures.tar.bz2.%s' % (
                    self.type, self.type, str(max(fileList)+1)))
                os.system('/bin/mv %s_structures.tar.bz2 %s_structures.tar.bz2.%s' %
                          (self.type, self.type, str(max(fileList)+1)))

        self.Logfile.insert('preparing tar archive...')
        os.system(
            'tar -cvf {0!s}_structures.tar *mmcif index.txt'.format(self.type))
        self.Logfile.insert('bzipping archive...')
        os.system('bzip2 {0!s}_structures.tar'.format(self.type))
        self.Logfile.insert(
            'removing all bound mmcif files and index.txt file from '+self.depositDir)
        os.system('/bin/rm -f *mmcif index.txt')
        self.Logfile.insert('done!')


#        # apo structures
#        TextIndex=''
#        toDeposit=self.db.execute_statement("select CrystalName,mmCIF_model_file,mmCIF_SF_file from depositTable where StructureType is 'apo';")
#        for item in sorted(toDeposit):
#            xtal=str(item[0])
#            mmcif=str(item[1])
#            mmcif_sf=str(item[2])
#            if os.path.isfile(mmcif) and os.path.isfile(mmcif_sf):
#                self.Logfile.insert('copying {0!s} to {1!s}'.format(mmcif, self.depositDir))
#                os.system('/bin/cp {0!s} .'.format(mmcif))
#                self.Logfile.insert('copying {0!s} to {1!s}'.format(mmcif_sf, self.depositDir))
#                os.system('/bin/cp {0!s} .'.format(mmcif_sf))
#            else:
#                self.Logfile.insert('cannot find apo mmcif file for '+xtal+' => ERROR')
#
#            text = (    'label: {0!s}-apo\n'.format(xtal)+
#                        'description: apo structure of {0!s}\n'.format(xtal)+
#                        'model: {0!s}\n'.format(mmcif[mmcif.rfind('/')+1:])+
#                        'sf: {0!s}\n\n'.format(mmcif_sf[mmcif_sf.rfind('/')+1:])          )
#            TextIndex+=text
#
#        f = open('index.txt','w')
#        f.write(TextIndex)
#        f.close()
#
#        self.Logfile.insert('preparing tar archive...')
#        os.system('tar -cvf apo_structures.tar *apo* index.txt')
#        self.Logfile.insert('bzipping archive...')
#        os.system('bzip2 apo_structures.tar')
#        self.Logfile.insert('removing all apo mmcif files and index.txt file from '+self.depositDir)
#        os.system('/bin/rm -f *apo*mmcif index.txt')
#        self.Logfile.insert('done!')
#


class import_PDB_IDs(QtCore.QThread):
    def __init__(self, pdbCodes, database, xce_logfile):
        QtCore.QThread.__init__(self)
        self.pdbCodes = pdbCodes
        self.database = database
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.db = XChemDB.data_source(database)

    def run(self):

        for line in self.pdbCodes.split('\n'):
            if len(line.split('/')) == 2 and '-ligand_bound' in line:
                xtal = line[:line.rfind('-ligand_bound')].replace(' ', '')
                pdbID = line.split('/')[1].replace(' ', '')
                self.Logfile.insert('setting PDB ID for '+xtal+' to '+pdbID)
                sqlite = "UPDATE mainTable SET Deposition_PDB_ID='{0!s}',RefinementOutcome='6 - Deposited' where CrystalName is '{1!s}';".format(
                    pdbID, xtal)
                self.db.execute_statement(sqlite)


class compare_smiles_in_db_with_ligand_in_pdb(QtCore.QThread):
    def __init__(self, projectDir, database, xce_logfile):
        QtCore.QThread.__init__(self)
        self.projectDir = projectDir
        self.database = database
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.db = XChemDB.data_source(database)
        self.ErrorDict = {}

    def update_ErrorDict(self, xtal, message):
        if xtal not in self.ErrorDict:
            self.ErrorDict[xtal] = []
        self.ErrorDict[xtal].append(message)

    def run(self):

        os.chdir(self.projectDir)

        progress_step = 1
        if len(glob.glob('*')) != 0:
            progress_step = 100/float(len(glob.glob('*')))
        else:
            progress_step = 1
        progress = 0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        for xtal in sorted(glob.glob('*')):
            if os.path.isfile(os.path.join(xtal, 'refine.pdb')):
                smiles = self.db.execute_statement(
                    "select CompoundSmiles,CompoundCode from mainTable where CrystalName is '{0!s}'".format(xtal))
                try:
                    LigandSmiles = str(smiles[0][0])
                    LigandCode = str(smiles[0][1])
                    elementDict_smiles = smilestools(
                        LigandSmiles).ElementDict()
                except IndexError:
                    self.Logfile.error(
                        "{0!s}: something is seems to be wrong with the CompoundCode or SMILES string: {1!s}".format(xtal, str(smiles)))
                    continue

                pdb = pdbtools(os.path.join(xtal, 'refine.pdb'))
                ligandList = pdb.ligand_details_as_list()
                for ligand in ligandList:
                    resname = ligand[0]
                    chainID = ligand[1]
                    resseq = ligand[2]
                    altLoc = ligand[3]
                    elementDict_ligand = pdb.ElementDict(
                        resname, chainID, resseq, altLoc)
                    for element in elementDict_ligand:
                        if elementDict_ligand[element] != elementDict_smiles[element]:
                            self.Logfile.error('{0!s}: {1!s} {2!s} {3!s} {4!s} contains different number of atoms than smiles in DB: {5!s} -> {6!s}'.format(
                                xtal, resname, chainID, resseq, altLoc, LigandSmiles, LigandCode))
                            self.update_ErrorDict(xtal, '{0!s} {1!s} {2!s} {3!s} contains different number of atoms than smiles in DB'.format(
                                resname, chainID, resseq, altLoc))
                            break

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.emit(QtCore.SIGNAL('show_error_dict'), self.ErrorDict)
