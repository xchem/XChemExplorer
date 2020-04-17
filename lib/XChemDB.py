# last edited: 21/07/2017, 15:00

import sqlite3
import os,sys,glob
import csv
from datetime import datetime
import getpass

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
#from XChemUtils import parse
import XChemLog


class data_source:

    def __init__(self,data_source_file):

        self.data_source_file=data_source_file

        #   [column_name in DB, column_name shown in XCE, SQLite column type (Integer,Text,PKEY...)]
        self.column_list=[
            # SQLite column name                    XCE column name                             SQLite type             display in overview tab
            # from Lab36
            ['ID',                                   'ID',                                      'INTEGER PRIMARY KEY',  0],
            ['LabVisit',                             'LabVisit',                                'TEXT',                 1],
            ['LibraryPlate',                         'LibraryPlate',                            'TEXT',                 0],
            ['SourceWell',                           'SourceWell',                              'TEXT',                 0],
            ['LibraryName',                          'LibraryName',                             'TEXT',                 1],
            ['CompoundSMILES',                       'Smiles',                                  'TEXT',                 1],
            ['CompoundSMILESproduct',                'Smiles - Product',                        'TEXT',                 1],
            ['CompoundCode',                         'Compound ID',                             'TEXT',                 1],
            ['CompoundStereo',                       'Compound Stereo',                         'TEXT',                 1],
            ['CompoundStereoSMILES',                 'Stereo SMILES',                           'TEXT',                 1],
            ['CompoundStereoCIFprogram',             'Stereo CIF program',                      'TEXT',                 1],
            ['CompoundStereoCIFs',                   'Compound Stereo IDs',                     'TEXT',                 1],
            ['CompoundAutofitCC',                    'autofit CC',                              'TEXT',                 1],
            ['CompoundAutofitprogram',               'autofit program',                         'TEXT',                 1],
            ['CrystalPlate',                         'CrystalPlate',                            'TEXT',                 1],
            ['CrystalWell',                          'CrystalWell',                             'TEXT',                 1],
            ['EchoX',                                'EchoX',                                   'TEXT',                 0],
            ['EchoY',                                'EchoY',                                   'TEXT',                 0],
            ['DropVolume',                           'DropVolume',                              'TEXT',                 0],
            ['ProteinName',                          'ProteinName',                             'TEXT',                 1],
            ['BatchNumber',                          'BatchNumber',                             'TEXT',                 0],
            ['CompoundStockConcentration',           'CompoundStockConcentration',              'TEXT',                 0],
            ['CompoundConcentration',                'CompoundConcentration',                   'TEXT',                 1],
            ['SolventFraction',                      'SolventFraction',                         'TEXT',                 1],
            ['SoakTransferVol',                      'SoakTransferVol',                         'TEXT',                 0],
            ['SoakStatus',                           'SoakStatus',                              'TEXT',                 0],
            ['SoakTimestamp',                        'SoakTimestamp',                           'TEXT',                 0],
            ['CryoStockFraction',                    'CryoStockFraction',                       'TEXT',                 0],
            ['CryoFraction',                         'CryoFraction',                            'TEXT',                 0],
            ['CryoWell',                             'CryoWell',                                'TEXT',                 0],
            ['CryoTransferVolume',                   'CryoTransferVolume',                      'TEXT',                 0],
            ['CryoStatus',                           'CryoStatus',                              'TEXT',                 0],
            ['CryoTimestamp',                        'CryoTimestamp',                           'TEXT',                 0],
            ['SoakingTime',                          'SoakingTime',                             'TEXT',                 1],
            ['HarvestStatus',                        'HarvestStatus',                           'TEXT',                 0],
            ['CrystalName',                          'Sample ID',                               'TEXT',                 0],
            ['Puck',                                 'Puck',                                    'TEXT',                 1],
            ['PuckPosition',                         'PuckPosition',                            'TEXT',                 1],
            ['PinBarcode',                           'SoakDB\nBarcode',                         'TEXT',                 1],
            ['MountingResult',                       'MountingResult',                          'TEXT',                 0],
            ['MountingArrivalTime',                  'MountingArrivalTime',                     'TEXT',                 0],
            ['MountedTimestamp',                     'MountedTimestamp',                        'TEXT',                 0],
            ['MountingTime',                         'MountingTime',                            'TEXT',                 0],
            ['ispybStatus',                          'ispybStatus',                             'TEXT',                 0],
            ['DataCollectionVisit',                  'Visit',                                   'TEXT',                 1],

            # from XChemExplorer

            ['ProjectDirectory',                     'ProjectDirectory',                        'TEXT',                 0],

            ['CrystalTag',                           'Tag',                                     'TEXT',                 0],
            ['CrystalFormName',                      'Crystal Form\nName',                      'TEXT',                 0],
            ['CrystalFormSpaceGroup',                'Space\nGroup',                            'TEXT',                 0],
            ['CrystalFormPointGroup',                'Point\nGroup',                            'TEXT',                 0],
            ['CrystalFormA',                         'a',                                       'TEXT',                 0],
            ['CrystalFormB',                         'b',                                       'TEXT',                 0],
            ['CrystalFormC',                         'c',                                       'TEXT',                 0],
            ['CrystalFormAlpha',                     'alpha',                                   'TEXT',                 0],
            ['CrystalFormBeta',                      'beta',                                    'TEXT',                 0],
            ['CrystalFormGamma',                     'gamma',                                   'TEXT',                 0],
            ['CrystalFormVolume',                    'Crystal Form\nVolume',                    'TEXT',                 0],

            ['DataCollectionBeamline',               'Beamline',                                'TEXT',                 0],
            ['DataCollectionDate',                   'Data Collection\nDate',                   'TEXT',                 1],
            ['DataCollectionOutcome',                'DataCollection\nOutcome',                 'TEXT',                 1],
            ['DataCollectionRun',                    'Run',                                     'TEXT',                 0],
            ['DataCollectionSubdir',                 'SubDir',                                  'TEXT',                 0],
            ['DataCollectionComment',                'DataCollection\nComment',                 'TEXT',                 0],
            ['DataCollectionWavelength',             'Wavelength',                              'TEXT',                 0],
            ['DataCollectionPinBarcode',             'GDA\nBarcode',                            'TEXT',                 1],

            ['DataCollectionCrystalImage1', 'img1', 'TEXT', 1],
            ['DataCollectionCrystalImage2', 'img2', 'TEXT', 1],
            ['DataCollectionCrystalImage3', 'img3', 'TEXT', 1],
            ['DataCollectionCrystalImage4', 'img4', 'TEXT', 1],

            ['DataProcessingPathToImageFiles',       'Path to diffraction\nimage files',        'TEXT',                 1],
            ['DataProcessingProgram',                'Program',                                 'TEXT',                 1],
            ['DataProcessingSpaceGroup',             'DataProcessing\nSpaceGroup',              'TEXT',                 1],
            ['DataProcessingUnitCell',               'DataProcessing\nUnitCell',                'TEXT',                 0],
            ['DataProcessingAutoAssigned',           'auto-assigned',                           'TEXT',                 0],

            ['DataProcessingA',                      'DataProcessing\nA',                       'TEXT',                 0],
            ['DataProcessingB',                      'DataProcessing\nB',                       'TEXT',                 0],
            ['DataProcessingC',                      'DataProcessing\nC',                       'TEXT',                 0],
            ['DataProcessingAlpha',                  'DataProcessing\nAlpha',                   'TEXT',                 0],
            ['DataProcessingBeta',                   'DataProcessing\nBeta',                    'TEXT',                 0],
            ['DataProcessingGamma',                  'DataProcessing\nGamma',                   'TEXT',                 0],

        # True or False depending on whether user changed automatic selection
#            ['DataProcessingUserSelected',                  'DataProcessingUserSelected',                   'TEXT', 0],

            ['DataProcessingResolutionOverall',             'Resolution\nOverall',                          'TEXT',                 0],
            ['DataProcessingResolutionLow',                 'Resolution\nLow',                              'TEXT',                 0],
            ['DataProcessingResolutionLowInnerShell',       'Resolution\nLow (Inner Shell)',                'TEXT',                 0],
            ['DataProcessingResolutionHigh',                'Resolution\nHigh',                             'TEXT',                 1],
            ['DataProcessingResolutionHigh15sigma',         'Resolution\n[Mn<I/sig(I)> = 1.5]',             'TEXT',                 1],
            ['DataProcessingResolutionHigh20sigma',         'Resolution\n[Mn<I/sig(I)> = 2.0]',             'TEXT',                 1],
            ['DataProcessingResolutionHighOuterShell',      'Resolution\nHigh (Outer Shell)',               'TEXT',                 0],
            ['DataProcessingRmergeOverall',                 'Rmerge\nOverall',                              'TEXT',                 1],
            ['DataProcessingRmergeLow',                     'Rmerge\nLow',                                  'TEXT',                 1],
            ['DataProcessingRmergeHigh',                    'Rmerge\nHigh',                                 'TEXT',                 1],
            ['DataProcessingIsigOverall',                   'Mn<I/sig(I)>\nOverall',                        'TEXT',                 1],
            ['DataProcessingIsigLow',                       'Mn<I/sig(I)>\nLow',                            'TEXT',                 1],
            ['DataProcessingIsigHigh',                      'Mn<I/sig(I)>\nHigh',                           'TEXT',                 1],
            ['DataProcessingCompletenessOverall',           'Completeness\nOverall',                        'TEXT',                 1],
            ['DataProcessingCompletenessLow',               'Completeness\nLow',                            'TEXT',                 1],
            ['DataProcessingCompletenessHigh',              'Completeness\nHigh',                           'TEXT',                 1],
            ['DataProcessingMultiplicityOverall',           'Multiplicity\nOverall',                        'TEXT',                 1],
            ['DataProcessingMultiplicityLow',               'Multiplicity\nLow',                            'TEXT',                 1],
            ['DataProcessingMultiplicityHigh',              'Multiplicity\nHigh',                           'TEXT',                 1],
            ['DataProcessingCChalfOverall',                 'CC(1/2)\nOverall',                             'TEXT',                 1],
            ['DataProcessingCChalfLow',                     'CC(1/2)\nLow',                                 'TEXT',                 1],
            ['DataProcessingCChalfHigh',                    'CC(1/2)\nHigh',                                'TEXT',                 1],

            # the data source is a bit exploding with entries like the ones below,
            # but the many different filenames and folder structures of Diamond autoprocessing makes it necessary
            ['DataProcessingPathToLogfile',                 'DataProcessingPathToLogfile',                  'TEXT',                 1],
            ['DataProcessingPathToMTZfile',                 'DataProcessingPathToMTZfile',                  'TEXT',                 1],
            ['DataProcessingLOGfileName',                   'DataProcessingLOGfileName',                    'TEXT',                 0],
            ['DataProcessingMTZfileName',                   'DataProcessingMTZfileName',                    'TEXT',                 0],
            ['DataProcessingDirectoryOriginal',             'DataProcessingDirectoryOriginal',              'TEXT',                 0],

            ['DataProcessingUniqueReflectionsOverall',      'Unique Reflections\nOverall',                  'TEXT',                 1],
            ['DataProcessingLattice',                       'DataProcessing\nLattice',                      'TEXT',                 0],
            ['DataProcessingPointGroup',                    'DataProcessing\nPointGroup',                   'TEXT',                 0],
            ['DataProcessingUnitCellVolume',                'DataProcessing\nUnit Cell Volume',             'TEXT',                 0],
            ['DataProcessingAlert',                         'DataProcessing\nAlert',                        'TEXT',                 0],
            ['DataProcessingScore',                         'DataProcessing\nScore',                        'TEXT',                 1],

            ['DataProcessingStatus',                        'DataProcessing\nStatus',                       'TEXT',                 1],

            ['DataProcessingRcryst',                        'DataProcessing\nRcryst',                       'TEXT',                 0],
            ['DataProcessingRfree',                         'DataProcessing\nRfree',                        'TEXT',                 0],
            ['DataProcessingPathToDimplePDBfile',           'DataProcessingPathToDimplePDBfile',            'TEXT',                 0],
            ['DataProcessingPathToDimpleMTZfile',           'DataProcessingPathToDimpleMTZfile',            'TEXT',                 0],
            ['DataProcessingDimpleSuccessful',              'DataProcessingDimpleSuccessful',               'TEXT',                 0],

            ['DimpleResolutionHigh',                        'Dimple\nResolution High',                      'TEXT',                 1],
            ['DimpleRcryst',                                'Dimple\nRcryst',                               'TEXT',                 1],
            ['DimpleRfree',                                 'Dimple\nRfree',                                'TEXT',                 1],
            ['DimplePathToPDB',                             'Dimple\nPath to PDB file',                     'TEXT',                 1],
            ['DimplePathToMTZ',                             'Dimple\nPath to MTZ file',                     'TEXT',                 1],
            ['DimpleReferencePDB',                          'Dimple\nReference PDB',                        'TEXT',                 1],
            ['DimpleStatus',                                'Dimple\nStatus',                               'TEXT',                 1],

            ['DimplePANDDAwasRun',                          'PanDDA\nlaunched?',                            'TEXT',                 1],
            ['DimplePANDDAhit',                             'PanDDA\nhit?',                                 'TEXT',                 1],
            ['DimplePANDDAreject',                          'PanDDA\nreject?',                              'TEXT',                 1],
            ['DimplePANDDApath',                            'PanDDA\npath?',                                'TEXT',                 1],
            ['PANDDAStatus',                                'PanDDA\nStatus',                               'TEXT',                 1],

            ['DatePanDDAModelCreated',                      'DatePanDDAModelCreated',                       'TEXT',                 0],

            ['RefinementResolution',                        'Refinement\nResolution',                       'TEXT',                 1],
            ['RefinementResolutionTL',                      'RefinementResolutionTL',                       'TEXT',                 0],
            ['RefinementRcryst',                            'Refinement\nRcryst',                           'TEXT',                 1],
            ['RefinementRcrystTraficLight',                 'RefinementRcrystTraficLight',                  'TEXT',                 0],
            ['RefinementRfree',                             'Refinement\nRfree',                            'TEXT',                 1],
            ['RefinementRfreeTraficLight',                  'RefinementRfreeTraficLight',                   'TEXT',                 0],
            ['RefinementSpaceGroup',                        'Refinement\nSpace Group',                      'TEXT',                 1],
            ['RefinementLigandCC',                          'RefinementLigandCC',                           'TEXT',                 0],
            ['RefinementRmsdBonds',                         'RefinementRmsdBonds',                          'TEXT',                 1],
            ['RefinementRmsdBondsTL',                       'RefinementRmsdBondsTL',                        'TEXT',                 0],
            ['RefinementRmsdAngles',                        'RefinementRmsdAngles',                         'TEXT',                 1],
            ['RefinementRmsdAnglesTL',                      'RefinementRmsdAnglesTL',                       'TEXT',                 0],
            ['RefinementOutcome',                           'Refinement\nOutcome',                          'TEXT',                 1],
            ['RefinementMTZfree',                           'RefinementMTZfree',                            'TEXT',                 1],
            ['RefinementCIF',                               'RefinementCIF',                                'TEXT',                 1],
            ['RefinementCIFStatus',                         'Compound\nStatus',                             'TEXT',                 1],
            ['RefinementCIFprogram',                        'RefinementCIFprogram',                         'TEXT',                 1],
            ['RefinementPDB_latest',                        'RefinementPDB_latest',                         'TEXT',                 1],
            ['RefinementMTZ_latest',                        'RefinementMTZ_latest',                         'TEXT',                 1],
            ['RefinementMatrixWeight',                      'RefinementMatrixWeight',                       'TEXT',                 0],
            ['RefinementComment',                           'RefinementComment',                            'TEXT',                 0],
            ['RefinementPathToRefinementFolder',            'RefinementPathToRefinementFolder',             'TEXT',                 0],
            ['RefinementLigandConfidence',                  'Ligand\nConfidence',                           'TEXT',                 0],
            ['RefinementLigandBoundConformation',           'RefinementLigandBoundConformation',            'TEXT',                 0],
            ['RefinementBoundConformation',                 'RefinementBoundConformation',                  'TEXT',                 0],
            ['RefinementMolProbityScore',                   'MolProbity Score',                             'TEXT',                 1],
            ['RefinementMolProbityScoreTL',                 'RefinementMolProbityScoreTL',                  'TEXT',                 0],
            ['RefinementRamachandranOutliers',              'Ramachandran\nOutliers',                       'TEXT',                 1],
            ['RefinementRamachandranOutliersTL',            'RefinementRamachandranOutliersTL',             'TEXT',                 0],
            ['RefinementRamachandranFavored',               'Ramachandran\nFavored',                        'TEXT',                 1],
            ['RefinementRamachandranFavoredTL',             'RefinementRamachandranFavoredTL',              'TEXT',                 0],
            ['RefinementProgram',                           'RefinementProgram',                            'TEXT',                 1],
            ['RefinementStatus',                            'Refinement\nStatus',                           'TEXT',                 1],

            ['Deposition_PDB_ID',                           'Deposition_PDB_ID',                            'TEXT',                 1],
            ['Deposition_PDB_file',                         'Deposition_PDB_file',                          'TEXT',                 0],
            ['Deposition_Date',                             'Deposition_Date',                              'TEXT',                 1],
            ['Deposition_mmCIF_model_file',                 'Deposition_mmCIF_model_file',                  'TEXT',                 0],
            ['Deposition_mmCIF_SF_file',                    'Deposition_mmCIF_SF_file',                     'TEXT',                 0],

            ['Label',                                       'Label',                                        'TEXT',                 0],
            ['table_one',                                   'table_one',                                    'TEXT',                 0],

            ['AssayIC50',                                   'AssayIC50',                                    'TEXT',                 0],
            ['LastUpdated',                                 'LastUpdated',                                  'TEXT',                 0],
            ['LastUpdated_by',                              'LastUpdated_by',                               'TEXT',                 0]
        ]

        self.pandda_table_columns = [
            ['ID',                                          'ID',                                       'INTEGER PRIMARY KEY'],
            ['CrystalName',                                 'Sample ID',                                'TEXT'],
            ['PANDDApath',                                  'PANDDApath',                               'TEXT'],
            ['PANDDA_site_index',                           'PANDDA_site_index',                        'TEXT'],
            ['PANDDA_site_name',                            'PANDDA_site_name',                         'TEXT'],
            ['PANDDA_site_comment',                         'PANDDA_site_comment',                      'TEXT'],
            ['PANDDA_site_event_index',                     'PANDDA_site_event_index',                  'TEXT'],
            ['PANDDA_site_event_comment',                   'PANDDA_site_event_comment',                'TEXT'],
            ['PANDDA_site_confidence',                      'PANDDA_site_confidence',                   'TEXT'],
            ['PANDDA_site_InspectConfidence',               'PANDDA_site_InspectConfidence',            'TEXT'],
            ['PANDDA_site_ligand_placed',                   'PANDDA_site_ligand_placed',                'TEXT'],
            ['PANDDA_site_viewed',                          'PANDDA_site_viewed',                       'TEXT'],
            ['PANDDA_site_interesting',                     'PANDDA_site_interesting',                  'TEXT'],
            ['PANDDA_site_z_peak',                          'PANDDA_site_z_peak',                       'TEXT'],
            ['PANDDA_site_x',                               'PANDDA_site_x',                            'TEXT'],
            ['PANDDA_site_y',                               'PANDDA_site_y',                            'TEXT'],
            ['PANDDA_site_z',                               'PANDDA_site_z',                            'TEXT'],
            ['PANDDA_site_ligand_id',                       'PANDDA_site_ligand_id',                    'TEXT'],

            ['PANDDA_site_ligand_resname',                  'PANDDA_site_ligand_resname',               'TEXT'],
            ['PANDDA_site_ligand_chain',                    'PANDDA_site_ligand_chain',                 'TEXT'],
            ['PANDDA_site_ligand_sequence_number',          'PANDDA_site_ligand_sequence_number',       'TEXT'],
            ['PANDDA_site_ligand_altLoc',                   'PANDDA_site_ligand_altLoc',                'TEXT'],

            ['PANDDA_site_event_map',                       'PANDDA_site_event_map',                    'TEXT'],
            ['PANDDA_site_event_map_mtz',                   'PANDDA_site_event_map_mtz',                'TEXT'],
            ['PANDDA_site_initial_model',                   'PANDDA_site_initial_model',                'TEXT'],
            ['PANDDA_site_initial_mtz',                     'PANDDA_site_initial_mtz',                  'TEXT'],
            ['PANDDA_site_spider_plot',                     'PANDDA_site_spider_plot',                  'TEXT'],
            ['PANDDA_site_occupancy',                       'PANDDA_site_occupancy',                    'TEXT'],
            ['PANDDA_site_B_average',                       'PANDDA_site_B_average',                    'TEXT'],
            ['PANDDA_site_B_ratio_residue_surroundings',    'PANDDA_site_B_ratio_residue_surroundings', 'TEXT'],
            ['PANDDA_site_RSCC',                            'PANDDA_site_RSCC',                         'TEXT'],
            ['PANDDA_site_RSR',                             'PANDDA_site_RSR',                          'TEXT'],
            ['PANDDA_site_RSZD',                            'PANDDA_site_RSZD',                         'TEXT'],
            ['PANDDA_site_rmsd',                            'PANDDA_site_rmsd',                         'TEXT'],
            ['RefinementOutcome',                           'RefinementOutcome',                        'TEXT'],
            ['ApoStructures',                               'ApoStructures',                            'TEXT'],
            ['LastUpdated',                                 'LastUpdated',                              'TEXT'],
            ['LastUpdated_by',                              'LastUpdated_by',                           'TEXT']
            ]

        self.deposition_table_columns = [
            ['ID',                                          'ID',                                       'INTEGER PRIMARY KEY'],

            ['CrystalName',                                 'Sample ID',                                'TEXT'],
            ['StructureType',                               'StructureType',                            'TEXT'],    # apo/model

            ['PDB_file',                                    'PDB_file',                                 'TEXT'],
            ['MTZ_file',                                    'MTZ_file',                                 'TEXT'],

            ['mmCIF_model_file',                            'mmCIF_model_file',                         'TEXT'],
            ['mmCIF_SF_file',                               'mmCIF_SF_file',                            'TEXT'],
            ['label',                                       'label',                                    'TEXT'],    # for index.txt
            ['description',                                 'description',                              'TEXT'],    # for index.txt

            ['DimplePANDDApath',                            'DimplePANDDApath',                         'TEXT'],

            ['contact_author_PI_salutation',                'contact_author_PI_salutation',             'TEXT'],
            ['contact_author_PI_first_name',                'contact_author_PI_first_name',             'TEXT'],
            ['contact_author_PI_last_name',                 'contact_author_PI_last_name',              'TEXT'],
            ['contact_author_PI_middle_name',               'contact_author_PI_middle_name',            'TEXT'],
            ['contact_author_PI_role',                      'contact_author_PI_role',                   'TEXT'],
            ['contact_author_PI_organization_type',         'contact_author_PI_organization_type',      'TEXT'],
            ['contact_author_PI_organization_name',         'contact_author_PI_organization_name',      'TEXT'],
            ['contact_author_PI_email',                     'contact_author_PI_email',                  'TEXT'],
            ['contact_author_PI_address',                   'contact_author_PI_address',                'TEXT'],
            ['contact_author_PI_city',                      'contact_author_PI_city',                   'TEXT'],
            ['contact_author_PI_State_or_Province',         'contact_author_PI_State_or_Province',      'TEXT'],
            ['contact_author_PI_Zip_Code',                  'contact_author_PI_Zip_Code',               'TEXT'],
            ['contact_author_PI_Country',                   'contact_author_PI_Country',                'TEXT'],
            ['contact_author_PI_fax_number',                'contact_author_PI_fax_number',             'TEXT'],
            ['contact_author_PI_phone_number',              'contact_author_PI_phone_number',           'TEXT'],
            ['contact_author_PI_ORCID',                     'contact_author_PI_ORCID',                  'TEXT'],

            ['contact_author_salutation',                   'contact_author_salutation',                'TEXT'],
            ['contact_author_first_name',                   'contact_author_first_name',                'TEXT'],
            ['contact_author_last_name',                    'contact_author_last_name',                 'TEXT'],
            ['contact_author_middle_name',                  'contact_author_middle_name',               'TEXT'],
            ['contact_author_role',                         'contact_author_role',                      'TEXT'],
            ['contact_author_organization_type',            'contact_author_organization_type',         'TEXT'],
            ['contact_author_organization_name',            'contact_author_organization_name',         'TEXT'],
            ['contact_author_email',                        'contact_author_email',                     'TEXT'],
            ['contact_author_address',                      'contact_author_address',                   'TEXT'],
            ['contact_author_city',                         'contact_author_city',                      'TEXT'],
            ['contact_author_State_or_Province',            'contact_author_State_or_Province',         'TEXT'],
            ['contact_author_Zip_Code',                     'contact_author_Zip_Code',                  'TEXT'],
            ['contact_author_Country',                      'contact_author_Country',                   'TEXT'],
            ['contact_author_fax_number',                   'contact_author_fax_number',                'TEXT'],
            ['contact_author_phone_number',                 'contact_author_phone_number',              'TEXT'],
            ['contact_author_ORCID',                        'contact_author_ORCID',                     'TEXT'],

            ['Release_status_for_coordinates',              'Release_status_for_coordinates',           'TEXT'],
            ['Release_status_for_structure_factor',         'Release_status_for_structure_factor',      'TEXT'],
            ['Release_status_for_sequence',                 'Release_status_for_sequence',              'TEXT'],

            ['group_deposition_title',                      'group_deposition_title',                   'TEXT'],
            ['group_description',                           'group_description',                        'TEXT'],
            ['structure_title',                             'structure_title',                          'TEXT'],
            ['structure_details',                           'structure_details',                        'TEXT'],
            ['group_deposition_title_apo',                  'group_deposition_title_apo',               'TEXT'],
            ['structure_title_apo',                         'structure_title_apo',                      'TEXT'],

            ['structure_author_name',                       'structure_author_name',                    'TEXT'],
            ['primary_citation_author_name',                'primary_citation_author_name',             'TEXT'],

            ['primary_citation_id',                         'primary_citation_id',                      'TEXT'],
            ['primary_citation_journal_abbrev',             'primary_citation_journal_abbrev',          'TEXT'],
            ['primary_citation_title',                      'primary_citation_title',                   'TEXT'],
            ['primary_citation_year',                       'primary_citation_year',                    'TEXT'],
            ['primary_citation_journal_volume',             'primary_citation_journal_volume',          'TEXT'],
            ['primary_citation_page_first',                 'primary_citation_page_first',              'TEXT'],
            ['primary_citation_page_last',                  'primary_citation_page_last',               'TEXT'],

            ['molecule_name',                               'molecule_name',                            'TEXT'],
            ['Fragment_name',                               'Fragment_name',                            'TEXT'],
            ['Specific_mutation',                           'Specific_mutation',                        'TEXT'],
            ['Enzyme_Comission_number',                     'Enzyme_Comission_number',                  'TEXT'],
            ['Source_organism_scientific_name',             'Source_organism_scientific_name',          'TEXT'],
            ['Source_organism_gene',                        'Source_organism_gene',                     'TEXT'],
            ['Source_organism_strain',                      'Source_organism_strain',                   'TEXT'],
            ['Expression_system_scientific_name',           'Expression_system_scientific_name',        'TEXT'],
            ['Expression_system_strain',                    'Expression_system_strain',                 'TEXT'],
            ['Expression_system_vector_type',               'Expression_system_vector_type',            'TEXT'],
            ['Expression_system_plasmid_name',              'Expression_system_plasmid_name',           'TEXT'],
            ['Manipulated_source_details',                  'Manipulated_source_details',               'TEXT'],
            ['fragment_name_one_specific_mutation',         'fragment_name_one_specific_mutation',      'TEXT'],
            ['molecule_chain_one',                          'molecule_chain_one',                       'TEXT'],

            ['molecule_name_two',                           'molecule_name_two',                        'TEXT'],
            ['Fragment_name_two',                           'Fragment_name_two',                        'TEXT'],
            ['Specific_mutation_two',                       'Specific_mutation_two',                    'TEXT'],
            ['Enzyme_Comission_number_two',                 'Enzyme_Comission_number_two',              'TEXT'],
            ['Source_organism_scientific_name_two',         'Source_organism_scientific_name_two',      'TEXT'],
            ['Source_organism_gene_two',                    'Source_organism_gene_two',                 'TEXT'],
            ['Source_organism_strain_two',                  'Source_organism_strain_two',               'TEXT'],
            ['Expression_system_scientific_name_two',       'Expression_system_scientific_name_two',    'TEXT'],
            ['Expression_system_strain_two',                'Expression_system_strain_two',             'TEXT'],
            ['Expression_system_vector_type_two',           'Expression_system_vector_type_two',        'TEXT'],
            ['Expression_system_plasmid_name_two',          'Expression_system_plasmid_name_two',       'TEXT'],
            ['Manipulated_source_details_two',              'Manipulated_source_details_two',           'TEXT'],
            ['fragment_name_two_specific_mutation',         'fragment_name_one_specific_mutation_two',  'TEXT'],
            ['molecule_chain_two',                          'molecule_chain_two',                       'TEXT'],

            ['structure_keywords',                          'structure_keywords',                       'TEXT'],
            ['biological_assembly_chain_number',            'biological_assembly_chain_number',         'TEXT'],


            ['crystallization_id',                          'crystallization_id',                       'TEXT'],
            ['crystallization_method',                      'crystallization_method',                   'TEXT'],
            ['crystallization_pH',                          'crystallization_pH',                       'TEXT'],
            ['crystallization_temperature',                 'crystallization_temperature',              'TEXT'],
            ['crystallization_details',                     'crystallization_details',                  'TEXT'],

            ['radiation_source',                            'radiation_source',                         'TEXT'],
            ['radiation_source_type',                       'radiation_source_type',                    'TEXT'],
            ['radiation_wavelengths',                       'radiation_wavelengths',                    'TEXT'],
            ['radiation_detector',                          'radiation_detector',                       'TEXT'],
            ['radiation_detector_type',                     'radiation_detector_type',                  'TEXT'],
            ['data_collection_date',                        'data_collection_date',                     'TEXT'],
            ['data_collection_temperature',                 'data_collection_temperature',              'TEXT'],
            ['data_collection_protocol',                    'data_collection_protocol',                 'TEXT'],

            ['SG_project_name',                             'SG_project_name',                          'TEXT'],
            ['full_name_of_SG_center',                      'full_name_of_SG_center',                   'TEXT'],

            ['molecule_one_letter_sequence',                'molecule_one_letter_sequence',             'TEXT'],
            ['molecule_one_letter_sequence_uniprot_id',     'molecule_one_letter_sequence_uniprot_id',  'TEXT'],
            ['molecule_two_letter_sequence',                'molecule_two_letter_sequence',             'TEXT'],
            ['molecule_two_letter_sequence_uniprot_id',     'molecule_two_letter_sequence_uniprot_id',  'TEXT'],

            ['CrystalName_of_pandda_input',                 'CrystalName_of_pandda_input',              'TEXT'],

            ['pdbx_starting_model',                         'pdbx_starting_model',                      'TEXT'],
            ['data_integration_software',                   'data_integration_software',                'TEXT'],
            ['phasing_software',                            'phasing_software',                         'TEXT'],

            ['LastUpdated',                                 'LastUpdated',                              'TEXT'],
            ['LastUpdated_by',                              'LastUpdated_by',                           'TEXT']
            ]



        self.data_collection_columns = [
            ['ID',                                          'ID',                                       'INTEGER PRIMARY KEY'],
            ['CrystalName',                          'Sample ID',                               'TEXT',                 0],
            ['ProteinName',                         'ProteinName',                              'TEXT', 1],

            ['DataCollectionVisit',                    'Visit',                                     'TEXT',                 0],
            ['DataCollectionRun',                    'Run',                                     'TEXT',                 0],
            ['DataCollectionSubdir',                    'SubDir',                                     'TEXT',                 0],
            ['DataCollectionBeamline',               'Beamline',                                'TEXT',                 0],
            ['DataCollectionOutcome',                'DataCollection\nOutcome',                 'TEXT',                 1],
            ['DataCollectionDate',                   'Data Collection\nDate',                   'TEXT',                 1],
            ['DataCollectionWavelength',             'Wavelength',                              'TEXT',                 0],
            ['DataCollectionPinBarcode',             'GDA\nBarcode',                            'TEXT',                 1],

            ['DataCollectionCrystalImage1',             'img1',                            'TEXT',                 1],
            ['DataCollectionCrystalImage2',             'img2',                            'TEXT',                 1],
            ['DataCollectionCrystalImage3',             'img3',                            'TEXT',                 1],
            ['DataCollectionCrystalImage4',             'img4',                            'TEXT',                 1],

            ['DataProcessingPathToImageFiles',       'Path to diffraction\nimage files',        'TEXT',                 1],
            ['DataProcessingProgram',                'Program',                                 'TEXT',                 1],
            ['DataProcessingSpaceGroup',             'DataProcessing\nSpaceGroup',              'TEXT',                 1],
            ['DataProcessingUnitCell',               'DataProcessing\nUnitCell',                'TEXT',                 0],
            ['DataProcessingAutoAssigned',           'auto-assigned',                           'TEXT',                 0],
            ['DataProcessingA',                      'DataProcessing\nA',                       'TEXT',                 0],
            ['DataProcessingB',                      'DataProcessing\nB',                       'TEXT',                 0],
            ['DataProcessingC',                      'DataProcessing\nC',                       'TEXT',                 0],
            ['DataProcessingAlpha',                  'DataProcessing\nAlpha',                   'TEXT',                 0],
            ['DataProcessingBeta',                   'DataProcessing\nBeta',                    'TEXT',                 0],
            ['DataProcessingGamma',                  'DataProcessing\nGamma',                   'TEXT',                 0],
            ['DataProcessingResolutionOverall',             'Resolution\nOverall',                          'TEXT',                 0],
            ['DataProcessingResolutionLow',                 'Resolution\nLow',                              'TEXT',                 0],
            ['DataProcessingResolutionLowInnerShell',       'Resolution\nLow (Inner Shell)',                'TEXT',                 0],
            ['DataProcessingResolutionHigh',                'Resolution\nHigh',                             'TEXT',                 1],
            ['DataProcessingResolutionHigh15sigma',         'Resolution\n[Mn<I/sig(I)> = 1.5]',             'TEXT',                 1],
            ['DataProcessingResolutionHigh20sigma',         'Resolution\n[Mn<I/sig(I)> = 2.0]',             'TEXT',                 1],
            ['DataProcessingResolutionHighOuterShell',      'Resolution\nHigh (Outer Shell)',               'TEXT',                 0],
            ['DataProcessingRmergeOverall',                 'Rmerge\nOverall',                              'TEXT',                 1],
            ['DataProcessingRmergeLow',                     'Rmerge\nLow',                                  'TEXT',                 1],
            ['DataProcessingRmergeHigh',                    'Rmerge\nHigh',                                 'TEXT',                 1],
            ['DataProcessingIsigOverall',                   'Mn<I/sig(I)>\nOverall',                        'TEXT',                 1],
            ['DataProcessingIsigLow',                       'Mn<I/sig(I)>\nLow',                            'TEXT',                 1],
            ['DataProcessingIsigHigh',                      'Mn<I/sig(I)>\nHigh',                           'TEXT',                 1],
            ['DataProcessingCompletenessOverall',           'Completeness\nOverall',                        'TEXT',                 1],
            ['DataProcessingCompletenessLow',               'Completeness\nLow',                            'TEXT',                 1],
            ['DataProcessingCompletenessHigh',              'Completeness\nHigh',                           'TEXT',                 1],
            ['DataProcessingMultiplicityOverall',           'Multiplicity\nOverall',                        'TEXT',                 1],
            ['DataProcessingMultiplicityLow',               'Multiplicity\nLow',                            'TEXT',                 1],
            ['DataProcessingMultiplicityHigh',              'Multiplicity\nHigh',                           'TEXT',                 1],
            ['DataProcessingCChalfOverall',                 'CC(1/2)\nOverall',                             'TEXT',                 1],
            ['DataProcessingCChalfLow',                     'CC(1/2)\nLow',                                 'TEXT',                 1],
            ['DataProcessingCChalfHigh',                    'CC(1/2)\nHigh',                                'TEXT',                 1],
            ['DataProcessingPathToLogfile',                 'DataProcessingPathToLogfile',                  'TEXT',                 1],
            ['DataProcessingPathToMTZfile',                 'DataProcessingPathToMTZfile',                  'TEXT',                 1],
            ['DataProcessingLOGfileName',                   'DataProcessingLOGfileName',                    'TEXT',                 0],
            ['DataProcessingMTZfileName',                   'DataProcessingMTZfileName',                    'TEXT',                 0],
            ['DataProcessingDirectoryOriginal',             'DataProcessingDirectoryOriginal',              'TEXT',                 0],
            ['DataProcessingUniqueReflectionsOverall',      'Unique Reflections\nOverall',                  'TEXT',                 1],
            ['DataProcessingLattice',                       'DataProcessing\nLattice',                      'TEXT',                 0],
            ['DataProcessingPointGroup',                    'DataProcessing\nPointGroup',                   'TEXT',                 0],
            ['DataProcessingUnitCellVolume',                'DataProcessing\nUnit Cell Volume',             'TEXT',                 0],
            ['DataProcessingAlert',                         'DataProcessing\nAlert',                        'TEXT',                 0],
            ['DataProcessingScore',                         'DataProcessing\nScore',                        'TEXT',                 1],
            ['DataProcessingStatus',                        'DataProcessing\nStatus',                       'TEXT',                 1],
            ['LastUpdated',                                 'LastUpdated',                                  'TEXT',                 0],
            ['LastUpdated_by',                              'LastUpdated_by',                               'TEXT',                 0]
        ]


        self.zenodo_table_columns = [
            ['ID',                          'ID',                           'INTEGER PRIMARY KEY'   ],
            ['DimplePANDDApath',            'DimplePANDDApath',             'TEXT'                  ],
            ['ZenodoTitle',                 'ZenodoTitle',                  'TEXT'                  ],
            ['ZenodoHTTPS',                 'ZenodoHTTPS',                  'TEXT'                  ],
            ['ZenodoDOI',                   'ZenodoDOI',                    'TEXT'                  ],
            ['LastUpdated',                 'LastUpdated',                  'TEXT'                  ],
            ['LastUpdated_by',              'LastUpdated_by',               'TEXT'                  ]
        ]

        self.label_table_columns = [
            ['ID',                          'ID',                           'INTEGER PRIMARY KEY'   ],
            ['Label',                       'Label',                        'TEXT'                  ],
            ['Description',                 'Description',                  'TEXT'                  ],
        ]


    def columns_not_to_display(self):
        do_not_display = []
        for column in self.column_list:
            if column[3]==0:
                do_not_display.append(column[1])
        return do_not_display
        
    def get_empty_db_dict(self):
        db_dict={}
        for column in self.column_list:
            if column[0] != 'ID':
                db_dict[column[0]]=''
        return db_dict

    def create_missing_columns(self):
        existing_columns=[]
        connect=sqlite3.connect(self.data_source_file)
        connect.row_factory = sqlite3.Row
        cursor = connect.cursor()

        tableDict = {   'mainTable':        self.column_list,
                        'panddaTable':      self.pandda_table_columns,
                        'depositTable':     self.deposition_table_columns,
                        'collectionTable':  self.data_collection_columns,
                        'zenodoTable':      self.zenodo_table_columns,
                        'labelTable':       self.label_table_columns    }

        for table in tableDict:
            cursor.execute("create table if not exists "+table+" (ID INTEGER);")
            existing_columns = []
            cursor.execute("SELECT * FROM "+table)
            for column in cursor.description:
                existing_columns.append(column[0])
            for column in tableDict[table]:
                if column[0] not in existing_columns:
                    cursor.execute("alter table " + table + " add column '" + column[0] + "' '" + column[2] + "'")
                    connect.commit()
            if table == 'labelTable':
                cursor.execute('select ID from labelTable')
                id = cursor.fetchall()
                if id == []:
                    for idx in range(5):
                        cursor.execute("insert into labelTable (ID) Values (%s)" %str(idx+1))

#        existing_columns=[]
#        connect=sqlite3.connect(self.data_source_file)
#        connect.row_factory = sqlite3.Row
#        cursor = connect.cursor()
#        cursor.execute("SELECT * FROM mainTable")
#        for column in cursor.description:
#            existing_columns.append(column[0])
#        for column in self.column_list:
#            if column[0] not in existing_columns:
#                cursor.execute("alter table mainTable add column '"+column[0]+"' '"+column[2]+"'")
#                connect.commit()
#        # create PANDDA table if not exists
#        cursor.execute("create table if not exists panddaTable (ID INTEGER);")
#        existing_columns=[]
#        cursor.execute("SELECT * FROM panddaTable")
#        for column in cursor.description:
#            existing_columns.append(column[0])
#        for column in self.pandda_table_columns:
#            if column[0] not in existing_columns:
#                cursor.execute("alter table panddaTable add column '"+column[0]+"' '"+column[2]+"'")
#                connect.commit()
#        # create DEPOSIT table if not exists
#        cursor.execute("create table if not exists depositTable (ID INTEGER);")
#        existing_columns=[]
#        cursor.execute("SELECT * FROM depositTable")
#        for column in cursor.description:
#            existing_columns.append(column[0])
#        for column in self.deposition_table_columns:
#            if column[0] not in existing_columns:
#                cursor.execute("alter table depositTable add column '"+column[0]+"' '"+column[2]+"'")
#                connect.commit()


    def return_column_list(self):
        return self.column_list


    def create_empty_data_source_file(self):

        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        with connect:
            cursor = connect.cursor()
            cursor.execute("CREATE TABLE mainTable("+self.column_list[0][0]+' '+self.column_list[0][2]+")")
            for i in range(1,len(self.column_list)):
                cursor.execute("alter table mainTable add column '"+self.column_list[i][0]+"' '"+self.column_list[i][2]+"'")
            connect.commit()
        # Don't need to create panddaTable at this point, because table will be created by create_missing_columns
        # which is called the first time a data source in specified in XCE
#        with connect:
#            cursor = connect.cursor()
#            cursor.execute("CREATE TABLE panddaTable("+self.pandda_table_columns[0][0]+' '+self.pandda_table_columns[0][2]+")")
#            for i in range(1,len(self.pandda_table_columns)):
#                cursor.execute("alter table mainTable add column '"+self.pandda_table_columns[i][0]+"' '"+self.pandda_table_columns[i][2]+"'")
#            connect.commit()

    def get_all_samples_in_data_source_as_list(self):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("SELECT CrystalName FROM mainTable")
        existing_samples_in_db=[]
        samples = cursor.fetchall()
        for sample in samples:
            existing_samples_in_db.append(str(sample[0]))
        return existing_samples_in_db

    def execute_statement(self,cmd):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute(cmd)
        output=cursor.fetchall()
        connect.commit()
        return output

    def get_db_dict_for_sample(self,sampleID):
        db_dict={}
        header=[]
        data=[]
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("select * from mainTable where CrystalName='{0!s}';".format(sampleID))
        for column in cursor.description:
            header.append(column[0])
        data = cursor.fetchall()
        for n,item in enumerate(data[0]):
            db_dict[header[n]]=str(item)
        return db_dict

    def get_deposit_dict_for_sample(self,sampleID):
        db_dict={}
        header=[]
        data=[]
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
#        if sampleID == 'ground-state':      # just select first row in depositTable
#            cursor.execute("SELECT * FROM depositTable ORDER BY ROWID ASC LIMIT 1;")
#        else:
#            cursor.execute("select * from depositTable where CrystalName='{0!s}';".format(sampleID))
        cursor.execute("select * from depositTable where CrystalName='{0!s}';".format(sampleID))
        for column in cursor.description:
            header.append(column[0])
        data = cursor.fetchall()
        try:
            for n,item in enumerate(data[0]):
                db_dict[header[n]]=str(item)
        except IndexError:
            pass
        return db_dict

    def get_zenodo_dict_for_pandda_analysis(self,panddaPath):
        db_dict={}
        header=[]
        data=[]
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        try:
            cursor.execute("select * from zenodoTable where DimplePANDDApath='{0!s}';".format(panddaPath))
            for column in cursor.description:
                header.append(column[0])
            data = cursor.fetchall()
            try:
                for n,item in enumerate(data[0]):
                    db_dict[header[n]]=str(item)
            except IndexError:
                pass
        except sqlite3.OperationalError:
            db_dict = {}
        return db_dict



    def get_db_pandda_dict_for_sample_and_site(self,sampleID,site_index):
        db_dict={}
        header=[]
        data=[]
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("select * from panddaTable where CrystalName='{0!s}' and PANDDA_site_index='{1!s}';".format(sampleID, site_index))
        for column in cursor.description:
            header.append(column[0])
        data = cursor.fetchall()
        for n,item in enumerate(data[0]):
            db_dict[header[n]]=str(item)
        return db_dict

    def get_db_pandda_dict_for_sample_and_site_and_event(self,sampleID,site_index,event_index):
        db_dict={}
        header=[]
        data=[]
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("select * from panddaTable where CrystalName='{0!s}' and PANDDA_site_index='{1!s}' and PANDDA_site_event_index='{2!s}' and PANDDA_site_ligand_placed='True';".format(sampleID, site_index, event_index))
        for column in cursor.description:
            header.append(column[0])
        data = cursor.fetchall()
        for n,item in enumerate(data[0]):
            db_dict[header[n]]=str(item)
        return db_dict

    def check_if_sample_exists_in_data_source(self,sampleID):
        sample_exists=False
        existing_samples_in_db=self.get_all_samples_in_data_source_as_list()
        if sampleID in existing_samples_in_db:
            sample_exists=True
        return sample_exists

    def import_csv_file(self,csv_file):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        available_columns=[]
        cursor.execute("SELECT * FROM mainTable")
        for column in cursor.description:           # only update existing columns in data source
            available_columns.append(column[0])
        with open(csv_file,'rb') as csv_import: # `with` statement available in 2.5+
            # csv.DictReader uses first line in file for column headings by default
            csv_dict = csv.DictReader(csv_import) # comma is default delimiter
            for line in csv_dict:
                sampleID=line['CrystalName']
                if str(sampleID).replace(' ','')=='':
                    continue
                if self.check_if_sample_exists_in_data_source(sampleID):
                    update_string=''
                    for key,value in line.iteritems():
                        if key=='ID' or key=='CrystalName':
                            continue
                        if key not in available_columns:
                            continue
                        # this is how I had it originally, so it would ignore empty fields
                        # the problem is that if the user wants to set all values to Null,
                        # if will ignore it and leave the inital value in the datasource
#                        if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
#                            update_string+=str(key)+'='+"'"+str(value)+"'"+','
#                            print "UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"';"
#                            cursor.execute("UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"';")
                        # now try this instead; not sure what will break now...
                        update_string+=str(key)+'='+"'"+str(value)+"'"+','
#                        print "UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"';"
                        cursor.execute("UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"';")
                else:
                    column_string=''
                    value_string=''
                    for key,value in line.iteritems():
                        if key=='ID':
                            continue
                        if key not in available_columns:
                            continue
                        if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
                            value_string+="'"+value+"'"+','
                            column_string+=key+','
#                    print sampleID
#                    print          "INSERT INTO mainTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");"
                    cursor.execute("INSERT INTO mainTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");")

        connect.commit()

    def update_data_source(self,sampleID,data_dict):
        print 'here'
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()

        # need to do this since some older sqlite files contain a columnn called
        # DataProcessingResolutionHigh1.5sigma
        # and this does not go down well with the SQLite statement below
        removeKey=''
        for key in data_dict:
            if '.5' in key:
                removeKey=key
                break
        if removeKey != '':
            del data_dict[removeKey]

        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        update_string=''
        for key in data_dict:
            value=data_dict[key]
            if key=='ID' or key=='CrystalName':
                continue
            if not str(value).replace(' ','')=='':  # ignore empty fields
                update_string+=str(key)+'='+"'"+str(value)+"'"+','
            else:
                update_string+=str(key)+' = null,'
        if update_string != '':
#            print "UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"'"
            cursor.execute("UPDATE mainTable SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"'")
            connect.commit()

    def update_panddaTable(self,sampleID,site_index,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        update_string=''
        for key in data_dict:
            value=data_dict[key]
            if key=='ID' or key=='CrystalName' or key=='PANDDA_site_index':
                continue
            if not str(value).replace(' ','')=='':  # ignore empty fields
                update_string+=str(key)+'='+"'"+str(value)+"'"+','
            else:
                update_string+=str(key)+' = null,'
        if update_string != '':
            cursor.execute("UPDATE panddaTable SET "+update_string[:-1]+" WHERE CrystalName='{0!s}' and PANDDA_site_index='{1!s}'".format(sampleID, site_index))
            connect.commit()

    def update_site_event_panddaTable(self,sampleID,site_index,event_index,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        update_string=''
        for key in data_dict:
            value=data_dict[key]
            if key=='ID' or key=='CrystalName' or key=='PANDDA_site_index' or key=='PANDDA_site_event_index':
                continue
            if not str(value).replace(' ','')=='':  # ignore empty fields
                update_string+=str(key)+'='+"'"+str(value)+"'"+','
            else:
                update_string+=str(key)+' = null,'
        if update_string != '':
            cursor.execute("UPDATE panddaTable SET "+update_string[:-1]+" WHERE CrystalName='{0!s}' and PANDDA_site_index='{1!s}' and PANDDA_site_event_index='{2!s}'".format(sampleID, site_index, event_index))
            connect.commit()

    def update_depositTable(self,sampleID,structure_type,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        update_string=''
        for key in data_dict:
            value=data_dict[key]
            if key=='ID' or key=='CrystalName' or key=='StructureType':
                continue
            if not str(value).replace(' ','')=='':  # ignore empty fields
#                update_string+=str(key)+'='+"'"+str(value)+"'"+','
#                update_string+=str(key)+'='+'"'+str(value)+'"'+','
                update_string+=str(key)+'='+'"'+str(value).replace('\r','').replace('\n','')+'"'+','
            else:
                update_string+=str(key)+' = null,'
        if update_string != '':
#            print '-->',"UPDATE depositTable SET "+update_string[:-1]+" WHERE CrystalName='{0!s}' and StructureType='{1!s}'".format(sampleID, structure_type)
#            cursor.execute("UPDATE depositTable SET "+update_string[:-1]+" WHERE CrystalName='{0!s}' and StructureType='{1!s}'".format(sampleID, structure_type))
            cursor.execute('UPDATE depositTable SET '+update_string[:-1]+' WHERE CrystalName="{0!s}" and StructureType="{1!s}"'.format(sampleID, structure_type))
            connect.commit()


    def update_specified_table(self,sampleID,data_dict,table):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        update_string=''
        for key in data_dict:
            value=data_dict[key]
            if key=='ID' or key=='CrystalName':
                continue
            if not str(value).replace(' ','')=='':  # ignore empty fields
                update_string+=str(key)+'='+"'"+str(value)+"'"+','
        if update_string != '':
            cursor.execute("UPDATE "+table+" SET "+update_string[:-1]+" WHERE CrystalName="+"'"+sampleID+"'")
            connect.commit()

    def update_insert_data_source(self,sampleID,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute('Select CrystalName FROM mainTable')
        available_columns=[]
        cursor.execute("SELECT * FROM mainTable")
        for column in cursor.description:           # only update existing columns in data source
            available_columns.append(column[0])
        if self.check_if_sample_exists_in_data_source(sampleID):
            for key in data_dict:
                value=data_dict[key]
                if key=='ID' or key=='CrystalName':
                    continue
                if not str(value).replace(' ','')=='':  # ignore empty fields
                    update_string=str(key)+'='+"'"+str(value)+"'"
                    cursor.execute("UPDATE mainTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"';")
        else:
            column_string='CrystalName'+','
            value_string="'"+sampleID+"'"+','
            for key in data_dict:
                value=data_dict[key]
                if key=='ID':
                    continue
                if key not in available_columns:
                    continue
                if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
                    value_string+="'"+str(value)+"'"+','
                    column_string+=key+','
            cursor.execute("INSERT INTO mainTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");")
        connect.commit()



    def update_insert_panddaTable(self,sampleID,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute('Select CrystalName,PANDDA_site_index FROM panddaTable')
#        available_columns=[]
#        cursor.execute("SELECT * FROM panddaTable")
#        for column in cursor.description:           # only update existing columns in data source
#            available_columns.append(column[0])
        samples_sites_in_table=[]
        tmp=cursor.fetchall()
        for item in tmp:
            line=[x.encode('UTF8') for x in list(item)]
            samples_sites_in_table.append(line)

        found_sample_site=False
        for entry in samples_sites_in_table:
            if entry[0]==sampleID and entry[1]==data_dict['PANDDA_site_index']:
                found_sample_site=True

        if found_sample_site:
            for key in data_dict:
                value=data_dict[key]
                if key=='ID' or key=='CrystalName' or key=='PANDDA_site_index':
                    continue
                if not str(value).replace(' ','')=='':  # ignore empty fields
                    update_string=str(key)+'='+"'"+str(value)+"'"
#                    print "UPDATE panddaTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"' and PANDDA_site_index is '"+data_dict['PANDDA_site_index']+"';"
                    cursor.execute("UPDATE panddaTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"' and PANDDA_site_index is '"+data_dict['PANDDA_site_index']+"';")
        else:
            column_string=''
            value_string=''
            for key in data_dict:
                value=data_dict[key]
                if key=='ID':
                    continue
#                if key not in available_columns:
#                    continue
                if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
                    value_string+="'"+str(value)+"'"+','
                    column_string+=key+','
            print "INSERT INTO panddaTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");"
            cursor.execute("INSERT INTO panddaTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");")
        connect.commit()

    def update_insert_any_table(self,table,data_dict,condition_dict):
        data_dict['LastUpdated'] = str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by'] = getpass.getuser()
        connect = sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()

        # columns
        columns = ''
        for c in condition_dict:
            columns+=c+','

        # condition
        condition_string = ''
        for key in condition_dict:
            condition = condition_dict[key]
            condition_string += str(key) + '=' + "'" + str(condition) + "' and "

        cursor.execute('Select %s FROM %s where %s' %(columns[:-1],table,condition_string[:-5]))

        tmp = cursor.fetchall()
        if not tmp:
            data_dict.update(condition_dict)
            value_string=''
            column_string=''
            for key in data_dict:
                value = data_dict[key]
                value_string += "'" + str(value) + "'" + ','
                column_string += key + ','
            cursor.execute("INSERT INTO "+table+" (" + column_string[:-1] + ") VALUES (" + value_string[:-1] + ");")
        else:
            update_string=''
            for key in data_dict:
                value = data_dict[key]
                update_string += str(key) + '=' + "'" + str(value) + "',"
            cursor.execute(
                "UPDATE " + table +
                " SET " + update_string[:-1] +
                " WHERE " + condition_string[:-5] + ";")
        connect.commit()


    def update_insert_site_event_panddaTable(self,sampleID,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute('Select CrystalName,PANDDA_site_index,PANDDA_site_event_index FROM panddaTable')
#        available_columns=[]
#        cursor.execute("SELECT * FROM panddaTable")
#        for column in cursor.description:           # only update existing columns in data source
#            available_columns.append(column[0])
        samples_sites_in_table=[]
        tmp=cursor.fetchall()
        for item in tmp:
            line=[x.encode('UTF8') for x in list(item)]
            samples_sites_in_table.append(line)

        found_sample_site=False
        for entry in samples_sites_in_table:
            if entry[0]==sampleID and entry[1]==data_dict['PANDDA_site_index'] and entry[2]==data_dict['PANDDA_site_event_index']:
                found_sample_site=True

        if found_sample_site:
            for key in data_dict:
                value=data_dict[key]
                if key=='ID' or key=='CrystalName' or key=='PANDDA_site_index' or key=='PANDDA_site_event_index':
                    continue
                if not str(value).replace(' ','')=='':  # ignore empty fields
                    update_string=str(key)+'='+"'"+str(value)+"'"
#                    print "UPDATE panddaTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"' and PANDDA_site_index is '"+data_dict['PANDDA_site_index']+"';"
                    cursor.execute("UPDATE panddaTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"' and PANDDA_site_index is '"+data_dict['PANDDA_site_index']+"' and PANDDA_site_event_index is '"+data_dict['PANDDA_site_event_index']+"';")
        else:
            column_string=''
            value_string=''
            for key in data_dict:
                value=data_dict[key]
                if key=='ID':
                    continue
#                if key not in available_columns:
#                    continue
                if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
                    value_string+="'"+str(value)+"'"+','
                    column_string+=key+','
            print "INSERT INTO panddaTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");"
            cursor.execute("INSERT INTO panddaTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");")
        connect.commit()



    def update_insert_depositTable(self,sampleID,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        available_columns=[]
        cursor.execute("SELECT * FROM depositTable")
        for column in cursor.description:           # only update existing columns in data source
            available_columns.append(column[0])


        cursor.execute('Select CrystalName FROM depositTable')
        samples_in_table=[]
        tmp=cursor.fetchall()
        for item in tmp:
            line=[x.encode('UTF8') for x in list(item)]
#            print 'a',item
#            print 'b',str(item)
#            print 'c',line
#            print 'd',line[0]
#            print 'e',str(line[0])
            if str(item) not in samples_in_table: samples_in_table.append(str(line[0]))
#            samples_in_table.append(line)



        if sampleID in samples_in_table:
            for key in data_dict:
                value=data_dict[key]
                if key=='ID' or key=='CrystalName':
                    continue
                if not str(value).replace(' ','')=='':  # ignore empty fields
                    update_string=str(key)+'='+"'"+str(value)+"'"
                    print "UPDATE depositTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"';"
                    cursor.execute("UPDATE depositTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"';")
        else:
            column_string='CrystalName'+','
            value_string="'"+sampleID+"'"+','
            for key in data_dict:
                value=data_dict[key]
                if key=='ID':
                    continue
                if key not in available_columns:
                    continue
                if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
                    value_string+="'"+str(value)+"'"+','
                    column_string+=key+','
            print "INSERT INTO depositTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");"
            cursor.execute("INSERT INTO depositTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");")
        connect.commit()



    def update_insert_not_null_fields_only(self,sampleID,data_dict):
        data_dict['LastUpdated']=str(datetime.now().strftime("%Y-%m-%d %H:%M"))
        data_dict['LastUpdated_by']=getpass.getuser()
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute('Select CrystalName FROM mainTable')
        available_columns=[]
        cursor.execute("SELECT * FROM mainTable")
        for column in cursor.description:           # only update existing columns in data source
            available_columns.append(column[0])
        if self.check_if_sample_exists_in_data_source(sampleID):
            for key in data_dict:
                value=data_dict[key]
#                print value
                if key=='ID' or key=='CrystalName':
                    continue
                if not str(value).replace(' ','')=='':  # ignore empty fields
                    update_string=str(key)+'='+"'"+str(value)+"'"
#                    cursor.execute("UPDATE mainTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"' and "+str(key)+" is null;")
                    cursor.execute("UPDATE mainTable SET "+update_string+" WHERE CrystalName="+"'"+sampleID+"' and ("+str(key)+" is null or "+str(key)+"='');")
        else:
            column_string='CrystalName'+','
            value_string="'"+sampleID+"'"+','
            for key in data_dict:
                value=data_dict[key]
                if key=='ID':
                    continue
                if key not in available_columns:
                    continue
                if not str(value).replace(' ','')=='':  # ignore if nothing in csv field
                    value_string+="'"+str(value)+"'"+','
                    column_string+=key+','
            cursor.execute("INSERT INTO mainTable ("+column_string[:-1]+") VALUES ("+value_string[:-1]+");")
        connect.commit()

    def get_value_from_field(self,sampleID,column):
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute("SELECT "+column+" FROM  mainTable WHERE CrystalName='"+sampleID+"';")
        return cursor.fetchone()

    def export_to_csv_file(self,csv_file):
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute("SELECT * FROM mainTable")
        header=()
        for column in cursor.description:
            header+=(column[0],)
        rows = cursor.fetchall()
        csvWriter = csv.writer(open(csv_file, "w"))
        csvWriter.writerows([header]+rows)


#    def load_samples_from_data_source(self):
#        header=[]
#        data=[]
#        connect=sqlite3.connect(self.data_source_file)
#        cursor = connect.cursor()
#        cursor.execute("SELECT * FROM mainTable")
#        for column in cursor.description:
#            header.append(column[0])
#        data = cursor.fetchall()
#        return ([header,data])

    def load_samples_from_data_source(self):
        header=[]
        data=[]
        connect=sqlite3.connect(self.data_source_file)
        connect.text_factory = str
        cursor = connect.cursor()
        cursor.execute("SELECT * FROM mainTable")
        for column in cursor.description:
            header.append(column[0])
        data = cursor.fetchall()
        return header,data





#    def get_samples_for_coot(self,RefinementOutcome,pandda_site):
#        sample_list_for_coot=[]
#        connect=sqlite3.connect(self.data_source_file)
#        cursor = connect.cursor()
#
#        if RefinementOutcome=='0 - All Datasets':
#            outcome = " not null "
#        else:
#            outcome = " '%s' " %RefinementOutcome
#
#        if int(pandda_site) > 0:
#            sqlite = (
#                "select"
#                " mainTable.CrystalName,"
#                " mainTable.CompoundCode,"
#                " mainTable.RefinementCIF,"
#                " mainTable.RefinementMTZfree,"
#                " mainTable.RefinementPathToRefinementFolder,"
#                " panddaTable.RefinementOutcome, "
#                " panddaTable.PANDDA_site_event_map,"
#                " panddaTable.PANDDA_site_confidence,"
#                " panddaTable.PANDDA_site_x,"
#                " panddaTable.PANDDA_site_y,"
#                " panddaTable.PANDDA_site_z,"
#                " panddaTable.PANDDA_site_initial_model,"
#                " panddaTable.PANDDA_site_initial_mtz,"
#                " panddaTable.PANDDA_site_spider_plot, "
#                " panddaTable.PANDDA_site_index "
#                "from mainTable inner join panddaTable on mainTable.CrystalName = panddaTable.CrystalName "
#                "where panddaTable.PANDDA_site_index is '%s'" %pandda_site+
#                " and panddaTable.PANDDA_site_ligand_placed is 'True'"
#                " and panddaTable.RefinementOutcome is %s;" %outcome
#                )
#
#        else:
#            sqlite = (
#                "select"
#                " CrystalName,"
#                " CompoundCode,"
#                " RefinementCIF,"
#                " RefinementMTZfree,"
#                " RefinementPathToRefinementFolder,"
#                " RefinementOutcome "
#                "from mainTable "
#                "where RefinementOutcome is %s;" %outcome
#                )
#
#        cursor.execute(sqlite)
#
#        tmp = cursor.fetchall()
#        for item in tmp:
#            tmpx=[]
#            for i in list(item):
#                if i==None:
#                    tmpx.append('None')
#                else:
#                    tmpx.append(i)
#            line=[x.encode('UTF8') for x in tmpx]
#            sample_list_for_coot.append(line)
#
#        return sample_list_for_coot

    def get_pandda_info_for_coot(self,xtalID,pandda_site):
        pandda_info_for_coot=[]
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()

        sqlite = (
            "select"
            " PANDDA_site_event_map,"
            " PANDDA_site_x,"
            " PANDDA_site_y,"
            " PANDDA_site_z,"
            " PANDDA_site_spider_plot "
            "from panddaTable  "
            "where "
            " CrystalName is '%s'" %xtalID+
            " and PANDDA_site_index is '{0!s}';".format(pandda_site)
            )

        cursor.execute(sqlite)

        tmp = cursor.fetchall()
        for item in tmp:
            tmpx=[]
            for i in list(item):
                if i is None:
                    tmpx.append('None')
                else:
                    tmpx.append(i)
            line=[x.encode('UTF8') for x in tmpx]
            pandda_info_for_coot.append(line)

        return pandda_info_for_coot



    def get_todo_list_for_coot(self,RefinementOutcome,pandda_site):
        sample_list_for_coot=[]
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()

        if RefinementOutcome=='0 - All Datasets':
            outcome = " not null "
        else:
            outcome = " '{0!s}' ".format(RefinementOutcome)

        if int(pandda_site) > 0:
#            sqlite = (
#                "select"
#                " mainTable.CrystalName,"
#                " mainTable.CompoundCode,"
#                " mainTable.RefinementCIF,"
#                " mainTable.RefinementMTZfree,"
#                " mainTable.RefinementPathToRefinementFolder,"
#                " panddaTable.RefinementOutcome, "
#                " panddaTable.PANDDA_site_confidence "
#                "from mainTable inner join panddaTable on mainTable.CrystalName = panddaTable.CrystalName "
#                "where panddaTable.PANDDA_site_index is '%s'" %pandda_site+
#                " and panddaTable.PANDDA_site_ligand_placed is 'True'"
#                " and panddaTable.RefinementOutcome is %s;" %outcome
#                )
            sqlite = (
                "select"
                " mainTable.CrystalName,"
                " mainTable.CompoundCode,"
                " mainTable.RefinementCIF,"
                " mainTable.RefinementMTZfree,"
                " mainTable.RefinementPathToRefinementFolder,"
                " panddaTable.RefinementOutcome, "
                " panddaTable.PANDDA_site_confidence "
                "from mainTable inner join panddaTable on mainTable.CrystalName = panddaTable.CrystalName "
                "where panddaTable.PANDDA_site_index is '%s'" %pandda_site+
                " and panddaTable.PANDDA_site_ligand_placed is 'True'"
                " and panddaTable.RefinementOutcome like "+outcome.split()[0]+"%';"
                )
#            sqlite = (
#                "select"
#                " mainTable.CrystalName,"
#                " mainTable.CompoundCode,"
#                " mainTable.RefinementCIF,"
#                " mainTable.RefinementMTZfree,"
#                " mainTable.RefinementPathToRefinementFolder,"
#                " panddaTable.RefinementOutcome, "
#                " panddaTable.PANDDA_site_confidence "
#                "from mainTable inner join panddaTable on mainTable.CrystalName = panddaTable.CrystalName "
#                "where panddaTable.PANDDA_site_index is '%s'" %pandda_site+
#                " and panddaTable.RefinementOutcome like "+outcome.split()[0]+"%';"
#                )
        else:
            sqlite = (
                "select"
                " CrystalName,"
                " CompoundCode,"
                " RefinementCIF,"
                " RefinementMTZfree,"
                " RefinementPathToRefinementFolder,"
                " RefinementOutcome,"
                " RefinementLigandConfidence "
                "from mainTable "
                "where RefinementOutcome is %s and DimpleRfree is not Null;" %outcome
                )
        cursor.execute(sqlite)

        tmp = cursor.fetchall()
        for item in tmp:
            tmpx=[]
            for i in list(item):
                if i is None:
                    tmpx.append('None')
                else:
                    tmpx.append(i)
            line=[x.encode('UTF8') for x in tmpx]
            sample_list_for_coot.append(line)

        return sample_list_for_coot


    def get_todoList_for_coot(self,RefinementOutcome):
        sample_list_for_coot=[]
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()

        if RefinementOutcome=='0 - All Datasets':
            outcome = " not null "
        else:
            outcome = " '{0!s}' ".format(RefinementOutcome)

        sqlite = (
            "select"
            " CrystalName,"
            " CompoundCode,"
            " RefinementCIF,"
            " RefinementMTZfree,"
            " RefinementPathToRefinementFolder,"
            " RefinementOutcome,"
            " RefinementLigandConfidence "
            "from mainTable "
            "where RefinementOutcome is %s and (DimpleRfree is not Null or RefinementRfree is not Null) order by CrystalName ASC;" %outcome
            )

        cursor.execute(sqlite)

        tmp = cursor.fetchall()
        for item in tmp:
            tmpx=[]
            for i in list(item):
                if i is None:
                    tmpx.append('None')
                else:
                    tmpx.append(i)
            line=[x.encode('UTF8') for x in tmpx]
            sample_list_for_coot.append(line)

        crystalDict={}
        for entry in sample_list_for_coot:
            if entry[0] not in crystalDict:
                sqlite = (  "select"
                            " PANDDA_site_event_map,"
                            " PANDDA_site_x,"
                            " PANDDA_site_y,"
                            " PANDDA_site_z,"
                            " PANDDA_site_spider_plot,"
                            " PANDDA_site_index,"
                            " PANDDA_site_event_index,"
                            " PANDDA_site_confidence,"
                            " PANDDA_site_name,"
                            " PANDDA_site_InspectConfidence,"
                            " PANDDA_site_interesting,"
                            " PANDDA_site_event_comment "
                            "from panddaTable  "
                            "where "
                            " CrystalName is '%s';" %entry[0]     )
                cursor.execute(sqlite)
                tmp = cursor.fetchall()
                if tmp:
                    crystalDict[entry[0]]=[]
                    for item in tmp:
                        print [entry[0], str(item[0]),str(item[1]),str(item[2]),str(item[3]),str(item[4]),str(item[5]),str(item[6]),str(item[7]),str(item[8]),str(item[9]),str(item[10]),str(item[11]) ]
                        crystalDict[entry[0]].append( [ str(item[0]),str(item[1]),str(item[2]),str(item[3]),str(item[4]),str(item[5]),str(item[6]),str(item[7]),str(item[8]),str(item[9]),str(item[10]),str(item[11]) ])

        return sample_list_for_coot,crystalDict




    def translate_xce_column_list_to_sqlite(self,column_list):
        out_list=[]
        for item in column_list:
            if item.startswith('Exclude'):
                out_list.append(['Exclude'])
            if item.startswith('Ignore'):
                out_list.append(['Ignore'])
#            if item.startswith('img'):
#                out_list.append([item,item])
#                continue
            if item.startswith('Show'):
                out_list.append([item,item])
                continue
            if item.startswith('Run\nDimple'):
                out_list.append([item,item])
                continue
            if item.startswith('Select'):
                out_list.append([item,item])
                continue
            if item.startswith('Run\nxia2'):
                out_list.append([item,item])
                continue
            if item.startswith('Dataset ID'):
                out_list.append([item,item])
                continue
            if item.startswith('Reference\nSpaceGroup'):
                out_list.append([item,item])
                continue
            if item.startswith('Difference\nUC Volume (%)'):
                out_list.append([item,item])
                continue
            if item.startswith('Reference File'):
                out_list.append([item,item])
                continue
            if item.startswith('PanDDA site details'):
                out_list.append([item,item])
                continue
            for entry in self.column_list:
                if entry[1]==item:
                    out_list.append([item,entry[0]])
                    break
        return out_list

    def get_list_of_pandda_sites_for_coot(self):
        site_list=[ ['0','any site'] ]
        sqlite = (
            'select distinct'
            ' panddaTable.PANDDA_site_index,'
            ' panddaTable.PANDDA_site_name '
            'from panddaTable '
            'order by cast (panddaTable.PANDDA_site_index as integer) ASC;'
                  )

        print 'data source',self.data_source_file
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute(sqlite)
        tmp=cursor.fetchall()
        for item in tmp:
#            print item
#            print str(item[0])
#            print str(item[1])
#            line=[x.encode('UTF8') for x in list(item)]
            site_list.append([str(item[0]),str(item[1])])

        return site_list

    def export_csv_for_WONKA(self):
        SQLite =    (   "select"
                        " panddaTable.PANDDA_site_ligand_resname,"
                        " panddaTable.PANDDA_site_ligand_sequence_number,"
                       	" panddaTable.PANDDA_site_ligand_chain,"
                        " panddaTable.RefinementOutcome,"
                        " mainTable.CompoundSMILES,"
                        " mainTable.RefinementBoundConformation,"
                        " panddaTable.PANDDA_site_initial_mtz,"
                        " panddaTable.PANDDA_site_event_map,"
                        " panddaTable.CrystalName,"
                        " mainTable.CompoundCode "
                        "from"
                        " mainTable inner join panddaTable on mainTable.CrystalName = panddaTable.CrystalName "
                        "where"
                        " (panddaTable.RefinementOutcome like '3%' or panddaTable.RefinementOutcome like '4%' or panddaTable.RefinementOutcome like '5%')" )
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute(SQLite)
        header=()
        for column in cursor.description:
            header+=(column[0],)
        rows = cursor.fetchall()
        csvWriter = csv.writer(open('for_wonka.csv', "w"))
        csvWriter.writerows([header]+rows)

    def create_missing_apo_records_in_depositTable(self,xce_logfile):
        Logfile=XChemLog.updateLog(xce_logfile)
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute("select panddaTable.ApoStructures from panddaTable where panddaTable.ApoStructures is not Null")
        tmp = cursor.fetchall()
#        print str(tmp)
        apoStructureList=[]
        for item in tmp:
            for xtal in str(item[0]).split(';'):
                if xtal not in apoStructureList: apoStructureList.append(xtal)
        if not apoStructureList:
            Logfile.insert('no apo structures were assigned to your pandda models')
        else:
            Logfile.insert('the following datasets were at some point used as apo structures for pandda.analyse: '+str(apoStructureList))
            apoInDB=[]
            cursor.execute("select CrystalName from depositTable where StructureType is 'apo'")
            tmp = cursor.fetchall()
            for xtal in tmp:
                Logfile.insert(str(xtal[0])+' exists as entry for apo structure in database')
                apoInDB.append(str(xtal[0]))
            newEntries=''
            for xtal in apoStructureList:
                if xtal not in apoInDB:
                    Logfile.insert('no entry for '+xtal+' in depositTable')
                    newEntries+="('{0!s}','apo'),".format(xtal)
            if newEntries != '':
                sqlite='insert into depositTable (CrystalName,StructureType) values {0!s};'.format(newEntries[:-1])
                Logfile.insert('creating new entries with the following SQLite command:\n'+sqlite)
                cursor.execute(sqlite)
                connect.commit()

    def create_missing_apo_records_for_all_structures_in_depositTable(self,projectDir,xce_logfile):
        Logfile=XChemLog.updateLog(xce_logfile)
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        cursor.execute("select CrystalName from depositTable where StructureType is 'apo'")
        tmp = cursor.fetchall()
#        print str(tmp)
        apoStructureList=[]
        for item in tmp:
            apoStructureList.append(str(item[0]))

        counter=0   # there is a SQLite limit which does not allow to insert more than 500 records at once
        newEntries=''
        for files in glob.glob(os.path.join(projectDir,'*')):
            xtal=files[files.rfind('/')+1:]
            if xtal not in apoStructureList:
                dimple_pdb=False
                dimple_mtz=False
                aimless_log=False
                if os.path.isfile(os.path.join(files,'dimple.pdb')):
                    dimple_pdb=True
                if os.path.isfile(os.path.join(files,'dimple.mtz')):
                    dimple_mtz=True
                if os.path.isfile(os.path.join(files,xtal+'.log')) or os.path.isfile(os.path.join('aimless.log')):
                    aimless_log=True

                if dimple_mtz==False or dimple_pdb==False or aimless_log==False:
                    Logfile.insert('{0!s}: missing files -> dimple.pdb: {1!s}, dimple.mtz: {2!s}, {3!s}.log: {4!s}'.format(xtal, str(dimple_pdb), str(dimple_mtz), xtal, str(aimless_log)))
                else:
                    newEntries+="('{0!s}','apo'),".format(xtal)
                    counter+=1

                if counter == 450:  # set to 450 to stay well below 500 records limit
                    sqlite='insert into depositTable (CrystalName,StructureType) values {0!s};'.format(newEntries[:-1])
                    Logfile.insert('creating new entries with the following SQLite command:\n'+sqlite)
                    cursor.execute(sqlite)
                    connect.commit()
                    counter=0
                    newEntries=''

        if newEntries != '':
            sqlite='insert into depositTable (CrystalName,StructureType) values {0!s};'.format(newEntries[:-1])
            Logfile.insert('creating new entries with the following SQLite command:\n'+sqlite)
            cursor.execute(sqlite)
            connect.commit()
        else:
            Logfile.insert('did not find any new apo structures')

#        if apoStructureList==[]:
#            Logfile.insert('no apo structures were assigned to your pandda models')
#        else:
#            Logfile.insert('the following datasets were at some point used as apo structures for pandda.analyse: '+str(apoStructureList))
#            apoInDB=[]
#            cursor.execute("select CrystalName from depositTable where StructureType is 'apo'")
#            tmp = cursor.fetchall()
#            for xtal in tmp:
#                Logfile.insert(str(xtal[0])+' exists as entry for apo structure in database')
#                apoInDB.append(str(xtal[0]))
#            newEntries=''
#            for xtal in apoStructureList:
#                if xtal not in apoInDB:
#                    Logfile.insert('no entry for '+xtal+' in depositTable')
#                    newEntries+="('%s','apo')," %xtal
#            if newEntries != '':
#                sqlite='insert into depositTable (CrystalName,StructureType) values %s;' %newEntries[:-1]
#                Logfile.insert('creating new entries with the following SQLite command:\n'+sqlite)
#                cursor.execute(sqlite)
#                connect.commit()

#    def get_deposit_dict_for_xtal(self,sampleID,structure_type):

    def create_or_remove_missing_records_in_depositTable(self,xce_logfile,xtal,type,db_dict):
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()
        Logfile=XChemLog.updateLog(xce_logfile)
        oldRefiStage = ''
        if type == 'ligand_bound':
            cursor.execute("select RefinementOutcome from mainTable where CrystalName is '{0!s}'".format(xtal))
            tmp = cursor.fetchall()
            oldRefiStage=str(tmp[0][0])
            Logfile.insert('setting RefinementOutcome field to "'+db_dict['RefinementOutcome']+'" for '+xtal)
            self.update_insert_data_source(xtal,db_dict)
        elif type == 'ground_state':
            cursor.execute("select DimplePANDDApath from depositTable where StructureType is '{0!s}' and DimplePANDDApath is '{1!s}'".format(type,db_dict['DimplePANDDApath']))
            tmp = cursor.fetchall()
            if not tmp:
                Logfile.insert('entry for ground-state model in depositTable does not exist')
            else:
                Logfile.warning('entry for ground-state model in depositTable does already exist')
                Logfile.warning('updating PDB, MTZ and DimplePANDDApath for ground-state entry')
                cursor.execute("update depositTable set PDB_file='%s', MTZ_file='%s', DimplePANDDApath='%s' where StructureType='ground_state'" %(db_dict['PDB_file'],db_dict['MTZ_file'],db_dict['DimplePANDDApath']))
                connect.commit()
                return

        cursor.execute("select CrystalName,StructureType from depositTable where CrystalName is '{0!s}' and StructureType is '{1!s}'".format(xtal, type))
        tmp = cursor.fetchall()
        if type == 'ligand_bound':
            if not tmp and int(db_dict['RefinementOutcome'].split()[0]) == 5:
                sqlite="insert into depositTable (CrystalName,StructureType) values ('{0!s}','{1!s}');".format(xtal, type)
                Logfile.insert('creating new entry for '+str(type)+' structure of '+xtal+' in depositTable')
                cursor.execute(sqlite)
                connect.commit()
            else:
                if int(db_dict['RefinementOutcome'].split()[0]) != 5:
                    sqlite="delete from depositTable where CrystalName is '{0!s}' and StructureType is '{1!s}'".format(xtal, type)
                    if oldRefiStage.startswith('5'):
                        Logfile.insert('removing entry for '+str(type)+' structure of '+xtal+' from depositTable')
                    cursor.execute(sqlite)
                    connect.commit()
        elif type == 'ground_state':
            sqlite = "insert into depositTable (CrystalName,StructureType,DimplePANDDApath,PDB_file,MTZ_file) values ('{0!s}','{1!s}','{2!s}','{3!s}','{4!s}');".format(xtal, type, db_dict['DimplePANDDApath'], db_dict['PDB_file'], db_dict['MTZ_file'])
            Logfile.insert('creating new entry for ' + str(type) + ' structure of ' + xtal + ' in depositTable')
            cursor.execute(sqlite)
            connect.commit()

    def remove_selected_apo_structures_from_depositTable(self,xce_logfile,xtalList):
        connect=sqlite3.connect(self.data_source_file)
        cursor = connect.cursor()

        Logfile=XChemLog.updateLog(xce_logfile)

        Logfile.insert('removing the following apo structures from depositTable')
        deleteStrg = "("
        for xtal in xtalList:
            Logfile.insert(xtal)
            deleteStrg+="Crystalname is '"+xtal+"' or "
        deleteStrg=deleteStrg[:-4]+")"

        sqlite="delete from depositTable where {0!s} and StructureType is 'apo'".format(deleteStrg)
        Logfile.insert('executing the following SQLite command:\n'+sqlite)
        cursor.execute(sqlite)
        connect.commit()




#        for item in tmp:
#            print item
#            for xtal in str(item[0]).split(';'):
#                if xtal not in apoStructureList: apoStructureList.append(xtal)
#        if apoStructureList==[]:
#            Logfile.insert('no apo structures were assigned to your pandda models')
#        else:
#            Logfile.insert('the following datasets were at some point used as apo structures for pandda.analyse: '+str(apoStructureList))
#            apoInDB=[]
#            cursor.execute("select CrystalName from depositTable where StructureType is 'apo'")
#            tmp = cursor.fetchall()
#            for xtal in tmp:
#                Logfile.insert(str(xtal[0])+' exists as entry for apo structure in database')
#                apoInDB.append(str(xtal[0]))
#            newEntries=''
#            for xtal in apoStructureList:
#                if xtal not in apoInDB:
#                    Logfile.insert('no entry for '+xtal+' in depositTable')
#                    newEntries+="('%s','apo')," %xtal
#            if newEntries != '':
#                sqlite='insert into depositTable (CrystalName,StructureType) values %s;' %newEntries[:-1]
#                Logfile.insert('creating new entries with the following SQLite command:\n'+sqlite)
#                cursor.execute(sqlite)
#                connect.commit()
#

    def collected_xtals_during_visit(self,visitID):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("select CrystalName from collectionTable where DataCollectionVisit = '{0!s}'".format(visitID))
        collectedXtals=[]
        samples = cursor.fetchall()
        for sample in samples:
            if str(sample[0]) not in collectedXtals:
                collectedXtals.append(str(sample[0]))
        return collectedXtals

    def collected_xtals_during_visit_for_scoring(self,visit):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        xtalList = []
#        if rescore == True:
#            cursor.execute("select CrystalName from mainTable where DataCollectionVisit = '%s'" %visit)
#        else:
#            cursor.execute("select CrystalName from mainTable where DataProcessingAutoAssigned = 'True' and DataCollectionVisit = '%s'" %visit)
        cursor.execute("select distinct CrystalName from collectionTable where DataCollectionVisit = '%s'" %visit)
        samples = cursor.fetchall()
        for sample in samples:
            if str(sample[0]) not in xtalList:
                xtalList.append(str(sample[0]))
        return xtalList

    def autoprocessing_result_user_assigned(self,sample):
        userassigned = False
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("select DataProcessingAutoAssigned from mainTable where CrystalName = '%s'" %sample)
        outcome = cursor.fetchall()
        try:
            if 'false' in str(outcome[0]).lower() or str(outcome[0]).encode('ascii', 'ignore') == '':
                userassigned = True        # a bit counterintuitive, but here we ask about userassigned
                                            # whereas DB records autoassigned
        except IndexError:
            pass                # this is the case when sample is encountered the first time
                                # and not yet in mainTable

        return userassigned


    def all_results_of_xtals_collected_during_visit_as_dict(self,visitID):
        # first get all collected xtals as list
        collectedXtals = self.collected_xtals_during_visit(visitID)
        xtalDict = {}
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        for xtal in sorted(collectedXtals):
            if xtal not in xtalDict:
                xtalDict[xtal] = []
            sqlite = ("select DataCollectionVisit,"
                      "       DataCollectionRun,"
                      "       DataProcessingProgram,"
                      "       DataCollectionDate,"
                      "       DataProcessingResolutionHigh,"
                      "       DataProcessingRmergeLow,"
                      "       DataCollectionCrystalImage1,"
                      "       DataCollectionCrystalImage2,"
                      "       DataCollectionCrystalImage3,"
                      "       DataCollectionCrystalImage4 "
                      "from collectionTable "
                      "where CrystalName = '%s'" % xtal)
#            print sqlite
            cursor.execute(sqlite)
            sq = cursor.fetchall()
#            print sq
            for s in sq:
#                print s
                t=[str(s[0]),str(s[1]),str(s[2]),str(s[3]),str(s[4]),str(s[5]),str(s[6]),str(s[7]),str(s[8]),str(s[9])]
                inDict = False
                for entry in xtalDict[xtal]:
                    if t[0] == entry[0] and t[1] == entry[1] and t[2] == entry[2]:
                        inDict = True
                        break
                if not inDict:
                    xtalDict[xtal].append(t)
        return xtalDict

    def all_autoprocessing_results_for_xtal_as_dict(self,xtal):
        dbList = []
        header=[]
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("select * from collectionTable where CrystalName='{0!s}';".format(xtal))
        for column in cursor.description:
            header.append(column[0])
        data = cursor.fetchall()
        for result in data:
            db_dict = {}
            for n,item in enumerate(result):
                db_dict[header[n]]=str(item)
            dbList.append(db_dict)
        return dbList

    def get_db_list_for_sample_in_panddaTable(self,xtal):
        dbList = []
        header=[]
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("select * from panddaTable where CrystalName='{0!s}';".format(xtal))
        for column in cursor.description:
            header.append(column[0])
        data = cursor.fetchall()
        for result in data:
            db_dict = {}
            for n,item in enumerate(result):
                db_dict[header[n]]=str(item)
            dbList.append(db_dict)
        return dbList

    def get_db_dict_for_visit_run_autoproc(self,xtal,visit,run,autoproc):
        db_dict = {}
        header=[]
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        sqlite = (  'select * '
                    'from collectionTable '
                    "where CrystalName ='%s' and" %xtal +
                    "      DataCollectionVisit = '%s' and" %visit +
                    "      DataCollectionRun = '%s' and" %run +
                    "      DataProcessingProgram = '%s'" %autoproc  )
        cursor.execute(sqlite)
        for column in cursor.description:
            header.append(column[0])
        data = cursor.fetchall()
        for n, item in enumerate(data[0]):
            db_dict[header[n]] = str(item)
        return db_dict

    def xtals_collected_during_visit_as_dict(self,visitID):
        # first get all collected xtals as list
        if isinstance(visitID,list):    # for Agamemnon data structure
            collectedXtals = []
            for visit in visitID:
                x = self.collected_xtals_during_visit(visit)
                for e in x:
                    collectedXtals.append(e)
        else:
            collectedXtals = self.collected_xtals_during_visit(visitID)
        xtalDict = {}
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        for xtal in sorted(collectedXtals):
            db_dict = self.get_db_dict_for_sample(xtal)
            xtalDict[xtal] = db_dict
        return xtalDict



    def getCrystalImageDict(self,visit,xtalList):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        imageDict = {}
        for xtal in xtalList:
            sqlite = (  "select CrystalName,"
		                "   min(DataCollectionDate),"
		                "   DataCollectionCrystalImage1,"
		                "   DataCollectionCrystalImage2,"
    		            "   DataCollectionCrystalImage3,"
	    	            "   DataCollectionCrystalImage4 "
                        "from collectionTable "
                        "where DataCollectionVisit = '%s'"  %visit+
                        " and CrystalName = '%s'"           % xtal      )
            cursor.execute(sqlite)
            img = cursor.fetchall()
            imageDict[xtal]=[str(img[1]),str(img[2]),str(img[3]),str(img[4])]
        return imageDict

    def samples_for_html_summary(self):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("select CrystalName from mainTable where RefinementOutcome like '4%' or "
                                                               "RefinementOutcome like '5%' or "
                                                               "RefinementOutcome like '6%' order by CrystalName ASC")
        xtalList=[]
        samples = cursor.fetchall()
        for sample in samples:
            if str(sample[0]) not in xtalList:
                xtalList.append(str(sample[0]))
        return xtalList

    def get_event_map_for_ligand(self,xtal,ligChain,ligNumber,ligName):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()

        sql = (
            'select '
            ' PANDDA_site_event_map '
            'from '
            ' panddaTable '
            'where '
            " CrystalName = '%s' and " %xtal +
            " PANDDA_site_ligand_chain='%s' and " %ligChain +
            " PANDDA_site_ligand_sequence_number='%s' and " %ligNumber +
            " PANDDA_site_ligand_resname='%s'"  %ligName
        )

        cursor.execute(sql)

        eventMap = ''
        maps = cursor.fetchall()
        for map in maps:
            eventMap = map[0]

        return eventMap

    def get_ligand_confidence_for_ligand(self,xtal,ligChain,ligNumber,ligName):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()

        sql = (
            'select '
            ' PANDDA_site_confidence '
            'from '
            ' panddaTable '
            'where '
            " CrystalName = '%s' and " %xtal +
            " PANDDA_site_ligand_chain='%s' and " %ligChain +
            " PANDDA_site_ligand_sequence_number='%s' and " %ligNumber +
            " PANDDA_site_ligand_resname='%s'"  %ligName
        )

        cursor.execute(sql)

        ligConfidence = 'not assigned'
        ligs = cursor.fetchall()
        for lig in ligs:
            ligConfidence = lig[0]

        return ligConfidence

    def get_labels_from_db(self):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute('select Label from labelTable')
        labels = cursor.fetchall()
        labelList = []
        for l in labels:
            labelList.append(l[0])
        return labelList

    def get_label_info_from_db(self):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute('select Label,Description from labelTable')
        labels = cursor.fetchall()
        labelList = []
        for l in labels:
            labelList.append([str(l[0]),str(l[1])])
        return labelList

    def get_label_of_sample(self,xtalID):
        connect=sqlite3.connect(self.data_source_file)     # creates sqlite file if non existent
        cursor = connect.cursor()
        cursor.execute("select label from mainTable where CrystalName is '%s'" %xtalID)
        label = None
        lab = cursor.fetchall()
        for l in lab:
            label = l[0]
            break
        return label
