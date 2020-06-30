# last edited: 10/08/2017, 10:25

import os, sys, glob, subprocess
from datetime import datetime
from PyQt4 import QtGui, QtCore
import math
#from XChemUtils import mtztools
import XChemDB
import XChemRefine
import XChemUtils
import XChemLog
import XChemToolTips
import csv

try:
    import gemmi
    import pandas
except ImportError:
    pass

#def get_names_of_current_clusters(xce_logfile,panddas_directory):
#    Logfile=XChemLog.updateLog(xce_logfile)
#    Logfile.insert('parsing {0!s}/cluster_analysis'.format(panddas_directory))
#    os.chdir('{0!s}/cluster_analysis'.format(panddas_directory))
#    cluster_dict={}
#    for out_dir in sorted(glob.glob('*')):
#        if os.path.isdir(out_dir):
#            cluster_dict[out_dir]=[]
#            found_first_pdb=False
#            for folder in glob.glob(os.path.join(out_dir,'pdbs','*')):
#                xtal=folder[folder.rfind('/')+1:]
#                if not found_first_pdb:
#                    if os.path.isfile(os.path.join(panddas_directory,'cluster_analysis',out_dir,'pdbs',xtal,xtal+'.pdb') ):
#                        cluster_dict[out_dir].append(os.path.join(panddas_directory,'cluster_analysis',out_dir,'pdbs',xtal,xtal+'.pdb'))
#                        found_first_pdb=True
#                cluster_dict[out_dir].append(xtal)
#    return cluster_dict



class export_and_refine_ligand_bound_models(QtCore.QThread):

    def __init__(self,PanDDA_directory,datasource,project_directory,xce_logfile,which_models):
        QtCore.QThread.__init__(self)
        self.PanDDA_directory = PanDDA_directory
        self.datasource = datasource
        self.db = XChemDB.data_source(self.datasource)
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.project_directory = project_directory
        self.which_models=which_models

#        self.initial_model_directory=initial_model_directory
#        self.db.create_missing_columns()
#        self.db_list=self.db.get_empty_db_dict()
#        self.external_software=XChemUtils.external_software(xce_logfile).check()
#        self.xce_logfile=xce_logfile

#        self.already_exported_models=[]

    def run(self):

        self.Logfile.warning(XChemToolTips.pandda_export_ligand_bound_models_only_disclaimer())

        # find all folders with *-pandda-model.pdb
        modelsDict = self.find_modeled_structures_and_timestamps()

        # if only NEW models shall be exported, check timestamps
#        if self.which_models != 'all':
#            modelsDict = self.find_new_models(modelsDict)

        # find pandda_inspect_events.csv and read in as pandas dataframe
        inspect_csv = None
        if os.path.isfile(os.path.join(self.PanDDA_directory,'analyses','pandda_inspect_events.csv')):
            inspect_csv = pandas.read_csv(os.path.join(self.PanDDA_directory,'analyses','pandda_inspect_events.csv'))

        progress = 0
        try:
            progress_step = float(1/len(modelsDict))
        except TypeError:
            self.Logfile.error('DID NOT FIND ANY MODELS TO EXPORT')
            return None

        for xtal in modelsDict:
            os.chdir(os.path.join(self.PanDDA_directory,'processed_datasets',xtal))
            pandda_model = os.path.join('modelled_structures',xtal + '-pandda-model.pdb')
            pdb = gemmi.read_structure(pandda_model)

            # find out ligand event map relationship
            ligandDict = XChemUtils.pdbtools_gemmi(pandda_model).center_of_mass_ligand_dict('LIG')
            self.show_ligands_in_model(xtal,ligandDict)
            emapLigandDict = self.find_ligands_matching_event_map(inspect_csv,xtal,ligandDict)

            # convert event map to SF
#            self.event_map_to_sf(pdb.resolution,emapLigandDict)
            break
            # move existing event maps in project directory to old folder

            # update database

            # copy files

            # refine models

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

    def show_ligands_in_model(self,xtal,ligandDict):
        self.Logfile.insert(xtal + ': found the following ligands...')
        for lig in ligandDict:
            self.Logfile.insert(lig + ' -> coordinates ' + str(ligandDict[lig]))


    def find_modeled_structures_and_timestamps(self):
        self.Logfile.insert('finding out modelled structures in ' + self.PanDDA_directory)
        modelsDict={}
        for model in glob.glob(os.path.join(self.PanDDA_directory,'processed_datasets','*','modelled_structures','*-pandda-model.pdb')):
            sample=model[model.rfind('/')+1:].replace('-pandda-model.pdb','')
            timestamp=datetime.fromtimestamp(os.path.getmtime(model)).strftime('%Y-%m-%d %H:%M:%S')
            self.Logfile.insert(sample+'-pandda-model.pdb was created on '+str(timestamp))
            modelsDict[sample]=timestamp
        return modelsDict

    def find_new_models(self,modelsDict):
        for xtal in modelsDict:
            timestamp_file = modelsDict[xtal]
            db_query=self.db.execute_statement("select DatePanDDAModelCreated from mainTable where CrystalName is '"+xtal+"' and (RefinementOutcome like '3%' or RefinementOutcome like '4%' or RefinementOutcome like '5%')")
            timestamp_db=str(db_query[0][0])
            try:
                difference=(datetime.strptime(timestamp_file,'%Y-%m-%d %H:%M:%S') - datetime.strptime(timestamp_db,'%Y-%m-%d %H:%M:%S')  )
                if difference.seconds != 0:
                    self.Logfile.insert('exporting '+sample+' -> was already refined, but newer PanDDA model available')
                    samples_to_export[sample]=fileModelsDict[sample]
            except (ValueError, IndexError), e:
                self.Logfile.error(str(e))

    def event_map_to_sf(self,resolution,emapLigandDict):
        emtz = emap.replace('.ccp4','.mtz')
        emtz_ligand = emtz
        for emap in glob.glob('*-BDC_*.ccp4'):
            if emap in emapLigandDict:
                emtz_ligand = emap.replace('.ccp4','_' + emapLigandDict[emap] + '.mtz')
            self.Logfile.insert('trying to convert %s to SF' %emap)
            self.Logfile.insert('>>> ' + emtz)
            XChemUtils.maptools_gemmi(emap).map_to_sf(resolution)
            if os.path.isfile(emtz):
                print('/bin/mv %s %s' %(emtz,emtz_ligand))
                self.Logfile.insert('success; %s exists' %emtz_ligand)
            else:
                self.Logfile.warning('something went wrong; %s could not be created...' %emtz_ligand)

    def find_ligands_matching_event_map(self,inspect_csv,xtal,ligandDict):
        emapLigandDict = {}
        for index, row in inspect_csv.iterrows():
            if row['dtag'] == xtal:
#                site = str(row['site_idx'])
                site = str(row['Unnamed: 0'])
                event = str(row['event_idx'])
                for emap in glob.glob('*-BDC_*.ccp4'):
                    site_emap = emap[emap.find('event')+6:emap.find('BDC')-1].split('_')[0]
                    event_emap = emap[emap.find('event')+6:emap.find('BDC')-1].split('_')[1]
                    self.Logfile.insert(emap+ ' site: '+site+' esite: '+site_emap+' event: '+event+' eevent: '+event_emap)
                    if site == site_emap and event == event_emap:
                        self.Logfile.insert('found event map for site/event')
                        x = float(row['x'])
                        y = float(row['y'])
                        z = float(row['z'])
                        matching_ligand = self.calculate_distance_to_ligands(ligandDict,x,y,z)
                        emapLigandDict[emap] = matching_ligand
                        self.Logfile.insert('found matching ligand (%s) for %s' %(matching_ligand,emap))
                        break
        if emapLigandDict == {}:
            self.Logfile.error('could not find ligands within 5A of PanDDA events')
        return emapLigandDict

    def calculate_distance_to_ligands(self,ligandDict,x,y,z):
        matching_ligand = ''
        p_event = gemmi.Position(x, y, z)
        for ligand in ligandDict:
            c = ligandDict[ligand]
            p_ligand = gemmi.Position(c[0], c[1], c[2])
            self.Logfile.insert('lig:   ' + str(c[0])+' '+ str(c[1])+' '+str(c[2]))
            self.Logfile.insert('event: ' + str(x)+' '+ str(y)+' '+str(z))
            distance = p_event.dist(p_ligand)
            if distance < 5:
                matching_ligand = ligand
                break
        return matching_ligand

# mArh-x0969-event_1_1-BDC_0.67_map.ccp4
# site = x[x.find('event')+6:x.find('BDC')-1].split('_')[0]
# event = x[x.find('event')+6:x.find('BDC')-1].split('_')[1]










class refine_bound_state_with_buster(QtCore.QThread):

    def __init__(self,panddas_directory,datasource,initial_model_directory,xce_logfile,which_models):
        QtCore.QThread.__init__(self)
        self.panddas_directory=panddas_directory
        self.datasource=datasource
        self.initial_model_directory=initial_model_directory
        self.db=XChemDB.data_source(self.datasource)
        self.db.create_missing_columns()
        self.db_list=self.db.get_empty_db_dict()
        self.external_software=XChemUtils.external_software(xce_logfile).check()
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.which_models=which_models
        self.already_exported_models=[]

    def run(self):

        samples_to_export=self.export_models()

        self.refine_exported_models(samples_to_export)


    def refine_exported_models(self,samples_to_export):
        self.Logfile.insert('will try to refine the following crystals:')
        for xtal in sorted(samples_to_export):
            self.Logfile.insert(xtal)
        for xtal in sorted(samples_to_export):
            self.Logfile.insert('%s: getting compound code from database' %xtal)
            query=self.db.execute_statement("select CompoundCode from mainTable where CrystalName='%s';" %xtal)
            compoundID=str(query[0][0])
            self.Logfile.insert('%s: compounds code = %s' %(xtal,compoundID))
#            compoundID=str(item[1])
            if os.path.isfile(os.path.join(self.initial_model_directory,xtal,xtal+'.free.mtz')):
                if os.path.isfile(os.path.join(self.initial_model_directory,xtal,xtal+'-pandda-model.pdb')):
                    self.Logfile.insert('running inital refinement on PANDDA model of '+xtal)
                    Serial=XChemRefine.GetSerial(self.initial_model_directory,xtal)
                    #######################################################
                    if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'cootOut')):
                        os.mkdir(os.path.join(self.initial_model_directory,xtal,'cootOut'))
                    # create folder for new refinement cycle
                    if os.path.isdir(os.path.join(self.initial_model_directory,xtal,'cootOut','Refine_'+str(Serial))):
                        os.chdir(os.path.join(self.initial_model_directory,xtal,'cootOut','Refine_'+str(Serial)))
                    else:
                        os.mkdir(os.path.join(self.initial_model_directory,xtal,'cootOut','Refine_'+str(Serial)))
                        os.chdir(os.path.join(self.initial_model_directory,xtal,'cootOut','Refine_'+str(Serial)))
                    os.system('/bin/cp %s in.pdb' %os.path.join(self.initial_model_directory,xtal,xtal+'-pandda-model.pdb'))
                    Refine=XChemRefine.Refine(self.initial_model_directory,xtal,compoundID,self.datasource)
                    Refine.RunBuster(str(Serial),self.external_software,self.xce_logfile,None)
                else:
                    self.Logfile.error('%s: cannot find %s-pandda-model.pdb; cannot start refinement...' %(xtal,xtal))

            elif xtal in samples_to_export and not os.path.isfile(
                    os.path.join(self.initial_model_directory, xtal, xtal + '.free.mtz')):
                self.Logfile.error('%s: cannot start refinement because %s.free.mtz is missing in %s' % (
                xtal, xtal, os.path.join(self.initial_model_directory, xtal)))
            else:
                self.Logfile.insert('%s: nothing to refine' % (xtal))


    def export_models(self):

        self.Logfile.insert('finding out which PanDDA models need to be exported')

        # first find which samples are in interesting datasets and have a model
        # and determine the timestamp
        fileModelsDict={}
        queryModels=''
        for model in glob.glob(os.path.join(self.panddas_directory,'processed_datasets','*','modelled_structures','*-pandda-model.pdb')):
            sample=model[model.rfind('/')+1:].replace('-pandda-model.pdb','')
            timestamp=datetime.fromtimestamp(os.path.getmtime(model)).strftime('%Y-%m-%d %H:%M:%S')
            self.Logfile.insert(sample+'-pandda-model.pdb was created on '+str(timestamp))
            queryModels+="'"+sample+"',"
            fileModelsDict[sample]=timestamp

        # now get these models from the database and compare the datestamps
        # Note: only get the models that underwent some form of refinement,
        #       because only if the model was updated in pandda.inspect will it be exported and refined
        dbModelsDict={}
        if queryModels != '':
            dbEntries=self.db.execute_statement("select CrystalName,DatePanDDAModelCreated from mainTable where CrystalName in ("+queryModels[:-1]+") and (RefinementOutcome like '3%' or RefinementOutcome like '4%' or RefinementOutcome like '5%')")
            for item in dbEntries:
                xtal=str(item[0])
                timestamp=str(item[1])
                dbModelsDict[xtal]=timestamp
                self.Logfile.insert('PanDDA model for '+xtal+' is in database and was created on '+str(timestamp))

        # compare timestamps and only export the ones where the timestamp of the file is newer than the one in the DB
        samples_to_export={}
        self.Logfile.insert('checking which PanDDA models were newly created or updated')
        if self.which_models=='all':
            self.Logfile.insert('Note: you chose to export ALL available PanDDA!')

        for sample in fileModelsDict:
            if self.which_models=='all':
                self.Logfile.insert('exporting '+sample)
                samples_to_export[sample]=fileModelsDict[sample]
            else:
                if sample in dbModelsDict:
                    try:
                        difference=(datetime.strptime(fileModelsDict[sample],'%Y-%m-%d %H:%M:%S') - datetime.strptime(dbModelsDict[sample],'%Y-%m-%d %H:%M:%S')  )
                        if difference.seconds != 0:
                            self.Logfile.insert('exporting '+sample+' -> was already refined, but newer PanDDA model available')
                            samples_to_export[sample]=fileModelsDict[sample]
                    except ValueError:
                        # this will be raised if timestamp is not properly formatted;
                        # which will usually be the case when respective field in database is blank
                        # these are hopefully legacy cases which are from before this extensive check was introduced (13/01/2017)
                        advice = (  'The pandda model of '+xtal+' was changed, but it was already refined! '
                                    'This is most likely because this was done with an older version of XCE. '
                                    'If you really want to export and refine this model, you need to open the database '
                                    'with DBbroweser (sqlitebrowser.org); then change the RefinementOutcome field '
                                    'of the respective sample to "2 - PANDDA model", save the database and repeat the export prodedure.'        )
                        self.Logfile.insert(advice)
                else:
                    self.Logfile.insert('exporting '+sample+' -> first time to be exported and refined')
                    samples_to_export[sample]=fileModelsDict[sample]

        # update the DB:
        # set timestamp to current timestamp of file and set RefinementOutcome to '2-pandda...'

        if samples_to_export != {}:
            select_dir_string=''
            select_dir_string_new_pannda=' '
            for sample in samples_to_export:
                self.Logfile.insert('changing directory to ' + os.path.join(self.initial_model_directory,sample))
                os.chdir(os.path.join(self.initial_model_directory,sample))
                self.Logfile.insert(sample + ': copying ' + os.path.join(self.panddas_directory,'processed_datasets',sample,'modelled_structures',sample+'-pandda-model.pdb'))
                os.system('/bin/cp %s .' %os.path.join(self.panddas_directory,'processed_datasets',sample,'modelled_structures',sample+'-pandda-model.pdb'))
                db_dict= {'RefinementOutcome': '2 - PANDDA model', 'DatePanDDAModelCreated': samples_to_export[sample]}
                for old_event_map in glob.glob('*-BDC_*.ccp4'):
                    if not os.path.isdir('old_event_maps'):
                        os.mkdir('old_event_maps')
                        self.Logfile.warning(sample + ': moving ' + old_event_map + ' to old_event_maps folder')
                        os.system('/bin/mv %s old_event_maps' %old_event_map)
                for event_map in glob.glob(os.path.join(self.panddas_directory,'processed_datasets',sample,'*-BDC_*.ccp4')):
                    self.Logfile.insert(sample + ': copying ' + event_map)
                    os.system('/bin/cp %s .' %event_map)
                select_dir_string+="select_dir={0!s} ".format(sample)
                select_dir_string_new_pannda+='{0!s} '.format(sample)
                self.Logfile.insert('updating database for '+sample+' setting time model was created to '+db_dict['DatePanDDAModelCreated']+' and RefinementOutcome to '+db_dict['RefinementOutcome'])
                self.db.update_data_source(sample,db_dict)


        return samples_to_export
















class run_pandda_export(QtCore.QThread):

    def __init__(self,panddas_directory,datasource,initial_model_directory,xce_logfile,update_datasource_only,which_models):
        QtCore.QThread.__init__(self)
        self.panddas_directory=panddas_directory
        self.datasource=datasource
        self.initial_model_directory=initial_model_directory
        self.db=XChemDB.data_source(self.datasource)
        self.db.create_missing_columns()
        self.db_list=self.db.get_empty_db_dict()
        self.external_software=XChemUtils.external_software(xce_logfile).check()
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.update_datasource_only=update_datasource_only
        self.which_models=which_models
        self.already_exported_models=[]

        self.RefmacParams={ 'HKLIN':            '',                 'HKLOUT': '',
                            'XYZIN':            '',                 'XYZOUT': '',
                            'LIBIN':            '',                 'LIBOUT': '',
                            'TLSIN':            '',                 'TLSOUT': '',
                            'TLSADD':           '',
                            'NCYCLES':          '10',
                            'MATRIX_WEIGHT':    'AUTO',
                            'BREF':             '    bref ISOT\n',
                            'TLS':              '',
                            'NCS':              '',
                            'TWIN':             ''    }

    def run(self):

        # v1.3.8.2 - removed option to update database only
#        if not self.update_datasource_only:
        samples_to_export=self.export_models()

        self.import_samples_into_datasouce(samples_to_export)

#        if not self.update_datasource_only:
        self.refine_exported_models(samples_to_export)


    def refine_exported_models(self,samples_to_export):
        self.Logfile.insert('will try to refine the following crystals:')
        for xtal in samples_to_export: self.Logfile.insert(xtal)
#        sample_list=self.db.execute_statement("select CrystalName,CompoundCode from mainTable where RefinementOutcome='2 - PANDDA model';")
#        for item in sample_list:
#            xtal=str(item[0])
        for xtal in sorted(samples_to_export):
            self.Logfile.insert('%s: getting compound code from database' %xtal)
            query=self.db.execute_statement("select CompoundCode from mainTable where CrystalName='%s';" %xtal)
            compoundID=str(query[0][0])
            self.Logfile.insert('%s: compounds code = %s' %(xtal,compoundID))
#            compoundID=str(item[1])
            if os.path.isfile(os.path.join(self.initial_model_directory,xtal,xtal+'.free.mtz')):
                if os.path.isfile(os.path.join(self.initial_model_directory,xtal,xtal+'-ensemble-model.pdb')):
                    self.Logfile.insert('running inital refinement on PANDDA model of '+xtal)
                    Serial=XChemRefine.GetSerial(self.initial_model_directory,xtal)
                    #######################################################
                    if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'cootOut')):
                        os.mkdir(os.path.join(self.initial_model_directory,xtal,'cootOut'))
                    # create folder for new refinement cycle
                    if os.path.isdir(os.path.join(self.initial_model_directory,xtal,'cootOut','Refine_'+str(Serial))):
                        os.chdir(os.path.join(self.initial_model_directory,xtal,'cootOut','Refine_'+str(Serial)))
                        try:
                            os.system('/bin/rm *-ensemble-model.pdb *restraints*')
                        except:
                            self.Logfile.error("Restraint files didn't exist to remove. Will try to continue")
                    else:
                        os.mkdir(os.path.join(self.initial_model_directory,xtal,'cootOut','Refine_'+str(Serial)))
                        os.chdir(os.path.join(self.initial_model_directory,xtal,'cootOut','Refine_'+str(Serial)))
                    Refine=XChemRefine.panddaRefine(self.initial_model_directory,xtal,compoundID,self.datasource)
                    os.symlink(os.path.join(self.initial_model_directory,xtal,xtal+'-ensemble-model.pdb'),xtal+'-ensemble-model.pdb')
                    Refine.RunQuickRefine(Serial,self.RefmacParams,self.external_software,self.xce_logfile,'pandda_refmac',None)

#            elif xtal in os.path.join(self.panddas_directory,'processed_datasets',xtal,'modelled_structures',
#                                      '{}-pandda-model.pdb'.format(xtal)):
#                self.Logfile.insert('{}: cannot start refinement because {}'.format(xtal,xtal) +
#                                   ' does not have a modelled structure. Check whether you expect this dataset to ' +
#                                   ' have a modelled structure, compare pandda.inspect and datasource,'
#                                   ' then tell XCHEMBB ')
                else:
                    self.Logfile.error('%s: cannot find %s-ensemble-model.pdb; cannot start refinement...' %(xtal,xtal))
                    self.Logfile.error('Please check terminal window for any PanDDA related tracebacks')

            elif xtal in samples_to_export and not os.path.isfile(
                    os.path.join(self.initial_model_directory, xtal, xtal + '.free.mtz')):
                self.Logfile.error('%s: cannot start refinement because %s.free.mtz is missing in %s' % (
                xtal, xtal, os.path.join(self.initial_model_directory, xtal)))
            else:
                self.Logfile.insert('%s: nothing to refine' % (xtal))

    def import_samples_into_datasouce(self,samples_to_export):

        # first make a note of all the datasets which were used in pandda directory
        os.chdir(os.path.join(self.panddas_directory,'processed_datasets'))
        for xtal in glob.glob('*'):
            self.db.execute_statement("update mainTable set DimplePANDDAwasRun = 'True',DimplePANDDAreject = 'False',DimplePANDDApath='{0!s}' where CrystalName is '{1!s}'".format(self.panddas_directory, xtal))
        # do the same as before, but look for rejected datasets

        try:
            os.chdir(os.path.join(self.panddas_directory,'rejected_datasets'))
            for xtal in glob.glob('*'):
                self.db.execute_statement("update mainTable set DimplePANDDAwasRun = 'True',DimplePANDDAreject = 'True',DimplePANDDApath='{0!s}',DimplePANDDAhit = 'False' where CrystalName is '{1!s}'".format(self.panddas_directory, xtal))
        except OSError:
            pass

        site_list = []
        pandda_hit_list=[]

        with open(os.path.join(self.panddas_directory,'analyses','pandda_inspect_sites.csv'),'rb') as csv_import:
            csv_dict = csv.DictReader(csv_import)
            self.Logfile.insert('reding pandda_inspect_sites.csv')
            for i,line in enumerate(csv_dict):
                self.Logfile.insert(str(line).replace('\n','').replace('\r',''))
                site_index=line['site_idx']
                name=line['Name'].replace("'","")
                comment=line['Comment']
                site_list.append([site_index,name,comment])
                self.Logfile.insert('add to site_list_:' + str([site_index,name,comment]))

        progress_step=1
        for i,line in enumerate(open(os.path.join(self.panddas_directory,'analyses','pandda_inspect_events.csv'))):
            n_lines=i
        if n_lines != 0:
            progress_step=100/float(n_lines)
        else:
            progress_step=0
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.Logfile.insert('reading '+os.path.join(self.panddas_directory,'analyses','pandda_inspect_events.csv'))
        with open(os.path.join(self.panddas_directory,'analyses','pandda_inspect_events.csv'),'rb') as csv_import:
            csv_dict = csv.DictReader(csv_import)

            for i,line in enumerate(csv_dict):
                db_dict={}
                sampleID=line['dtag']
                if sampleID not in samples_to_export:
                    self.Logfile.warning('%s: not to be exported; will not add to panddaTable...' %sampleID)
                    continue
                if sampleID not in pandda_hit_list:
                    pandda_hit_list.append(sampleID)
                site_index=str(line['site_idx']).replace('.0','')
                event_index=str(line['event_idx']).replace('.0','')
                self.Logfile.insert(str(line))
                self.Logfile.insert('reading {0!s} -> site {1!s} -> event {2!s}'.format(sampleID, site_index, event_index))

                for entry in site_list:
                    if entry[0]==site_index:
                        site_name=entry[1]
                        site_comment=entry[2]
                        break

                # check if EVENT map exists in project directory
                event_map=''
                for file in glob.glob(os.path.join(self.initial_model_directory,sampleID,'*ccp4')):
                    filename=file[file.rfind('/')+1:]
                    if filename.startswith(sampleID+'-event_'+event_index) and filename.endswith('map.native.ccp4'):
                        event_map=file
                        self.Logfile.insert('found respective event maps in {0!s}: {1!s}'.format(self.initial_model_directory, event_map))
                        break

                # initial pandda model and mtz file
                pandda_model=''
                for file in glob.glob(os.path.join(self.initial_model_directory,sampleID,'*pdb')):
                    filename=file[file.rfind('/')+1:]
                    if filename.endswith('-ensemble-model.pdb'):
                        pandda_model=file
                        if sampleID not in self.already_exported_models:
                            self.already_exported_models.append(sampleID)
                        break
                inital_mtz=''
                for file in glob.glob(os.path.join(self.initial_model_directory,sampleID,'*mtz')):
                    filename=file[file.rfind('/')+1:]
                    if filename.endswith('pandda-input.mtz'):
                        inital_mtz=file
                        break

                db_dict['CrystalName']                      =   sampleID
                db_dict['PANDDApath']                       =   self.panddas_directory
                db_dict['PANDDA_site_index']                =   site_index
                db_dict['PANDDA_site_name']                 =   site_name
                db_dict['PANDDA_site_comment']              =   site_comment
                db_dict['PANDDA_site_event_index']          =   event_index
                db_dict['PANDDA_site_event_comment']        =   line['Comment'].replace("'","")
                db_dict['PANDDA_site_confidence']           =   line['Ligand Confidence']
                db_dict['PANDDA_site_InspectConfidence']    =   line['Ligand Confidence']
                db_dict['PANDDA_site_ligand_placed']        =   line['Ligand Placed']
                db_dict['PANDDA_site_viewed']               =   line['Viewed']
                db_dict['PANDDA_site_interesting']          =   line['Interesting']
                db_dict['PANDDA_site_z_peak']               =   line['z_peak']
                db_dict['PANDDA_site_x']                    =   line['x']
                db_dict['PANDDA_site_y']                    =   line['y']
                db_dict['PANDDA_site_z']                    =   line['z']
                db_dict['PANDDA_site_ligand_id']            =   ''
                db_dict['PANDDA_site_event_map']            =   event_map
                db_dict['PANDDA_site_initial_model']        =   pandda_model
                db_dict['PANDDA_site_initial_mtz']          =   inital_mtz
                db_dict['PANDDA_site_spider_plot']          =   ''

                # find apo structures which were used
                # XXX missing XXX

                self.db.update_insert_site_event_panddaTable(sampleID,db_dict)

                # this is necessary, otherwise RefinementOutcome will be reset for samples that are actually already in refinement
                self.db.execute_statement("update panddaTable set RefinementOutcome = '2 - PANDDA model' where CrystalName is '{0!s}' and RefinementOutcome is null".format(sampleID))
                self.db.execute_statement("update mainTable set RefinementOutcome = '2 - PANDDA model' where CrystalName is '{0!s}' and (RefinementOutcome is null or RefinementOutcome is '1 - Analysis Pending')".format(sampleID))
                self.db.execute_statement("update mainTable set DimplePANDDAhit = 'True' where CrystalName is '{0!s}'".format(sampleID))
                progress += progress_step
                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.Logfile.insert('done reading pandda_inspect_sites.csv')

        # finally find all samples which do not have a pandda hit
        os.chdir(os.path.join(self.panddas_directory,'processed_datasets'))
        self.Logfile.insert('check which datasets are not interesting')
        # DimplePANDDAhit
#        for xtal in glob.glob('*'):
#            if xtal not in pandda_hit_list:
#                self.Logfile.insert(xtal+': not in interesting_datasets; updating database...')
#                self.db.execute_statement("update mainTable set DimplePANDDAhit = 'False' where CrystalName is '{0!s}'".format(xtal))


    def export_models(self):

        self.Logfile.insert('finding out which PanDDA models need to be exported')

        # first find which samples are in interesting datasets and have a model
        # and determine the timestamp
        fileModelsDict={}
        queryModels=''
        for model in glob.glob(os.path.join(self.panddas_directory,'processed_datasets','*','modelled_structures','*-pandda-model.pdb')):
            sample=model[model.rfind('/')+1:].replace('-pandda-model.pdb','')
            timestamp=datetime.fromtimestamp(os.path.getmtime(model)).strftime('%Y-%m-%d %H:%M:%S')
            self.Logfile.insert(sample+'-pandda-model.pdb was created on '+str(timestamp))
            queryModels+="'"+sample+"',"
            fileModelsDict[sample]=timestamp

        # now get these models from the database and compare the datestamps
        # Note: only get the models that underwent some form of refinement,
        #       because only if the model was updated in pandda.inspect will it be exported and refined
        dbModelsDict={}
        if queryModels != '':
            dbEntries=self.db.execute_statement("select CrystalName,DatePanDDAModelCreated from mainTable where CrystalName in ("+queryModels[:-1]+") and (RefinementOutcome like '3%' or RefinementOutcome like '4%' or RefinementOutcome like '5%')")
            for item in dbEntries:
                xtal=str(item[0])
                timestamp=str(item[1])
                dbModelsDict[xtal]=timestamp
                self.Logfile.insert('PanDDA model for '+xtal+' is in database and was created on '+str(timestamp))

        # compare timestamps and only export the ones where the timestamp of the file is newer than the one in the DB
        samples_to_export={}
        self.Logfile.insert('checking which PanDDA models were newly created or updated')
        if self.which_models=='all':
            self.Logfile.insert('Note: you chose to export ALL available PanDDA!')

        for sample in fileModelsDict:
            if self.which_models=='all':
                self.Logfile.insert('exporting '+sample)
                samples_to_export[sample]=fileModelsDict[sample]
            else:
                if sample in dbModelsDict:
                    try:
                        difference=(datetime.strptime(fileModelsDict[sample],'%Y-%m-%d %H:%M:%S') - datetime.strptime(dbModelsDict[sample],'%Y-%m-%d %H:%M:%S')  )
                        if difference.seconds != 0:
                            self.Logfile.insert('exporting '+sample+' -> was already refined, but newer PanDDA model available')
                            samples_to_export[sample]=fileModelsDict[sample]
                    except ValueError:
                        # this will be raised if timestamp is not properly formatted;
                        # which will usually be the case when respective field in database is blank
                        # these are hopefully legacy cases which are from before this extensive check was introduced (13/01/2017)
                        advice = (  'The pandda model of '+xtal+' was changed, but it was already refined! '
                                    'This is most likely because this was done with an older version of XCE. '
                                    'If you really want to export and refine this model, you need to open the database '
                                    'with DBbroweser (sqlitebrowser.org); then change the RefinementOutcome field '
                                    'of the respective sample to "2 - PANDDA model", save the database and repeat the export prodedure.'        )
                        self.Logfile.insert(advice)
                else:
                    self.Logfile.insert('exporting '+sample+' -> first time to be exported and refined')
                    samples_to_export[sample]=fileModelsDict[sample]

        # update the DB:
        # set timestamp to current timestamp of file and set RefinementOutcome to '2-pandda...'

        if samples_to_export != {}:
            select_dir_string=''
            select_dir_string_new_pannda=' '
            for sample in samples_to_export:
                db_dict= {'RefinementOutcome': '2 - PANDDA model', 'DatePanDDAModelCreated': samples_to_export[sample]}
                select_dir_string+="select_dir={0!s} ".format(sample)
                select_dir_string_new_pannda+='{0!s} '.format(sample)
                self.Logfile.insert('updating database for '+sample+' setting time model was created to '+db_dict['DatePanDDAModelCreated']+' and RefinementOutcome to '+db_dict['RefinementOutcome'])
                self.db.update_data_source(sample,db_dict)


            if os.path.isdir(os.path.join(self.panddas_directory,'rejected_datasets')):

                Cmds = (

                    'pandda.export'
                    ' pandda_dir=%s' %self.panddas_directory+
                    ' export_dir={0!s}'.format(self.initial_model_directory)+
                    ' {0!s}'.format(select_dir_string)+
                    ' export_ligands=False'
                    ' generate_occupancy_groupings=True\n'
                    )

            else:

                Cmds = (
                    'source /dls/science/groups/i04-1/software/pandda-update/ccp4/ccp4-7.0/bin/ccp4.setup-sh\n'
#                    'source '+os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh')+'\n'
                    'pandda.export'
                    ' pandda_dir=%s' %self.panddas_directory+
                    ' export_dir={0!s}'.format(self.initial_model_directory)+
                    ' {0!s}'.format(select_dir_string_new_pannda)+
                    ' generate_restraints=True\n'
                    )



            self.Logfile.insert('running pandda.export with the following settings:\n'+Cmds)
            os.system(Cmds)

        return samples_to_export


class run_pandda_analyse(QtCore.QThread):

    def __init__(self,pandda_params,xce_logfile,datasource):
        QtCore.QThread.__init__(self)
        self.data_directory=pandda_params['data_dir']
        self.panddas_directory=pandda_params['out_dir']
        self.submit_mode=pandda_params['submit_mode']

        self.pandda_analyse_data_table = pandda_params['pandda_table']
        self.nproc=pandda_params['nproc']
        self.min_build_datasets=pandda_params['min_build_datasets']
        self.pdb_style=pandda_params['pdb_style']
        self.mtz_style=pandda_params['mtz_style']
        self.sort_event=pandda_params['sort_event']
        self.number_of_datasets=pandda_params['N_datasets']
        self.max_new_datasets=pandda_params['max_new_datasets']
        self.grid_spacing=pandda_params['grid_spacing']
        self.reference_dir=pandda_params['reference_dir']
        self.filter_pdb=os.path.join(self.reference_dir,pandda_params['filter_pdb'])
        self.wilson_scaling = pandda_params['perform_diffraction_data_scaling']
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.datasource=datasource
        self.db=XChemDB.data_source(datasource)
        self.appendix=pandda_params['appendix']
        self.write_mean_maps=pandda_params['write_mean_map']
        self.calc_map_by = pandda_params['average_map']
        self.select_ground_state_model=''
        projectDir = self.data_directory.replace('/*', '')
        self.make_ligand_links='$CCP4/bin/ccp4-python %s %s %s\n' %(os.path.join(os.getenv('XChemExplorer_DIR'),
                                                                                 'helpers',
                                                                                 'make_ligand_links_after_pandda.py')
                                                                    ,projectDir,self.panddas_directory)
        self.use_remote = pandda_params['use_remote']
        self.remote_string = pandda_params['remote_string']

        if self.appendix != '':
            self.panddas_directory=os.path.join(self.reference_dir,'pandda_'+self.appendix)
            if os.path.isdir(self.panddas_directory):
                os.system('/bin/rm -fr %s' %self.panddas_directory)
            os.mkdir(self.panddas_directory)
            self.select_ground_state_model='$CCP4/bin/ccp4-python %s %s\n' %(os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','select_ground_state_dataset.py'),self.panddas_directory)
            self.make_ligand_links=''

    def run(self):

        # print self.reference_dir
        # print self.filter_pdb

        # how to run pandda.analyse on large datasets
        #
        # 1) Run the normal pandda command, with the new setting, e.g.
        # pandda.analyse data_dirs=... max_new_datasets=500
        # This will do the analysis on the first 500 datasets and build the statistical maps - just as normal.
        #
        # 2) Run pandda with the same command:
        # pandda.analyse data_dirs=... max_new_datasets=500
        # This will add 500 new datasets, and process them using the existing statistical maps
        # (this will be quicker than the original analysis). It will then merge the results of the two analyses.
        #
        # 3) Repeat 2) until you don't add any "new" datasets. Then you can build the models as normal.


        number_of_cyles=int(self.number_of_datasets)/int(self.max_new_datasets)
        if int(self.number_of_datasets) % int(self.max_new_datasets) != 0:  # modulo gives remainder after integer division
            number_of_cyles+=1
        self.Logfile.insert('will run %s rounds of pandda.analyse' %str(number_of_cyles))

        if os.path.isfile(os.path.join(self.panddas_directory,'pandda.running')):
            self.Logfile.insert('it looks as if a pandda.analyse job is currently running in: '+self.panddas_directory)
            msg = ( 'there are three possibilities:\n'
                    '1.) choose another PANDDA directory\n'
                    '2.) - check if the job is really running either on the cluster (qstat) or on your local machine\n'
                    '    - if so, be patient and wait until the job has finished\n'
                    '3.) same as 2., but instead of waiting, kill the job and remove at least the pandda.running file\n'
                    '   (or all the contents in the directory if you want to start from scratch)\n'   )
            self.Logfile.insert(msg)
            return None
        else:
#            if os.getenv('SHELL') == '/bin/tcsh' or os.getenv('SHELL') == '/bin/csh':
#                source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-csh\n')
#            elif os.getenv('SHELL') == '/bin/bash' or self.use_remote:
#                source_file='export XChemExplorer_DIR="'+os.getenv('XChemExplorer_DIR')+'"\n'
#                source_file+='source %s\n' %os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh\n')
#            else:
#                source_file=''
            # v1.2.1 - pandda.setup files should be obsolete now that pandda is part of ccp4
            source_file='source /dls/science/groups/i04-1/software/pandda_0.2.12/ccp4/ccp4-7.0/bin/ccp4.setup-sh\n'
            source_file += 'export XChemExplorer_DIR="' + os.getenv('XChemExplorer_DIR') + '"\n'

            if os.path.isfile(self.filter_pdb + '.pdb'):
                print('filter pdb located')
                filter_pdb=' filter.pdb='+self.filter_pdb+'.pdb'
                print('will use ' + filter_pdb + 'as a filter for pandda.analyse')
            else:
                if self.use_remote:
                    stat_command = self.remote_string.replace("qsub'", str('stat ' + self.filter_pdb + "'"))
                    output = subprocess.Popen(stat_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    out, err = output.communicate()
                    print out
                    if 'cannot stat' in out:
                        filter_pdb = ''
                    else:
                        filter_pdb = ' filter.pdb=' + self.filter_pdb + '.pdb'

                else:
                    filter_pdb=''

            os.chdir(self.panddas_directory)

            # note: copied latest pandda.setup-sh from XCE2 installation (08/08/2017)

            dls = ''
            if self.data_directory.startswith('/dls'):
                dls = (
                    source_file +
                    '\n' + 'module load pymol/1.8.2.0' + ' \n'
                )

            Cmds = (
                '#!'+os.getenv('SHELL')+'\n' +
                '\n' +
                dls +
                'cd ' + self.panddas_directory + '\n' +
                '\n'
                )

            ignore = []
            char = []
            zmap = []

            for i in range(0, self.pandda_analyse_data_table.rowCount()):
                ignore_all_checkbox = self.pandda_analyse_data_table.cellWidget(i, 6)
                ignore_characterisation_checkbox = self.pandda_analyse_data_table.cellWidget(i, 7)
                ignore_zmap_checkbox = self.pandda_analyse_data_table.cellWidget(i, 8)

                if ignore_all_checkbox.isChecked():
                    ignore.append(str(self.pandda_analyse_data_table.item(i, 0).text()))
                if ignore_characterisation_checkbox.isChecked():
                    char.append(str(self.pandda_analyse_data_table.item(i, 0).text()))
                if ignore_zmap_checkbox.isChecked():
                    zmap.append(str(self.pandda_analyse_data_table.item(i, 0).text()))

            print ignore

            def append_to_ignore_string(datasets_list, append_string):
                if len(datasets_list)==0:
                    append_string = ''
                for i in range(0, len(datasets_list)):
                    if i < len(datasets_list)-1:
                        append_string += str(datasets_list[i] + ',')
                    else:
                        append_string += str(datasets_list[i] +'"')
                print(append_string)
                return append_string

            ignore_string = 'ignore_datasets="'
            ignore_string = append_to_ignore_string(ignore, ignore_string)

            char_string = 'exclude_from_characterisation="'
            char_string = append_to_ignore_string(char, char_string)

            zmap_string = 'exclude_from_z_map_analysis="'
            zmap_string = append_to_ignore_string(zmap, zmap_string)

            for i in range(number_of_cyles):
                Cmds += (
                    'pandda.analyse '+
                    ' data_dirs="'+self.data_directory.replace('/*','')+'/*"'+
                    ' out_dir="'+self.panddas_directory+'"'
                    ' min_build_datasets='+self.min_build_datasets+
                    ' max_new_datasets='+self.max_new_datasets+
                    ' grid_spacing='+self.grid_spacing+
                    ' cpus='+self.nproc+
                    ' events.order_by='+self.sort_event+
                    filter_pdb+
                    ' pdb_style='+self.pdb_style+
                    ' mtz_style='+self.mtz_style+
                    ' lig_style=/compound/*.cif'+
                    ' apply_b_factor_scaling='+self.wilson_scaling+
                    ' write_average_map='+self.write_mean_maps +
                    ' average_map=' + self.calc_map_by +
                    ' ' +
                    ignore_string +' '+
                    char_string +' '+
                    zmap_string +' '+
                    '\n'
                    )

            Cmds += self.select_ground_state_model
            Cmds += self.make_ligand_links
            Cmds += '\n'

            data_dir_string = self.data_directory.replace('/*', '')

            Cmds += str(
                        'find ' + data_dir_string +
                        '/*/compound -name "*.cif" | while read line; do  echo ${line//"' +
                        data_dir_string + '"/"' + self.panddas_directory +
                        '/processed_datasets/"}| while read line2; do cp $line ${line2//compound/ligand_files} > /dev/null 2>&1; '
                        'done; done;')

            Cmds += '\n'



            Cmds += str(
                        'find ' + data_dir_string +
                        '/*/compound -name "*.pdb" | while read line; do  echo ${line//"' +
                        data_dir_string + '"/"' + self.panddas_directory +
                        '/processed_datasets/"}| while read line2; do cp $line ${line2//compound/ligand_files} > /dev/null 2>&1; '
                        'done; done;')

            self.Logfile.insert('running pandda.analyse with the following command:\n'+Cmds)



            f = open('pandda.sh','w')
            f.write(Cmds)
            f.close()

#            #>>> for testing
#            self.submit_mode='local machine'

            self.Logfile.insert('trying to run pandda.analyse on ' + str(self.submit_mode))

            if self.submit_mode=='local machine':
                self.Logfile.insert('running PANDDA on local machine')
                os.system('chmod +x pandda.sh')
                os.system('./pandda.sh &')
            elif self.use_remote:
                # handles remote submission of pandda.analyse jobs
                submission_string = self.remote_string.replace("qsub'",
                                                               str('cd ' +
                                                                   self.panddas_directory +
                                                                   '; ' +
                                                                   "qsub -P labxchem -q medium.q -N pandda 5 -l exclusive,m_mem_free=100G pandda.sh'"))
                os.system(submission_string)
                self.Logfile.insert(str('running PANDDA remotely, using: ' + submission_string))
            else:
                self.Logfile.insert('running PANDDA on cluster, using qsub...')
                os.system('qsub -P labxchem -q medium.q -N pandda -l exclusive,m_mem_free=100G pandda.sh')

        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))



class giant_cluster_datasets(QtCore.QThread):

    def __init__(self,initial_model_directory,pandda_params,xce_logfile,datasource,):
        QtCore.QThread.__init__(self)
        self.panddas_directory=pandda_params['out_dir']
        self.pdb_style=pandda_params['pdb_style']
        self.mtz_style=pandda_params['mtz_style']
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.initial_model_directory=initial_model_directory
        self.db=XChemDB.data_source(datasource)


    def run(self):

        self.emit(QtCore.SIGNAL('update_progress_bar'), 0)

        if self.pdb_style.replace(' ','') == '':
            self.Logfile.insert('PDB style is not set in pandda.analyse!')
            self.Logfile.insert('cannot start pandda.analyse')
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'PDB style is not set in pandda.analyse!')
            return None

        if self.mtz_style.replace(' ','') == '':
            self.Logfile.insert('MTZ style is not set in pandda.analyse!')
            self.Logfile.insert('cannot start pandda.analyse')
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'MTZ style is not set in pandda.analyse!')
            return None

        # 1.) prepare output directory
        os.chdir(self.panddas_directory)
        if os.path.isdir('cluster_analysis'):
            self.Logfile.insert('removing old cluster_analysis directory in {0!s}'.format(self.panddas_directory))
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'removing old cluster_analysis directory in {0!s}'.format(self.panddas_directory))
            os.system('/bin/rm -fr cluster_analysis 2> /dev/null')
        self.Logfile.insert('creating cluster_analysis directory in {0!s}'.format(self.panddas_directory))
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating cluster_analysis directory in {0!s}'.format(self.panddas_directory))
        os.mkdir('cluster_analysis')
        self.emit(QtCore.SIGNAL('update_progress_bar'), 10)

        # 2.) go through project directory and make sure that all pdb files really exist
        # broken links derail the giant.cluster_mtzs_and_pdbs script
        self.Logfile.insert('cleaning up broken links of {0!s} and {1!s} in {2!s}'.format(self.pdb_style, self.mtz_style, self.initial_model_directory))
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'cleaning up broken links of {0!s} and {1!s} in {2!s}'.format(self.pdb_style, self.mtz_style, self.initial_model_directory))
        os.chdir(self.initial_model_directory)
        for xtal in glob.glob('*'):
            if not os.path.isfile(os.path.join(xtal,self.pdb_style)):
                self.Logfile.insert('missing {0!s} and {1!s} for {2!s}'.format(self.pdb_style, self.mtz_style, xtal))
                os.system('/bin/rm {0!s}/{1!s} 2> /dev/null'.format(xtal, self.pdb_style))
                os.system('/bin/rm {0!s}/{1!s} 2> /dev/null'.format(xtal, self.mtz_style))
        self.emit(QtCore.SIGNAL('update_progress_bar'), 20)

        # 3.) giant.cluster_mtzs_and_pdbs
        self.Logfile.insert("running giant.cluster_mtzs_and_pdbs {0!s}/*/{1!s} pdb_regex='{2!s}/(.*)/{3!s}' out_dir='{4!s}/cluster_analysis'".format(self.initial_model_directory, self.pdb_style, self.initial_model_directory, self.pdb_style, self.panddas_directory))
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'running giant.cluster_mtzs_and_pdbs')

        if os.getenv('SHELL') == '/bin/tcsh' or os.getenv('SHELL') == '/bin/csh':
            source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-csh')
        elif os.getenv('SHELL') == '/bin/bash':
            source_file=os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh')
        else:
            source_file=''

        Cmds = (
                '#!'+os.getenv('SHELL')+'\n'
                'unset PYTHONPATH\n'
                'source '+source_file+'\n'
                "giant.datasets.cluster %s/*/%s pdb_regex='%s/(.*)/%s' out_dir='%s/cluster_analysis'" %(self.initial_model_directory,self.pdb_style,self.initial_model_directory,self.pdb_style,self.panddas_directory)
            )

#        os.system("giant.cluster_mtzs_and_pdbs %s/*/%s pdb_regex='%s/(.*)/%s' out_dir='%s/cluster_analysis'" %(self.initial_model_directory,self.pdb_style,self.initial_model_directory,self.pdb_style,self.panddas_directory))
        os.system(Cmds)
        self.emit(QtCore.SIGNAL('update_progress_bar'), 80)

        # 4.) analyse output
        self.Logfile.insert('parsing {0!s}/cluster_analysis'.format(self.panddas_directory))
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'parsing {0!s}/cluster_analysis'.format(self.panddas_directory))
        os.chdir('{0!s}/cluster_analysis'.format(self.panddas_directory))
        cluster_dict={}
        for out_dir in sorted(glob.glob('*')):
            if os.path.isdir(out_dir):
                cluster_dict[out_dir]=[]
                for folder in glob.glob(os.path.join(out_dir,'pdbs','*')):
                    xtal=folder[folder.rfind('/')+1:]
                    cluster_dict[out_dir].append(xtal)
        self.emit(QtCore.SIGNAL('update_progress_bar'), 90)

        # 5.) update datasource
        self.Logfile.insert('updating datasource with results from giant.cluster_mtzs_and_pdbs')
        if cluster_dict != {}:
            for key in cluster_dict:
                for xtal in cluster_dict[key]:
                    db_dict= {'CrystalFormName': key}
                    self.db.update_data_source(xtal,db_dict)

        # 6.) finish
        self.emit(QtCore.SIGNAL('update_progress_bar'), 100)
        self.Logfile.insert('finished giant.cluster_mtzs_and_pdbs')
        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))

class check_if_pandda_can_run:

    # reasons why pandda cannot be run
    # - there is currently a job running in the pandda directory
    # - min datasets available is too low
    # - required input paramters are not complete
    # - map amplitude and phase labels don't exist

    def __init__(self,pandda_params,xce_logfile,datasource):
        self.data_directory=pandda_params['data_dir']
        self.panddas_directory=pandda_params['out_dir']
        self.min_build_datasets=pandda_params['min_build_datasets']
        self.pdb_style=pandda_params['pdb_style']
        self.mtz_style=pandda_params['mtz_style']
        self.input_dir_structure=pandda_params['pandda_dir_structure']
        self.problem_found=False
        self.error_code=-1
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.db=XChemDB.data_source(datasource)

    def number_of_available_datasets(self):
        counter=0
        for file in glob.glob(os.path.join(self.input_dir_structure,self.pdb_style)):
            if os.path.isfile(file):
                counter+=1
        self.Logfile.insert('pandda.analyse: found {0!s} useable datasets'.format(counter))
        return counter

    def get_first_dataset_in_project_directory(self):
        first_dataset=''
        for file in glob.glob(os.path.join(self.input_dir_structure,self.pdb_style)):
            if os.path.isfile(file):
                first_dataset=file
                break
        return first_dataset

    def compare_number_of_atoms_in_reference_vs_all_datasets(self,refData,dataset_list):
        mismatched_datasets=[]
        pdbtools=XChemUtils.pdbtools(refData)
        refPDB=refData[refData.rfind('/')+1:]
        refPDBlist=pdbtools.get_init_pdb_as_list()
        n_atom_ref=len(refPDBlist)
        for n_datasets,dataset in enumerate(dataset_list):
            if os.path.isfile(os.path.join(self.data_directory.replace('*',''),dataset,self.pdb_style)):
                n_atom=len(pdbtools.get_pdb_as_list(os.path.join(self.data_directory.replace('*',''),dataset,self.pdb_style)))
                if n_atom_ref == n_atom:
                    self.Logfile.insert('{0!s}: atoms in PDB file ({1!s}): {2!s}; atoms in Reference file: {3!s} ===> OK'.format(dataset, self.pdb_style, str(n_atom), str(n_atom_ref)))
                if n_atom_ref != n_atom:
                    self.Logfile.insert('{0!s}: atoms in PDB file ({1!s}): {2!s}; atoms in Reference file: {3!s} ===> ERROR'.format(dataset, self.pdb_style, str(n_atom), str(n_atom_ref)))
                    mismatched_datasets.append(dataset)
        return n_datasets,mismatched_datasets

    def get_datasets_which_fit_to_reference_file(self,ref,reference_directory,cluster_dict,allowed_unitcell_difference_percent):
        refStructure=XChemUtils.pdbtools(os.path.join(reference_directory,ref+'.pdb'))
        symmRef=refStructure.get_spg_number_from_pdb()
        ucVolRef=refStructure.calc_unitcell_volume_from_pdb()
        cluster_dict[ref]=[]
        cluster_dict[ref].append(os.path.join(reference_directory,ref+'.pdb'))
        for dataset in glob.glob(os.path.join(self.data_directory,self.pdb_style)):
            datasetStructure=XChemUtils.pdbtools(dataset)
            symmDataset=datasetStructure.get_spg_number_from_pdb()
            ucVolDataset=datasetStructure.calc_unitcell_volume_from_pdb()

            if symmDataset == symmRef:
                try:
                    difference=math.fabs(1-(float(ucVolRef)/float(ucVolDataset)))*100
                    if difference < allowed_unitcell_difference_percent:
                        sampleID=dataset.replace('/'+self.pdb_style,'')[dataset.replace('/'+self.pdb_style,'').rfind('/')+1:]
                        cluster_dict[ref].append(sampleID)
                except ZeroDivisionError:
                    continue
        return cluster_dict

    def remove_dimple_files(self,dataset_list):
        for n_datasets,dataset in enumerate(dataset_list):
            db_dict={}
            if os.path.isfile(os.path.join(self.data_directory.replace('*',''),dataset,self.pdb_style)):
                os.system('/bin/rm '+os.path.join(self.data_directory.replace('*',''),dataset,self.pdb_style))
                self.Logfile.insert('{0!s}: removing {1!s}'.format(dataset, self.pdb_style))
                db_dict['DimplePathToPDB']=''
                db_dict['DimpleRcryst']=''
                db_dict['DimpleRfree']=''
                db_dict['DimpleResolutionHigh']=''
                db_dict['DimpleStatus']='pending'
            if os.path.isfile(os.path.join(self.data_directory.replace('*',''),dataset,self.mtz_style)):
                os.system('/bin/rm '+os.path.join(self.data_directory.replace('*',''),dataset,self.mtz_style))
                self.Logfile.insert('{0!s}: removing {1!s}'.format(dataset, self.mtz_style))
                db_dict['DimplePathToMTZ']=''
            if db_dict != {}:
                self.db.update_data_source(dataset,db_dict)


    def analyse_pdb_style(self):
        pdb_found=False
        for file in glob.glob(os.path.join(self.data_directory,self.pdb_style)):
            if os.path.isfile(file):
                pdb_found=True
                break
        if not pdb_found:
            self.error_code=1
            message=self.warning_messages()
        return message

    def analyse_mtz_style(self):
        mtz_found=False
        for file in glob.glob(os.path.join(self.data_directory,self.mtz_style)):
            if os.path.isfile(file):
                mtz_found=True
                break
        if not mtz_found:
            self.error_code=2
            message=self.warning_messages()
        return message

    def analyse_min_build_dataset(self):
        counter=0
        for file in glob.glob(os.path.join(self.data_directory,self.mtz_style)):
            if os.path.isfile(file):
                counter+=1
        if counter <= self.min_build_datasets:
            self.error_code=3
            message=self.warning_messages()
        return message


    def warning_messages(self):
        message=''
        if self.error_code==1:
            message='PDB file does not exist'
        if self.error_code==2:
            message='MTZ file does not exist'
        if self.error_code==3:
            message='Not enough datasets available'

        return message


class convert_all_event_maps_in_database(QtCore.QThread):

    def __init__(self,initial_model_directory,xce_logfile,datasource):
        QtCore.QThread.__init__(self)
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.initial_model_directory=initial_model_directory
        self.datasource=datasource
        self.db=XChemDB.data_source(datasource)

    def run(self):
        sqlite = (
            'select'
            ' CrystalName,'
            ' PANDDA_site_event_map,'
            ' PANDDA_site_ligand_resname,'
            ' PANDDA_site_ligand_chain,'
            ' PANDDA_site_ligand_sequence_number,'
            ' PANDDA_site_ligand_altLoc '
            'from panddaTable '
            'where PANDDA_site_event_map not like "event%"'
        )
        print sqlite
        query=self.db.execute_statement(sqlite)
        print query
        progress_step=1
        if len(query) != 0:
            progress_step=100/float(len(query))
        else:
            progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        for item in query:
            print item
            xtalID=str(item[0])
            event_map=str(item[1])
            resname=str(item[2])
            chainID=str(item[3])
            resseq=str(item[4])
            altLoc=str(item[5])

            if os.path.isfile(os.path.join(self.initial_model_directory,xtalID,'refine.pdb')):
                os.chdir(os.path.join(self.initial_model_directory,xtalID))
                self.Logfile.insert('extracting ligand ({0!s},{1!s},{2!s},{3!s}) from refine.pdb'.format(str(resname), str(chainID), str(resseq), str(altLoc)))
                XChemUtils.pdbtools(os.path.join(self.initial_model_directory,xtalID,'refine.pdb')).save_specific_ligands_to_pdb(resname,chainID,resseq,altLoc)
                if os.path.isfile('ligand_{0!s}_{1!s}_{2!s}_{3!s}.pdb'.format(str(resname), str(chainID), str(resseq), str(altLoc))):
                    ligand_pdb='ligand_{0!s}_{1!s}_{2!s}_{3!s}.pdb'.format(str(resname), str(chainID), str(resseq), str(altLoc))
                    print os.path.join(self.initial_model_directory,xtalID,ligand_pdb)
                else:
                    self.Logfile.insert('could not extract ligand; trying next...')
                    continue
            else:
                self.Logfile.insert('directory: '+os.path.join(self.initial_model_directory,xtalID)+' -> cannot find refine.pdb; trying next')
                continue

            if os.path.isfile(os.path.join(self.initial_model_directory,xtalID,'refine.mtz')):
                resolution=XChemUtils.mtztools(os.path.join(self.initial_model_directory,xtalID,'refine.mtz')).get_high_resolution_from_mtz()
            else:
                self.Logfile.insert('directory: '+os.path.join(self.initial_model_directory,xtalID)+' -> cannot find refine.mtz; trying next')
                continue

            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'eventMap -> SF for '+event_map)
            convert_event_map_to_SF(self.initial_model_directory,xtalID,event_map,ligand_pdb,self.xce_logfile,self.datasource,resolution).run()

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)



class convert_event_map_to_SF:
    def __init__(self,project_directory,xtalID,event_map,ligand_pdb,xce_logfile,db_file,resolution):
        self.Logfile=XChemLog.updateLog(xce_logfile)

        self.event_map=event_map
        if not os.path.isfile(self.event_map):
            self.Logfile.insert('cannot find Event map: '+self.event_map)
            self.Logfile.insert('cannot convert event_map to structure factors!')
            return None

        self.project_directory=project_directory
        self.xtalID=xtalID
        self.event_map=event_map
        self.ligand_pdb=ligand_pdb
        self.event=event_map[event_map.rfind('/')+1:].replace('.map','').replace('.ccp4','')
        self.db=XChemDB.data_source(db_file)
        self.resolution=resolution

    def run(self):
        os.chdir(os.path.join(self.project_directory,self.xtalID))


        # remove exisiting mtz file
        if os.path.isfile(self.event+'.mtz'):
            self.Logfile.insert('removing existing '+self.event+'.mtz')
            os.system('/bin/rm '+self.event+'.mtz')

        # event maps generated with pandda v0.2 or higher have the same symmetry as the crystal
        # but phenix.map_to_structure_facors only accepts maps in spg P1
        # therefore map is first expanded to full unit cell and spg of map then set tp p1
        # other conversion option like cinvfft give for whatever reason  uninterpretable maps
        self.convert_map_to_p1()

        # run phenix.map_to_structure_factors
        self.run_phenix_map_to_structure_factors()

        self.remove_and_rename_column_labels()

        # check if output files exist
        if not os.path.isfile('{0!s}.mtz'.format(self.event)):
            self.Logfile.insert('cannot find {0!s}.mtz'.format(self.event))
        else:
            self.Logfile.insert('conversion successful, {0!s}.mtz exists'.format(self.event))
            # update datasource with event_map_mtz information
            self.update_database()



    def calculate_electron_density_map(self,mtzin):
        missing_columns=False
        column_dict=XChemUtils.mtztools(mtzin).get_all_columns_as_dict()
        if 'FWT' in column_dict['F'] and 'PHWT' in column_dict['PHS']:
            labin=' labin F1=FWT PHI=PHWT\n'
        elif '2FOFCWT' in column_dict['F'] and 'PH2FOFCWT' in column_dict['PHS']:
            labin=' labin F1=2FOFCWT PHI=PH2FOFCWT\n'
        else:
            missing_columns=True

        if not missing_columns:
            os.chdir(os.path.join(self.project_directory,self.xtalID))
            cmd = (
                'fft hklin '+mtzin+' mapout 2fofc.map << EOF\n'
                +labin+
                'EOF\n'
                    )
            self.Logfile.insert('calculating 2fofc map from '+mtzin)
            os.system(cmd)
        else:
            self.Logfile.insert('cannot calculate 2fofc.map; missing map coefficients')

    def prepare_conversion_script(self):

        os.chdir(os.path.join(self.project_directory, self.xtalID))

        # see also:
        # http://www.phaser.cimr.cam.ac.uk/index.php/Using_Electron_Density_as_a_Model

        if os.getcwd().startswith('/dls'):
            phenix_module='module_load_phenix\n'
        else:
            phenix_module=''

        cmd = (
            '#!'+os.getenv('SHELL')+'\n'
            '\n'
            +phenix_module+
            '\n'
            'pdbset XYZIN %s XYZOUT mask_ligand.pdb << eof\n' %self.ligand_pdb+
            ' SPACEGROUP {0!s}\n'.format(self.space_group)+
            ' CELL {0!s}\n'.format((' '.join(self.unit_cell)))+
            ' END\n'
            'eof\n'
            '\n'
            'ncsmask XYZIN mask_ligand.pdb MSKOUT mask_ligand.msk << eof\n'
            ' GRID %s\n' %(' '.join(self.gridElectronDensityMap))+
            ' RADIUS 10\n'
            ' PEAK 1\n'
            'eof\n'
            '\n'
            'mapmask MAPIN %s MAPOUT onecell_event_map.map << eof\n' %self.event_map+
            ' XYZLIM CELL\n'
            'eof\n'
            '\n'
            'maprot MAPIN onecell_event_map.map MSKIN mask_ligand.msk WRKOUT masked_event_map.map << eof\n'
            ' MODE FROM\n'
            ' SYMMETRY WORK %s\n' %self.space_group_numberElectronDensityMap+
            ' AVERAGE\n'
            ' ROTATE EULER 0 0 0\n'
            ' TRANSLATE 0 0 0\n'
            'eof\n'
            '\n'
            'mapmask MAPIN masked_event_map.map MAPOUT masked_event_map_fullcell.map << eof\n'
            ' XYZLIM CELL\n'
            ' PAD 0.0\n'
            'eof\n'
            '\n'
            'sfall HKLOUT %s.mtz MAPIN masked_event_map_fullcell.map << eof\n' %self.event+
            ' LABOUT FC=FC_event PHIC=PHIC_event\n'
            ' MODE SFCALC MAPIN\n'
            ' RESOLUTION %s\n' %self.resolution+
            ' END\n'
        )

        self.Logfile.insert('preparing script for conversion of Event map to SF')
        f = open('eventMap2sf.sh','w')
        f.write(cmd)
        f.close()
        os.system('chmod +x eventMap2sf.sh')


    def run_conversion_script(self):
        self.Logfile.insert('running conversion script...')
        os.system('./eventMap2sf.sh')

    def convert_map_to_p1(self):
        self.Logfile.insert('running mapmask -> converting map to p1...')
        cmd = (     '#!'+os.getenv('SHELL')+'\n'
                    '\n'
                    'mapmask mapin %s mapout %s_p1.map << eof\n' %(self.event_map,self.event) +
                    'xyzlin cell\n'
                    'symmetry p1\n' )
        self.Logfile.insert('mapmask command:\n%s' %cmd)
        os.system(cmd)

    def run_phenix_map_to_structure_factors(self):
        if float(self.resolution) < 1.21:   # program complains if resolution is 1.2 or higher
            self.resolution='1.21'
        self.Logfile.insert('running phenix.map_to_structure_factors {0!s}_p1.map d_min={1!s} output_file_name={2!s}_tmp.mtz'.format(self.event, self.resolution, self.event))
        os.system('phenix.map_to_structure_factors {0!s}_p1.map d_min={1!s} output_file_name={2!s}_tmp.mtz'.format(self.event, self.resolution, self.event))

    def run_cinvfft(self,mtzin):
        # mtzin is usually refine.mtz
        self.Logfile.insert('running cinvfft -mapin {0!s} -mtzin {1!s} -mtzout {2!s}_tmp.mtz -colout event'.format(self.event_map, mtzin, self.event))
        os.system('cinvfft -mapin {0!s} -mtzin {1!s} -mtzout {2!s}_tmp.mtz -colout event'.format(self.event_map, mtzin, self.event))


    def remove_and_rename_column_labels(self):

        cmd = (     '#!'+os.getenv('SHELL')+'\n'
                    '\n'
                    'cad hklin1 %s_tmp.mtz hklout %s.mtz << eof\n' %(self.event,self.event)+
                    ' labin file_number 1 E1=F-obs E2=PHIF\n'
                    ' labout file_number 1 E1=F_ampl E2=PHIF\n'
                    'eof\n'
                    '\n' )
        self.Logfile.insert('running CAD: new column labels F_ampl,PHIF')
        os.system(cmd)

    def remove_and_rename_column_labels_after_cinvfft(self):

        cmd = (     '#!'+os.getenv('SHELL')+'\n'
                    '\n'
                    'cad hklin1 %s_tmp.mtz hklout %s.mtz << eof\n' %(self.event,self.event)+
                    ' labin file_number 1 E1=event.F_phi.F E2=event.F_phi.phi\n'
                    ' labout file_number 1 E1=F_ampl E2=PHIF\n'
                    'eof\n'
                    '\n' )
        self.Logfile.insert('running CAD: renaming event.F_phi.F -> F_ampl and event.F_phi.phi -> PHIF')
        os.system(cmd)



    def update_database(self):
        sqlite = ( "update panddaTable set "
                   " PANDDA_site_event_map_mtz = '%s' " %os.path.join(self.project_directory,self.xtalID,self.event+'.mtz')+
                   " where PANDDA_site_event_map is '{0!s}' ".format(self.event_map)
                    )
        self.db.execute_statement(sqlite)
        self.Logfile.insert('updating data source: '+sqlite)

    def clean_output_directory(self):
        os.system('/bin/rm mask_targetcell.pdb')
        os.system('/bin/rm mask_targetcell.msk')
        os.system('/bin/rm onecell.map')
        os.system('/bin/rm masked_targetcell.map')
        os.system('/bin/rm masked_fullcell.map')
        os.system('/bin/rm eventMap2sf.sh')
        os.system('/bin/rm '+self.ligand_pdb)


class run_pandda_inspect_at_home(QtCore.QThread):

    def __init__(self,panddaDir,xce_logfile):
        QtCore.QThread.__init__(self)
        self.panddaDir=panddaDir
        self.Logfile=XChemLog.updateLog(xce_logfile)

    def run(self):
        os.chdir(os.path.join(self.panddaDir,'processed_datasets'))

        progress_step=1
        if len(glob.glob('*')) != 0:
            progress_step=100/float(len(glob.glob('*')))
        else:
            progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)


        self.Logfile.insert('parsing '+self.panddaDir)
        for xtal in sorted(glob.glob('*')):
            for files in glob.glob(xtal+'/ligand_files/*'):
                if os.path.islink(files):
                    self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'replacing symlink for {0!s} with real file'.format(files))
                    self.Logfile.insert('replacing symlink for {0!s} with real file'.format(files))
                    os.system('cp --remove-destination {0!s} {1!s}/ligand_files'.format(os.path.realpath(files), xtal))
            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        XChemToolTips.run_pandda_inspect_at_home(self.panddaDir)


class convert_apo_structures_to_mmcif(QtCore.QThread):

    def __init__(self,panddaDir,xce_logfile):
        QtCore.QThread.__init__(self)
        self.panddaDir=panddaDir
        self.Logfile=XChemLog.updateLog(xce_logfile)

    def sf_convert_environment(self):
        pdb_extract_init = ''
        if os.path.isdir('/dls'):
            pdb_extract_init = 'source /dls/science/groups/i04-1/software/pdb-extract-prod/setup.sh\n'
            pdb_extract_init += '/dls/science/groups/i04-1/software/pdb-extract-prod/bin/sf_convert'
        else:
            pdb_extract_init = 'source ' + os.path.join(os.getenv('XChemExplorer_DIR'),
                                                            'pdb_extract/pdb-extract-prod/setup.sh') + '\n'
            pdb_extract_init += +os.path.join(os.getenv('XChemExplorer_DIR'),
                                                              'pdb_extract/pdb-extract-prod/bin/sf_convert')
        return pdb_extract_init


    def run(self):
        self.Logfile.insert('converting apo structures in pandda directory to mmcif files')
        self.Logfile.insert('chanfing to '+self.panddaDir)
        progress_step=1
        if len(glob.glob('*')) != 0:
            progress_step=100/float(len(glob.glob(os.path.join(self.panddaDir,'processed_datasets','*'))))
        else:
            progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        pdb_extract_init = self.sf_convert_environment()

        self.Logfile.insert('parsing '+self.panddaDir)
        for dirs in glob.glob(os.path.join(self.panddaDir,'processed_datasets','*')):
            xtal = dirs[dirs.rfind('/')+1:]
            self.Logfile.insert('%s: converting %s to mmcif' %(xtal,xtal+'-pandda-input.mtz'))
            if os.path.isfile(os.path.join(dirs,xtal+'-pandda-input.mtz')):
                if os.path.isfile(os.path.join(dirs,xtal+'_sf.mmcif')):
                    self.Logfile.insert('%s: %s_sf.mmcif exists; skipping...' %(xtal,xtal))
                else:
                    os.chdir(dirs)
                    Cmd = (pdb_extract_init +
                       ' -o mmcif'
                       ' -sf %s' % xtal+'-pandda-input.mtz' +
                       ' -out {0!s}_sf.mmcif  > {1!s}.sf_mmcif.log'.format(xtal, xtal))
                    self.Logfile.insert('running command: '+Cmd)
                    os.system(Cmd)
            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)










class check_number_of_modelled_ligands(QtCore.QThread):

    def __init__(self,project_directory,xce_logfile,db_file):
        QtCore.QThread.__init__(self)
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.project_directory=project_directory
        self.db=XChemDB.data_source(db_file)
        self.errorDict={}

    def update_errorDict(self,xtal,message):
        if xtal not in self.errorDict:
            self.errorDict[xtal]=[]
        self.errorDict[xtal].append(message)

    def insert_new_row_in_panddaTable(self,xtal,ligand,site,dbDict):
        resname=    site[0]
        chain=      site[1]
        seqnum=     site[2]
        altLoc=     site[3]
        x_site=     site[5][0]
        y_site=     site[5][1]
        z_site=     site[5][2]

        resnameSimilarSite= ligand[0]
        chainSimilarSite=   ligand[1]
        seqnumSimilarSite=  ligand[2]

        siteList=[]
        for entry in dbDict[xtal]:
            siteList.append(str(entry[0]))
            if entry[4] == resnameSimilarSite and entry[5] == chainSimilarSite and entry[6] == seqnumSimilarSite:
                eventMap=       str(entry[7])
                eventMap_mtz=   str(entry[8])
                initialPDB=     str(entry[9])
                initialMTZ=     str(entry[10])
                event_id=       str(entry[12])
                PanDDApath=     str(entry[13])

        db_dict={
            'PANDDA_site_index':                    str(int(max(siteList))+1),
            'PANDDApath':                           PanDDApath,
            'PANDDA_site_ligand_id':                resname+'-'+chain+'-'+seqnum,
            'PANDDA_site_ligand_resname':           resname,
            'PANDDA_site_ligand_chain':             chain,
            'PANDDA_site_ligand_sequence_number':   seqnum,
            'PANDDA_site_ligand_altLoc':            'D',
            'PANDDA_site_event_index':              event_id,
            'PANDDA_site_event_map':                eventMap,
            'PANDDA_site_event_map_mtz':            eventMap_mtz,
            'PANDDA_site_initial_model':            initialPDB,
            'PANDDA_site_initial_mtz':              initialMTZ,
            'PANDDA_site_ligand_placed':            'True',
            'PANDDA_site_x':                        x_site,
            'PANDDA_site_y':                        y_site,
            'PANDDA_site_z':                        z_site          }

        print xtal,db_dict





    def run(self):

        self.Logfile.insert('reading modelled ligands from panddaTable')
        dbDict={}

        sqlite = (  "select "
                    " CrystalName,"
                    " PANDDA_site_index,"
                    " PANDDA_site_x,"
                    " PANDDA_site_y,"
                    " PANDDA_site_z,"
                    " PANDDA_site_ligand_resname,"
                    " PANDDA_site_ligand_chain,"
                    " PANDDA_site_ligand_sequence_number,"
                    " PANDDA_site_event_map,"
                    " PANDDA_site_event_map_mtz,"
                    " PANDDA_site_initial_model,"
                    " PANDDA_site_initial_mtz,"
                    " RefinementOutcome,"
                    " PANDDA_site_event_index,"
                    " PANDDApath "
                    "from panddaTable "    )


        dbEntries=self.db.execute_statement(sqlite)
        for item in dbEntries:
            xtal=           str(item[0])
            site=           str(item[1])
            x=              str(item[2])
            y=              str(item[3])
            z=              str(item[4])
            resname=        str(item[5])
            chain=          str(item[6])
            seqnum=         str(item[7])
            eventMap=       str(item[8])
            eventMap_mtz=   str(item[9])
            initialPDB=     str(item[10])
            initialMTZ=     str(item[11])
            outcome=        str(item[12])
            event=          str(item[13])
            PanDDApath=     str(item[14])

            if xtal not in dbDict:
                dbDict[xtal]=[]
            dbDict[xtal].append([site,x,y,z,resname,chain,seqnum,eventMap,eventMap_mtz,initialPDB,initialMTZ,outcome,event,PanDDApath])


        os.chdir(self.project_directory)

        progress_step=1
        if len(glob.glob('*')) != 0:
            progress_step=100/float(len(glob.glob('*')))
        else:
            progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        for xtal in sorted(glob.glob('*')):
            if os.path.isfile(os.path.join(xtal,'refine.pdb')):
                ligands=XChemUtils.pdbtools(os.path.join(xtal,'refine.pdb')).ligand_details_as_list()
                self.Logfile.insert('{0!s}: found file refine.pdb'.format(xtal))
                if ligands:
                    if os.path.isdir(os.path.join(xtal,'xceTmp')):
                        os.system('/bin/rm -fr {0!s}'.format(os.path.join(xtal,'xceTmp')))
                    os.mkdir(os.path.join(xtal,'xceTmp'))
                else:
                    self.Logfile.warning('{0!s}: cannot find ligand molecule in refine.pdb; skipping...'.format(xtal))
                    continue

                made_sym_copies=False
                ligands_not_in_panddaTable=[]
                for n,item in enumerate(ligands):
                    resnameLIG=     item[0]
                    chainLIG=       item[1]
                    seqnumLIG=      item[2]
                    altLocLIG=      item[3]
                    occupancyLig=   item[4]
                    if altLocLIG.replace(' ','') == '':
                        self.Logfile.insert(xtal+': found a ligand not modelled with pandda.inspect -> {0!s} {1!s} {2!s}'.format(resnameLIG, chainLIG, seqnumLIG))
                    residue_xyz = XChemUtils.pdbtools(os.path.join(xtal,'refine.pdb')).get_center_of_gravity_of_residue_ish(item[1],item[2])
                    ligands[n].append(residue_xyz)
                    foundLigand=False
                    if xtal in dbDict:
                        for entry in dbDict[xtal]:
                            resnameTable=entry[4]
                            chainTable=entry[5]
                            seqnumTable=entry[6]
                            self.Logfile.insert('panddaTable: {0!s} {1!s} {2!s} {3!s}'.format(xtal, resnameTable, chainTable, seqnumTable))
                            if resnameLIG == resnameTable and chainLIG == chainTable and seqnumLIG == seqnumTable:
                                self.Logfile.insert('{0!s}: found ligand in database -> {1!s} {2!s} {3!s}'.format(xtal, resnameTable, chainTable, seqnumTable))
                                foundLigand=True
                        if not foundLigand:
                            self.Logfile.error('{0!s}: did NOT find ligand in database -> {1!s} {2!s} {3!s}'.format(xtal, resnameLIG, chainLIG, seqnumLIG))
                            ligands_not_in_panddaTable.append([resnameLIG,chainLIG,seqnumLIG,altLocLIG,occupancyLig,residue_xyz])
                    else:
                        self.Logfile.warning('ligand in PDB file, but dataset not listed in panddaTable: {0!s} -> {1!s} {2!s} {3!s}'.format(xtal, item[0], item[1], item[2]))

                for entry in ligands_not_in_panddaTable:
                    self.Logfile.error('{0!s}: refine.pdb contains a ligand that is not assigned in the panddaTable: {1!s} {2!s} {3!s} {4!s}'.format(xtal, entry[0], entry[1], entry[2], entry[3]))

                for site in ligands_not_in_panddaTable:

                    for files in glob.glob(os.path.join(self.project_directory,xtal,'xceTmp','ligand_*_*.pdb')):
                        mol_xyz = XChemUtils.pdbtools(files).get_center_of_gravity_of_molecule_ish()
                        # now need to check if there is a unassigned entry in panddaTable that is close
                        for entry in dbDict[xtal]:
                            distance = XChemUtils.misc().calculate_distance_between_coordinates(mol_xyz[0], mol_xyz[1],mol_xyz[2],entry[1],entry[2], entry[3])
                            self.Logfile.insert('{0!s}: {1!s} {2!s} {3!s} <---> {4!s} {5!s} {6!s}'.format(xtal, mol_xyz[0], mol_xyz[1], mol_xyz[2], entry[1], entry[2], entry[3]))
                            self.Logfile.insert('{0!s}: symm equivalent molecule: {1!s}'.format(xtal, files))
                            self.Logfile.insert('{0!s}: distance: {1!s}'.format(xtal, str(distance)))

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        if self.errorDict != {}:
            self.update_errorDict('General','The aforementioned PDB files were automatically changed by XCE!\nPlease check and refine them!!!')
        self.emit(QtCore.SIGNAL('show_error_dict'), self.errorDict)



class find_event_map_for_ligand(QtCore.QThread):

    def __init__(self,project_directory,xce_logfile,external_software):
        QtCore.QThread.__init__(self)
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.project_directory=project_directory
        self.external_software=external_software
        try:
            import gemmi
            self.Logfile.insert('found gemmi library in ccp4-python')
        except ImportError:
            self.external_software['gemmi'] = False
            self.Logfile.warning('cannot import gemmi; will use phenix.map_to_structure_factors instead')

    def run(self):
        self.Logfile.insert('======== checking ligand CC in event maps ========')
        for dirs in sorted(glob.glob(os.path.join(self.project_directory,'*'))):
            xtal = dirs[dirs.rfind('/')+1:]
            if os.path.isfile(os.path.join(dirs,'refine.pdb')) and \
               os.path.isfile(os.path.join(dirs,'refine.mtz')):
                self.Logfile.insert('%s: found refine.pdb' %xtal)
                os.chdir(dirs)
                try:
                    p = gemmi.read_structure('refine.pdb')
                except:
                    self.Logfile.error('gemmi library not available')
                    self.external_software['gemmi'] = False
                reso = XChemUtils.mtztools('refine.mtz').get_dmin()
                ligList = XChemUtils.pdbtools('refine.pdb').save_residues_with_resname(dirs,'LIG')
                self.Logfile.insert('%s: found %s ligands of type LIG in refine.pdb' %(xtal,str(len(ligList))))

                for maps in glob.glob(os.path.join(dirs,'*event*.native.ccp4')):
                    if self.external_software['gemmi']:
                        self.convert_map_to_sf_with_gemmi(maps,p)
                    else:
                        self.expand_map_to_p1(maps)
                        self.convert_map_to_sf(maps.replace('.ccp4','.P1.ccp4'),reso)

                summary = ''
                for lig in sorted(ligList):
                    if self.external_software['gemmi']:
                        for mtz in sorted(glob.glob(os.path.join(dirs,'*event*.native.mtz'))):
                            self.get_lig_cc(mtz,lig)
                            cc = self.check_lig_cc(mtz.replace('.mtz', '_CC.log'))
                            summary += '%s: %s LIG CC = %s (%s)\n' %(xtal,lig,cc,mtz[mtz.rfind('/')+1:])
                    else:
                        for mtz in sorted(glob.glob(os.path.join(dirs,'*event*.native*P1.mtz'))):
                            self.get_lig_cc(mtz,lig)
                            cc = self.check_lig_cc(mtz.replace('.mtz', '_CC.log'))
                            summary += '%s: %s LIG CC = %s (%s)\n' %(xtal,lig,cc,mtz[mtz.rfind('/')+1:])
                self.Logfile.insert('\nsummary of CC analysis:\n======================:\n'+summary)

    def expand_map_to_p1(self,emap):
        self.Logfile.insert('expanding map to P1: %s' %emap)
        if os.path.isfile(emap.replace('.ccp4','.P1.ccp4')):
            self.Logfile.warning('P1 map exists; skipping...')
            return
        cmd = ( 'mapmask MAPIN %s MAPOUT %s << eof\n' %(emap,emap.replace('.ccp4','.P1.ccp4'))+
                ' XYZLIM CELL\n'
                ' PAD 0.0\n'
                ' SYMMETRY 1\n'
                'eof\n' )
        os.system(cmd)

    def convert_map_to_sf(self,emap,reso):
        self.Logfile.insert('converting ccp4 map to mtz with phenix.map_to_structure_factors: %s' %emap)
        if os.path.isfile(emap.replace('.ccp4','.mtz')):
            self.Logfile.warning('mtz file of event map exists; skipping...')
            return
        cmd = ( 'module load phenix\n'
                'phenix.map_to_structure_factors %s d_min=%s\n' %(emap,reso)+
                '/bin/mv map_to_structure_factors.mtz %s' %emap.replace('.ccp4', '.mtz') )
        os.system(cmd)

    def get_lig_cc(self,mtz,lig):
        self.Logfile.insert('calculating CC for %s in %s' %(lig,mtz))
        if os.path.isfile(mtz.replace('.mtz', '_CC.log')):
            self.Logfile.warning('logfile of CC analysis exists; skipping...')
            return
        cmd = ( 'module load phenix\n'
                'phenix.get_cc_mtz_pdb %s %s > %s' % (mtz, lig, mtz.replace('.mtz', '_CC.log')) )
        os.system(cmd)

    def check_lig_cc(self,log):
        cc = 'n/a'
        if os.path.isfile(log):
            for line in open(log):
                if line.startswith('local'):
                    cc = line.split()[len(line.split()) - 1]
        else:
            self.Logfile.error('logfile does not exist: %s' %log)
        return cc

    def convert_map_to_sf_with_gemmi(self,emap,p):
        self.Logfile.insert('converting ccp4 map to mtz with gemmi map2sf: %s' %emap)
        if os.path.isfile(emap.replace('.ccp4','.mtz')):
            self.Logfile.warning('mtz file of event map exists; skipping...')
            return
        cmd = 'gemmi map2sf %s %s FWT PHWT --dmin=%s' %(emap,emap.replace('.ccp4','.mtz'),p.resolution)
        self.Logfile.insert('converting map with command:\n' + cmd)
        os.system(cmd)