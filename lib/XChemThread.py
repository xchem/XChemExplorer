# last edited: 31/07/2017, 10:30

import os, sys, glob
from datetime import datetime
from PyQt4 import QtGui, QtCore

import time
import pickle
import cPickle
import base64
import math
import subprocess
from datetime import datetime
import time
import getpass
import csv

from iotbx import mtz
from iotbx.reflection_file_reader import any_reflection_file

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
from XChemUtils import process
from XChemUtils import parse
from XChemUtils import queue
from XChemUtils import mtztools
from XChemUtils import pdbtools
from XChemUtils import helpers
from XChemUtils import reference
from XChemUtils import misc
import XChemDB
import XChemLog
import XChemMain

import iotbx.mtz


#class update_datasource_from_file_system(QtCore.QThread):
#    def __init__(self,initial_model_directory,datasource,panddas_directory,xce_logfile):
#        QtCore.QThread.__init__(self)
#        self.initial_model_directory=initial_model_directory
#        self.datasource=datasource
#        self.db=XChemDB.data_source(self.datasource)
#        self.panddas_directory=panddas_directory
#        self.Logfile=XChemLog.updateLog(xce_logfile)
#
#    def run(self):
#        self.Logfile.insert('new project directory: '+self.initial_model_directory)
#        self.Logfile.insert('updating data source from file system')
#        progress_step=1
#        if len(glob.glob(os.path.join(self.initial_model_directory,'*'))) != 0:
#            progress_step=100/float(len(glob.glob(os.path.join(self.initial_model_directory,'*'))))
#        else:
#            progress_step=1
#        progress=0
#        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
#
#        all_samples_in_datasource=self.db.get_all_samples_in_data_source_as_list()
#
#        for directory in sorted(glob.glob(os.path.join(self.initial_model_directory,'*'))):
#            try:
#                os.chdir(directory)
#            except OSError:
#                # this could happen if the user accidentaly left a file in the project directory
#                continue
#            xtal=directory[directory.rfind('/')+1:]
#            if xtal not in all_samples_in_datasource:
#                self.Logfile.insert('inserting '+xtal+' into data source')
#                self.db.execute_statement("insert into mainTable (CrystalName) values ('{0!s}');".format(xtal))
#                all_samples_in_datasource.append(xtal)
#            compoundID=str(self.db.get_value_from_field(xtal,'CompoundCode')[0])
#            db_dict={}
#            sample_dict=self.db.get_db_dict_for_sample(xtal)
#
#            dimple_path=''  # will be set to correct path if dimple.pdb is present;
#            if os.path.isfile('dimple.pdb'):
#                print ''
##                if not os.path.isfile('refine.pdb'):
##                    os.system('/bin/rm refine.pdb')         # this removes broken links that could trip the symlink
##                    os.symlink('dimple.pdb', 'refine.pdb')
#            else:
#                os.system('/bin/rm dimple.pdb 2> /dev/null')    # this makes sure that any broken link which could rerail PANDDA gets removed
#            if os.path.isfile('dimple.mtz'):
#                db_dict['DimplePathToMTZ']=os.path.realpath(os.path.join(directory,'dimple.mtz'))
#                dimple_mtz=db_dict['DimplePathToMTZ']
#                dimple_path=dimple_mtz[:dimple_mtz.rfind('/')]
#                if not os.path.isfile('refine.mtz'):
#                    os.system('/bin/rm refine.mtz')
#                    os.symlink('dimple.mtz', 'refine.mtz')
#            else:
#                os.system('/bin/rm dimple.mtz 2> /dev/null')    # this makes sure that any broken link which could rerail PANDDA gets removed
#            # this should not really be the case, but if a user does not provide an aimless logfile then that's all we can do
#            if not os.path.isfile(xtal+'.log'):
#                if os.path.isfile(xtal+'.mtz'):
#                    mtz_info=mtztools(xtal+'.mtz').get_information_for_datasource()
#                    db_dict.update(mtz_info)
#                    db_dict['DataCollectionOutcome']='success'
#            if not os.path.isfile(xtal+'.free.mtz'):
#                os.system('/bin/rm '+xtal+'.free.mtz 2> /dev/null')      # remove possible broken link
#                if os.path.isfile(os.path.join(dimple_path,'prepared2.mtz')):
#                    os.symlink(os.path.join(dimple_path,'prepared2.mtz'),xtal+'.free.mtz')
#                    db_dict['RefinementMTZfree']=xtal+'.free.mtz'
#                elif os.path.isfile(os.path.join(dimple_path,'prepared.mtz')):
#                    os.symlink(os.path.join(dimple_path,'prepared.mtz'),xtal+'.free.mtz')
#                    db_dict['RefinementMTZfree']=xtal+'.free.mtz'
#            if os.path.isfile('refine.mtz'):
#                if sample_dict['RefinementOutcome']=='None' or sample_dict['RefinementOutcome']=='':
#                    refinement_in_progress=False
#                    for dirs in glob.glob('*'):
#                        if dirs.startswith('Refine_') and os.path.isdir(dirs):
#                            db_dict['RefinementOutcome']='3 - In Refinement'
#                            refinement_in_progress=True
#                            break
#                    if not refinement_in_progress:
#                        db_dict['RefinementOutcome']='1 - Analysis Pending'
#            if os.path.isdir('compound'):
#                if sample_dict['CompoundCode']=='None' or sample_dict['CompoundCode']=='':
#                    for smiles in glob.glob('compound/*'):
#                        if smiles.endswith('smiles'):
#                            for line in open(smiles):
#                                if len(line.split()) >= 1:
#                                    db_dict['CompoundSMILES']=line.split()[0]
#                                    db_dict['CompoundCode']=smiles[smiles.rfind('/')+1:smiles.rfind('.')]
#                                    compoundID=db_dict['CompoundCode']
#                                    break
#
#            for file in glob.glob('*'):
#                if file==xtal+'.log' and os.path.isfile(file):
#                    db_dict['DataProcessingPathToLogfile']=os.path.join(directory,xtal+'.log')
#                    db_dict['DataProcessingLOGfileName']=xtal+'.log'
#                    if sample_dict['DataCollectionOutcome']=='None' or sample_dict['DataCollectionOutcome']=='':
#                        db_dict['DataCollectionOutcome']='success'
#                    aimless_results=parse().read_aimless_logfile(file)
#                    db_dict.update(aimless_results)
#                if file==xtal+'.mtz' and os.path.isfile(file):
#                    db_dict['DataProcessingPathToMTZfile']=os.path.join(directory,xtal+'.mtz')
#                    db_dict['DataProcessingMTZfileName']=xtal+'.mtz'
#                if file==xtal+'.free.mtz' and os.path.isfile(file):
#                    db_dict['RefinementMTZfree']=os.path.join(directory,xtal+'.free.mtz')
#                if file==compoundID+'.cif' and os.path.isfile(file):
#                    db_dict['RefinementCIF']=os.path.join(directory,compoundID+'.cif')
#                if file=='dimple.pdb' and os.path.isfile(file):
#                    db_dict['DimplePathToPDB']=os.path.realpath(os.path.join(directory,'dimple.pdb'))
#                    pdb_info=parse().PDBheader(file)
#                    db_dict['DimpleRcryst']=pdb_info['Rcryst']
#                    db_dict['DimpleRfree']=pdb_info['Rfree']
#                    db_dict['DimpleResolutionHigh']=pdb_info['ResolutionHigh']
#                if file=='refine.pdb' and os.path.isfile(file):
##                    if sample_dict['RefinementOutcome']=='None' or sample_dict['RefinementOutcome']=='':
##                        db_dict['RefinementOutcome']='3 - In Refinement'
#                    db_dict['RefinementPDB_latest']=os.path.realpath(os.path.join(directory,'refine.pdb'))
#                    pdb_info=parse().PDBheader(file)
#                    db_dict['RefinementRcryst']=pdb_info['Rcryst']
#                    db_dict['RefinementRfree']=pdb_info['Rfree']
#                    db_dict['RefinementSpaceGroup']=pdb_info['SpaceGroup']
#                    db_dict['RefinementRmsdBonds']=pdb_info['rmsdBonds']
#                    db_dict['RefinementRmsdAngles']=pdb_info['rmsdAngles']
#                if file=='refine.mtz' and os.path.isfile(file):
#                    db_dict['RefinementMTZ_latest']=os.path.realpath(os.path.join(directory,'refine.mtz'))
#
#            if db_dict != {}:
#                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'updating datasource for '+xtal)
#                self.db.update_data_source(xtal,db_dict)
#
#            # also need to update PANDDA table...
#            pandda_models=self.db.execute_statement("select CrystalName,PANDDA_site_index,PANDDA_site_spider_plot,PANDDA_site_event_map from panddaTable where CrystalName='{0!s}'".format(xtal))
#            if not pandda_models == []:
#                for entry in pandda_models:
#                    db_pandda_dict={}
#                    db_pandda_dict['PANDDA_site_index']=entry[1]
#                    db_pandda_dict['PANDDApath']=self.panddas_directory
#                    if entry[3] is not None:
#                        event_map=os.path.join(self.initial_model_directory,xtal,entry[3].split('/')[len(entry[3].split('/'))-1])
#                        if os.path.isfile(event_map):
#                            db_pandda_dict['PANDDA_site_event_map']=event_map
#                    if entry[2] is not None:
#                        spider_plot=os.path.join(self.initial_model_directory,xtal,entry[2].split('/')[len(entry[2].split('/'))-3],entry[2].split('/')[len(entry[2].split('/'))-2],entry[2].split('/')[len(entry[2].split('/'))-1])
#                        if os.path.isfile(spider_plot):
#                            db_pandda_dict['PANDDA_site_spider_plot']=spider_plot
##                            db_pandda_dict['RefinementOutcome']='3 - In Refinement'    # just in case; presence of a spider plot definitely signals that refinement happened
#                                                                                        # should probably be not updated! Will overwrite CompChem ready
#                    self.Logfile.insert('updating panddaTable for xtal: {0!s}, site: {1!s}'.format(entry[0], entry[1]))
#                    self.db.update_insert_panddaTable(xtal,db_pandda_dict)
#
#            progress += progress_step
#            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
#
#        self.Logfile.insert('datasource update finished')
#

class synchronise_db_and_filesystem(QtCore.QThread):

    '''
        - remove broken links
        - insert new samples in DB
        - update data for existing samples
    '''

    def __init__(self,initial_model_directory,datasource,panddas_directory,xce_logfile,mode):
        QtCore.QThread.__init__(self)
        self.initial_model_directory=initial_model_directory
        self.datasource=datasource
        self.db=XChemDB.data_source(self.datasource)
        self.all_samples_in_datasource=self.db.get_all_samples_in_data_source_as_list()
        self.panddas_directory=panddas_directory
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.mode=mode

    def run(self):
        self.Logfile.insert('synchronising database and filesystem')
        self.Logfile.insert('current project directory: '+self.initial_model_directory)

        #
        # get list of xtals
        #

        self.xtal_list = []
        progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
        if self.mode != 'project_directory':
            # if only a single xtal is synched, then self.mode==xtalID
            self.Logfile.insert('synchronising '+self.mode+' only')
            self.xtal_list.append(self.mode)
        else:
            if len(glob.glob(os.path.join(self.initial_model_directory,'*'))) != 0:
                progress_step=100/float(len(glob.glob(os.path.join(self.initial_model_directory,'*'))))
            self.Logfile.insert('found '+str(len(glob.glob(os.path.join(self.initial_model_directory,'*'))))+' samples in project directory')
            for directory in sorted(glob.glob(os.path.join(self.initial_model_directory,'*'))):
                try:
                    os.chdir(directory)
                except OSError:
                    # this could happen if the user accidentaly left a file in the project directory
                    continue
                xtal=directory[directory.rfind('/')+1:]
                self.xtal_list.append(xtal)

        #
        # go through list
        #

        for xtal in self.xtal_list:
            self.Logfile.insert('directory name: '+xtal+' = sampleID in database')
            os.chdir(os.path.join(self.initial_model_directory,xtal))
            if xtal not in self.all_samples_in_datasource:
                self.Logfile.insert('sampleID not found in database: inserting '+xtal)
                self.db.execute_statement("insert into mainTable (CrystalName) values ('{0!s}');".format(xtal))
                self.all_samples_in_datasource.append(xtal)

            db_dict=self.db.get_db_dict_for_sample(xtal)

            db_dict['ProjectDirectory'] = self.initial_model_directory

            db_dict=self.sync_data_processing(xtal,db_dict)

            db_dict=self.sync_dimple_results(xtal,db_dict)

            db_dict=self.sync_compound_information(db_dict)

            db_dict=self.sync_refinement_results(xtal,db_dict)

            if db_dict != {}:
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'updating datasource for '+xtal)
                self.db.update_data_source(xtal,db_dict)

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.Logfile.insert('database mainTable update finished')
        self.Logfile.insert('updating panddaTable')
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'updating panddaTable')
#        self.sync_pandda_table()
        self.sync_pandda_table_NEW()
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'database panddaTable update finished')
        self.Logfile.insert('database panddaTable update finished')

        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))


    def change_absolute_to_relative_links(self,target,filename):
        os.unlink(filename)
        os.symlink(os.path.relpath(target),filename)
        self.Logfile.insert('{0!s} -> {1!s}'.format(os.readlink(filename), target))

    def find_file(self,filename,xtal):
        found_file=False
        for files in glob.glob('*'):
            if files == filename:
                # at this point, this could ba a file or a (broken) symbolic link
                try:
                    link=os.readlink(filename)
                    if link.startswith('/'):
                        if os.path.isfile(link.replace(link[:link.find(xtal)],self.initial_model_directory+'/')):
                            target=link.replace(link[:link.find(xtal)],self.initial_model_directory+'/')
                            found_file=True
                            self.change_absolute_to_relative_links(target,filename)
                        elif os.path.isfile(files):
                            # this will leave the absolute path
                            found_file=True
                        else:
                            self.Logfile.insert('removing broken link: '+filename)
                            os.system('/bin/rm {0!s} 2> /dev/null'.format(filename))
                    else:
                        if not os.path.isfile(link):
                            self.Logfile.insert('removing broken link: '+filename)
                            os.system('/bin/rm {0!s} 2> /dev/null'.format(filename))
                        else:
                            found_file=True
                except OSError:
                    if os.path.isfile(filename):
                        found_file=True
                break
        return found_file




    def sync_data_processing(self,xtal,db_dict):

        # AIMLESS logfile

        # in case the MTZ file which is used for refinement is different to the one used for refinement
        if os.path.isfile('refine.mtz'):
            if os.path.isfile(xtal+'.free.mtz'):
                freeMTZ=mtztools(xtal+'.free.mtz')
                nREFfree=freeMTZ.get_number_measured_reflections()
                if os.path.isfile(xtal+'.mtz'):
                    procMTZ=mtztools(xtal+'.mtz')
                    nREF=procMTZ.get_number_measured_reflections()
                    CC,errorMessage=freeMTZ.calculate_correlaton_between_intensities_in_mtzfiles(xtal+'.mtz')
                    self.Logfile.insert('%s: calculating CC between %s.free.mtz (%s refl) and %s.mtz (%s refl): %s' %(xtal,xtal,str(nREFfree),xtal,str(nREF),str(CC)))
                    if errorMessage != '':
                        self.Logfile.insert('pointless failed with the following error: %s' %errorMessage)

                    try:
                        if float(CC) < 0.999:
                            self.Logfile.insert('correlation coefficient between the two files is below 0.999; will try to understand from dimple.log which one was used for initial map calculation')
                            if os.path.isfile('dimple/dimple_rerun_on_selected_file/dimple/dimple.log'):
                                foundLine=False
                                mtzin=''
                                for line in open('dimple/dimple_rerun_on_selected_file/dimple/dimple.log'):
                                    if foundLine:
                                        mtzin=line.replace(' ','').replace('\n','').replace('\r','')
                                        self.Logfile.insert('%s was used for inital map calculation' %mtzin)
                                        break
                                    if line.startswith(' --no-cleanup'):
                                        foundLine=True

                                if os.path.isfile(mtzin):
                                    self.Logfile.insert('%s: mtzfile used for refinement is not the same as the one chosen from autoprocessing' %xtal)
                                    self.Logfile.insert('%s: current mtzfile after autoprocessing: %s' %(xtal,os.path.realpath(xtal+'.mtz')))
                                    self.Logfile.insert('%s: removing links for %s.mtz/%s.log' %(xtal,xtal,xtal))
                                    os.system('/bin/rm %s.mtz 2> /dev/null' %xtal)
                                    os.system('/bin/rm %s.log 2> /dev/null' %xtal)
                                    self.Logfile.insert('linking %s to %s.mtz' %(os.path.relpath(mtzin),xtal))
                                    os.symlink(os.path.relpath(mtzin),xtal+'.mtz')
                                    for logfile in glob.glob(os.path.join(mtzin[:mtzin.rfind('/')],'*log')):
                                        self.Logfile.insert('linking %s to %s.log' %(os.path.relpath(logfile),xtal))
                                        os.symlink(os.path.relpath(logfile),xtal+'.log')
                                        break

                    except ValueError:
                        self.Logfile.insert('something went wrong: calculated CC value does not seem to be a floating point number')

        found_logfile=False
        if os.path.isfile(xtal+'.log'):
            found_logfile=True
#            db_dict['DataProcessingPathToLogfile']=os.path.realpath(xtal+'.log').replace(os.getcwd()+'/','')
            db_dict['DataProcessingPathToLogfile']=os.path.realpath(xtal+'.log')
            db_dict['DataProcessingLOGfileName']=xtal+'.log'
            if db_dict['DataCollectionOutcome']=='None' or db_dict['DataCollectionOutcome']=='':
                db_dict['DataCollectionOutcome']='success'
            aimless_results=parse().read_aimless_logfile(xtal+'.log')
            db_dict.update(aimless_results)
        else:
            db_dict['DataProcessingPathToLogfile']=''
            db_dict['DataProcessingLOGfileName']=''

        # MTZ file

#        found_mtzfile=self.find_file(xtal+'.mtz',xtal)
#        if found_mtzfile:
        if os.path.isfile(xtal+'.mtz'):
            db_dict['DataProcessingPathToMTZfile']=os.path.realpath(xtal+'.mtz')
#            db_dict['DataProcessingPathToMTZfile']=os.path.realpath(xtal+'.mtz').replace(os.getcwd()+'/','')
            db_dict['DataProcessingMTZfileName']=xtal+'.mtz'
            if not found_logfile:
                mtz_info=mtztools(xtal+'.mtz').get_information_for_datasource()
                db_dict.update(mtz_info)
                db_dict['DataCollectionOutcome']='success'
        else:
            db_dict['DataProcessingPathToMTZfile']=''
            db_dict['DataProcessingMTZfileName']=''

        return db_dict


    def sync_dimple_results(self,xtal,db_dict):

        # DIMPLE pdb

#        found_dimple_pdb=self.find_file('dimple.pdb',xtal)
#        if found_dimple_pdb:
        if os.path.isfile('dimple.pdb'):
#            db_dict['DimplePathToPDB']=os.path.realpath('dimple.pdb').replace(os.getcwd()+'/','')
            db_dict['DimplePathToPDB']=os.path.realpath('dimple.pdb')
            pdb_info=parse().PDBheader('dimple.pdb')
            db_dict['DimpleRcryst']=pdb_info['Rcryst']
            db_dict['DimpleRfree']=pdb_info['Rfree']
            db_dict['DimpleResolutionHigh']=pdb_info['ResolutionHigh']
            db_dict['DimpleStatus']='finished'
        else:
            db_dict['DimplePathToPDB']=''
            db_dict['DimpleRcryst']=''
            db_dict['DimpleRfree']=''
            db_dict['DimpleResolutionHigh']=''
            db_dict['DimpleStatus']='pending'
        # DIMPLE mtz

        dimple_path=''
#        found_dimple_mtz=self.find_file('dimple.mtz',xtal)
#        if found_dimple_mtz:
        if os.path.isfile('dimple.mtz'):
#            db_dict['DimplePathToMTZ']=os.path.realpath('dimple.mtz').replace(os.getcwd()+'/','')
            db_dict['DimplePathToMTZ']=os.path.realpath('dimple.mtz')
            dimple_mtz=db_dict['DimplePathToMTZ']
            dimple_path=dimple_mtz[:dimple_mtz.rfind('/')]
        else:
            db_dict['DimplePathToMTZ']=''
            db_dict['DimpleStatus']='pending'

        if os.path.isfile(os.path.join(dimple_path,'dimple','dimple_rerun_on_selected_file','dimple_run_in_progress')):
            db_dict['DimpleStatus']='running'

        # MTZ free file

#        found_free_mtz=self.find_file(xtal+'.free.mtz',xtal)
#        if found_free_mtz:
        if os.path.isfile(xtal+'.free.mtz'):
            db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz')
#            db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz').replace(os.getcwd()+'/','')
        else:
            db_dict['RefinementMTZfree']=''
            os.system('/bin/rm %s.free.mtz 2> /dev/null' %xtal)
            if os.path.isfile(os.path.join(dimple_path,'prepared2.mtz')):
                os.symlink(os.path.relpath(os.path.join(dimple_path,'prepared2.mtz')),xtal+'.free.mtz')
#                db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz').replace(os.getcwd()+'/','')
                db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz')
#                db_dict['RefinementMTZfree']=xtal+'.free.mtz'
            elif os.path.isfile(os.path.join(dimple_path,'prepared.mtz')):
                os.symlink(os.path.relpath(os.path.join(dimple_path,'prepared.mtz')),xtal+'.free.mtz')
                db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz')
#                db_dict['RefinementMTZfree']=os.path.realpath(xtal+'.free.mtz').replace(os.getcwd()+'/','')
#               db_dict['RefinementMTZfree']=xtal+'.free.mtz'

        return db_dict


    def sync_compound_information(self,db_dict):
        # only update database if SMILES or compoundID field is blank!

        compoundID=db_dict['CompoundCode']
        if compoundID=='None' or compoundID=='':
            if os.path.isdir('compound'):
                for smiles in glob.glob('compound/*'):
                    if smiles.endswith('smiles'):
                        for line in open(smiles):
                            if len(line.split()) >= 1:
                                db_dict['CompoundCode']=smiles[smiles.rfind('/')+1:smiles.rfind('.')]
                                compoundID=db_dict['CompoundCode']
                                break

        if os.path.isfile(compoundID+'.cif') and os.path.getsize(compoundID+'.cif') > 20:
            db_dict['RefinementCIF']=os.path.realpath(compoundID+'.cif').replace(os.getcwd()+'/','')
            db_dict['RefinementCIFStatus']='restraints generated'
        else:
            os.system('/bin/rm {0!s}.cif 2> /dev/null'.format(compoundID))
            os.system('/bin/rm compound/{0!s}.cif 2> /dev/null'.format(compoundID))
            db_dict['RefinementCIF']=''
            db_dict['RefinementCIFStatus']='pending'

        smilesDB=db_dict['CompoundSMILES']
        smiles_found=True
        if smilesDB=='None' or smilesDB=='':
            smiles_found=False
            if os.path.isdir('compound'):
                for smiles in glob.glob('compound/*'):
                    if smiles.endswith('smiles'):
                        for line in open(smiles):
                            if len(line.split()) >= 1:
                                db_dict['CompoundSMILES']=line.split()[0]
                                smilesDB=db_dict['CompoundSMILES']
                                smiles_found=True
                                break

        if not smiles_found:
            db_dict['RefinementCIFStatus']='missing smiles'

        if not os.path.isfile(compoundID+'.pdb') or  os.path.getsize(compoundID+'.pdb') < 20:
            os.system('/bin/rm {0!s}.pdb 2> /dev/null'.format(compoundID))
            os.system('/bin/rm compound/{0!s}.pdb 2> /dev/null'.format(compoundID))

        if not os.path.isfile(compoundID+'.png') or  os.path.getsize(compoundID+'.png') < 20:
            os.system('/bin/rm {0!s}.png 2> /dev/null'.format(compoundID))
            os.system('/bin/rm compound/{0!s}.png 2> /dev/null'.format(compoundID))

        return db_dict


    def sync_refinement_results(self,xtal,db_dict):

        #
        # REFINE pdb
        #

#        found_refine_pdb=self.find_file('refine.pdb',xtal)
#        if found_refine_pdb and not 'dimple' in os.path.realpath('refine.pdb'):

#        if os.path.isfile('refine.pdb') and not 'dimple' in os.path.realpath('refine.pdb'):
        if os.path.isfile('refine.pdb'):
            db_dict['RefinementPDB_latest']=os.path.realpath('refine.pdb')
#            db_dict['RefinementPDB_latest']=os.path.realpath('refine.pdb').replace(os.getcwd()+'/','')
            db_dict['RefinementStatus']='finished'
            pdb_info=parse().dict_for_datasource_update('refine.pdb')
            db_dict.update(pdb_info)
            if db_dict['RefinementOutcome']=='None' or db_dict['RefinementOutcome']=='':
                db_dict['RefinementOutcome']='3 - In Refinement'
            elif str(db_dict['RefinementOutcome']).startswith('1'):
                db_dict['RefinementOutcome']='3 - In Refinement'
            elif str(db_dict['RefinementOutcome']).startswith('2'):
                db_dict['RefinementOutcome']='3 - In Refinement'
        else:
            db_dict['RefinementPDB_latest']=''
            db_dict['RefinementStatus']='pending'
            db_dict['RefinementOutcome']='1 - Analysis Pending'
            os.system('/bin/rm refine.pdb 2> /dev/null')

        if os.path.isfile('REFINEMENT_IN_PROGRESS'):
            db_dict['RefinementStatus']='running'

        #
        # REFINE bound pdb
        #

        if os.path.isfile('refine.bound.pdb'):
            db_dict['RefinementBoundConformation']=os.path.realpath('refine.bound.pdb')
#            db_dict['RefinementBoundConformation']='refine.bound.pdb'
        else:
#            if os.path.isfile('refine.pdb'):
#                os.system("giant.strip_conformations pdb=refine.pdb suffix='.bound.pdb'")
#                if os.path.isfile('refine.bound.pdb'):
#                    db_dict['RefinementBoundConformation']=os.path.realpath('refine.bound.pdb')
#                else:
#                    db_dict['RefinementBoundConformation']=''
#            else:
#                db_dict['RefinementBoundConformation']=''
            db_dict['RefinementBoundConformation']=''

        #
        # REFINE mtz
        #

#        found_refine_mtz=self.find_file('refine.mtz',xtal)
#        if found_refine_mtz and not 'dimple' in os.path.realpath('refine.mtz'):
#        if os.path.isfile('refine.mtz') and not 'dimple' in os.path.realpath('refine.mtz'):
        if os.path.isfile('refine.mtz'):
            db_dict['RefinementMTZ_latest']=os.path.realpath('refine.mtz')
#            db_dict['RefinementMTZ_latest']=os.path.realpath('refine.mtz').replace(os.getcwd()+'/','')
        else:
            db_dict['RefinementMTZ_latest']=''
            os.system('/bin/rm refine.mtz 2> /dev/null')

        return db_dict


    def find_apo_structures_for_PanDDA(self,panddaPATH):

        # first check if structure is already present in DB and if so if all the
        # information concur

        # need to update pandda directory for every exported structure so that
        # we know where to look for the pandda.log file that contains the relevant information

        # update CrystalName_of_pandda_input in DB

        # in DB: update StructureType field accordingly

        # newer pandda versions seem to have severl copies of pandda.log with names like
        # pandda-2016-09-01-2139.log
        panddaLog=glob.glob(os.path.join(panddaPATH,'pandda*log'))
        panddaLog.sort(key=os.path.getmtime)

        panddaVersion='unknown'
        readindApoStructures = False
        apoStructures = []
        apoStructureDict = {}
        apoString=''
        for files in panddaLog:
            for line in open(files):
                if line.startswith('-  Pandda Version'):
                    if len(line.split()) >= 4:
                        panddaVersion=line.split()[3]
                if 'No Statistical Maps Found:' in line:
                    readindApoStructures=True
                if readindApoStructures:
                    if 'Pickling Object: processed_datasets' in line:
                        if line.split() >= 2:
                            # e.g. line.split() -> ['Pickling', 'Object:', 'processed_datasets/NUDT22A-x0055/pickles/dataset.pickle']
                            xtal=line.split()[2].split('/')[1]
                            if os.path.isfile(os.path.join(panddaPATH,'processed_datasets',xtal,xtal+'-pandda-input.pdb')):
                                apoStructures.append(xtal)
                                apoString+=xtal+';'
                if 'Pre-existing statistical maps (from previous runs) have been found and will be reused:' in line:
                    readindApoStructures=False
            apoStructureDict[panddaPATH]=apoStructures

        return apoString[:-1]

#    def sync_pandda_table(self):
#        progress_step=1
#        progress=0
#        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
#
#        # also need to update PANDDA table...
#        pandda_models=self.db.execute_statement("select CrystalName,PANDDA_site_index,PANDDA_site_event_index,PANDDA_site_x,PANDDA_site_y,PANDDA_site_z,PANDDApath,ApoStructures from panddaTable")
#        if len(pandda_models) > 0:
#            progress_step=100/float(len(pandda_models))
#
#        if pandda_models != []:
#            for entry in pandda_models:
#                db_pandda_dict={}
#                xtal=entry[0]
#                site_index=entry[1]
#                event_index=entry[2]
#                panddaPATH=entry[6]
#                apoStructures=entry[7]
#                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'checking %s -> site %s -> event %s ' %(xtal,site_index,event_index))
#                try:
#                    event_x = float(str(entry[3]))
#                    event_y = float(str(entry[4]))
#                    event_z = float(str(entry[5]))
#                except ValueError:
#                    pass
#
#                # do not update pandda path since this one is updated during pandda export!
#                # instead try to get apo semi-colon separated list of apo structures that were used to
#                # calculate event maps; but only if field is blank!
##                db_pandda_dict['PANDDApath']=self.panddas_directory
#                if str(apoStructures)=='None' or apoStructures=='':
#                    if panddaPATH != 'None' or panddaPATH != '':
#                        self.Logfile.insert('trying to find which apo structures were used to calculate the event maps in '+panddaPATH)
#                        db_pandda_dict['ApoStructures']=self.find_apo_structures_for_PanDDA(panddaPATH)
#                    else:
#                        self.Logfile.insert('pandda path for '+xtal+' is empty in database')
#
#                # event map
#
#                found_event_map=False
#                db_pandda_dict['PANDDA_site_event_map']=''
#                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*ccp4')):
#                    filename=file[file.rfind('/')+1:]
#                    if filename.startswith(xtal+'-event_'+event_index) and filename.endswith('map.native.ccp4'):
#                        event_map=file
#                        db_pandda_dict['PANDDA_site_event_map']=os.path.realpath(event_map)
##                        db_pandda_dict['PANDDA_site_event_map']=os.path.realpath(event_map).replace(os.getcwd()+'/','')
#                        found_event_map=True
#                        break
#                if not found_event_map:
#                    db_pandda_dict['PANDDA_site_event_map']=''
#
#                db_pandda_dict['PANDDA_site_initial_model']=''
#                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*pdb')):
#                    filename=file[file.rfind('/')+1:]
#                    if filename.endswith('-ensemble-model.pdb'):
#                        db_pandda_dict['PANDDA_site_initial_model']=os.path.realpath(file).replace(os.getcwd()+'/','')
#                        break
#
#                db_pandda_dict['PANDDA_site_initial_mtz']=''
#                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*mtz')):
#                    filename=file[file.rfind('/')+1:]
#                    if filename.endswith('pandda-input.mtz'):
#                        db_pandda_dict['PANDDA_site_initial_mtz']=os.path.realpath(file).replace(os.getcwd()+'/','')
#                        break
#
#
#                db_pandda_dict['PANDDA_site_ligand_resname'] = ''
#                db_pandda_dict['PANDDA_site_ligand_chain'] = ''
#                db_pandda_dict['PANDDA_site_ligand_sequence_number'] = ''
#                db_pandda_dict['PANDDA_site_ligand_altLoc'] = ''
#                db_pandda_dict['PANDDA_site_ligand_placed'] = 'False'
#                db_pandda_dict['PANDDA_site_spider_plot'] = ''
#                db_pandda_dict['PANDDA_site_ligand_id'] = ''
#
#                db_pandda_dict['PANDDA_site_occupancy'] = ''
#                db_pandda_dict['PANDDA_site_B_average'] = ''
#                db_pandda_dict['PANDDA_site_B_ratio_residue_surroundings'] = ''
#                db_pandda_dict['PANDDA_site_rmsd'] = ''
#                db_pandda_dict['PANDDA_site_RSCC'] = ''
#                db_pandda_dict['PANDDA_site_RSR'] = ''
#                db_pandda_dict['PANDDA_site_RSZD'] = ''
#
#                if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'refine.pdb')):
#                    ligands_in_file=pdbtools(os.path.join(self.initial_model_directory,xtal,'refine.pdb')).find_xce_ligand_details()
#
#                    for ligand in ligands_in_file:
#                        residue_name=   ligand[0]
#                        residue_chain=  ligand[1]
#                        residue_number= ligand[2]
#                        residue_altLoc= ligand[3]
#                        residue_xyz = pdbtools(os.path.join(self.initial_model_directory,xtal,'refine.pdb')).get_center_of_gravity_of_residue_ish(residue_chain, residue_number)
#                        distance = misc().calculate_distance_between_coordinates(residue_xyz[0], residue_xyz[1],residue_xyz[2],event_x, event_y,event_z)
#                        # if coordinate of ligand and event are closer than 7A, then we assume they belong together
#                        if distance < 7:
#                            db_pandda_dict['PANDDA_site_ligand_resname'] = residue_name
#                            db_pandda_dict['PANDDA_site_ligand_chain'] = residue_chain
#                            db_pandda_dict['PANDDA_site_ligand_sequence_number'] = residue_number
#                            db_pandda_dict['PANDDA_site_ligand_altLoc'] = residue_altLoc
#                            db_pandda_dict['PANDDA_site_ligand_placed'] = 'True'
#                            db_pandda_dict['PANDDA_site_ligand_id']=residue_name+'-'+residue_chain+'-'+residue_number
#
#                            if xtal+'/Refine_' in os.path.realpath(os.path.join(self.initial_model_directory,xtal,'refine.pdb')):
#                                tmp=os.path.realpath(os.path.join(self.initial_model_directory,xtal,'refine.pdb'))
#                                spider_plot=os.path.join(tmp[:tmp.rfind('/')],'residue_plots',residue_name+'-'+residue_chain+'-'+residue_number+'.png')
#                                if os.path.isfile(spider_plot):
#                                    db_pandda_dict['PANDDA_site_spider_plot']=os.path.realpath(spider_plot)
##                                    db_pandda_dict['PANDDA_site_spider_plot']=os.path.realpath(spider_plot).replace(os.getcwd()+'/','')
#                                if os.path.isfile(os.path.join(tmp[:tmp.rfind('/')],'residue_scores.csv')):
#
#                                    with open(os.path.join(tmp[:tmp.rfind('/')],'residue_scores.csv'), 'rb') as csv_import:
#                                        csv_dict = csv.DictReader(csv_import)
#                                        for i, line in enumerate(csv_dict):
#                                            residueNameChainNumber = line['']
#                                            if residueNameChainNumber == residue_name+'-'+residue_chain+'-'+residue_number:
#                                                db_pandda_dict['PANDDA_site_occupancy'] = line['Occupancy']
#                                                db_pandda_dict['PANDDA_site_B_average'] = line['Average B-factor (Residue)']
#                                                db_pandda_dict['PANDDA_site_B_ratio_residue_surroundings'] = line['Surroundings B-factor Ratio']
#                                                db_pandda_dict['PANDDA_site_rmsd'] = line['Model RMSD']
#                                                db_pandda_dict['PANDDA_site_RSCC'] = line['RSCC']
#                                                db_pandda_dict['PANDDA_site_RSR'] = line['RSR']
#                                                db_pandda_dict['PANDDA_site_RSZD'] = line['RSZD']
#                            break
#
#                if db_pandda_dict != {}:
##                    self.db.update_panddaTable(xtal, site_index, db_pandda_dict)
#                    self.db.update_site_event_panddaTable(xtal, site_index, event_index, db_pandda_dict)
#                    self.Logfile.insert('updating panddaTable for xtal: %s, site: %s' %(xtal,site_index))
#
#                progress += progress_step
#                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
#


    def sync_pandda_table_NEW(self):
        progress_step=1
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        # also need to update PANDDA table...
        pandda_models=self.db.execute_statement("select CrystalName,PANDDA_site_index,PANDDA_site_event_index,PANDDA_site_x,PANDDA_site_y,PANDDA_site_z,PANDDApath,ApoStructures from panddaTable")
        if len(pandda_models) > 0:
            progress_step=100/float(len(pandda_models))
        else:
            self.Logfile.warning('panddaTable seems to be empty!')

        if pandda_models:
            for entry in pandda_models:
                db_pandda_dict={}
                xtal=entry[0]
                site_index=entry[1]
                event_index=entry[2]
                panddaPATH=entry[6]
                apoStructures=entry[7]
                self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'checking {0!s} -> site {1!s} -> event {2!s} '.format(xtal, site_index, event_index))
                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)
                progress += progress_step

                try:
                    event_x = float(str(entry[3]))
                    event_y = float(str(entry[4]))
                    event_z = float(str(entry[5]))
                except ValueError:
                    pass

                # do not update pandda path since this one is updated during pandda export!
                # instead try to get apo semi-colon separated list of apo structures that were used to
                # calculate event maps; but only if field is blank!
#                db_pandda_dict['PANDDApath']=self.panddas_directory
                if str(apoStructures)=='None' or apoStructures=='':
                    if panddaPATH != 'None' or panddaPATH != '':
                        self.Logfile.insert('trying to find which apo structures were used to calculate the event maps in '+panddaPATH)
                        db_pandda_dict['ApoStructures']=self.find_apo_structures_for_PanDDA(panddaPATH)
                    else:
                        self.Logfile.insert('pandda path for '+xtal+' is empty in database')

                # event map

                found_event_map=False
                db_pandda_dict['PANDDA_site_event_map']=''
                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*ccp4')):
                    filename=file[file.rfind('/')+1:]
                    if filename.startswith(xtal+'-event_'+event_index) and filename.endswith('map.native.ccp4'):
                        event_map=file
                        db_pandda_dict['PANDDA_site_event_map']=os.path.realpath(event_map)
#                        db_pandda_dict['PANDDA_site_event_map']=os.path.realpath(event_map).replace(os.getcwd()+'/','')
                        found_event_map=True
                        break
                if not found_event_map:
                    db_pandda_dict['PANDDA_site_event_map']=''

                db_pandda_dict['PANDDA_site_initial_model']=''
                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*pdb')):
                    filename=file[file.rfind('/')+1:]
                    if filename.endswith('-ensemble-model.pdb'):
                        db_pandda_dict['PANDDA_site_initial_model']=os.path.realpath(file).replace(os.getcwd()+'/','')
                        break

                db_pandda_dict['PANDDA_site_initial_mtz']=''
                for file in glob.glob(os.path.join(self.initial_model_directory,xtal,'*mtz')):
                    filename=file[file.rfind('/')+1:]
                    if filename.endswith('pandda-input.mtz'):
                        db_pandda_dict['PANDDA_site_initial_mtz']=os.path.realpath(file).replace(os.getcwd()+'/','')
                        break


                db_pandda_dict['PANDDA_site_ligand_resname'] = ''
                db_pandda_dict['PANDDA_site_ligand_chain'] = ''
                db_pandda_dict['PANDDA_site_ligand_sequence_number'] = ''
                db_pandda_dict['PANDDA_site_ligand_altLoc'] = ''
                db_pandda_dict['PANDDA_site_ligand_placed'] = 'False'
                db_pandda_dict['PANDDA_site_spider_plot'] = ''
                db_pandda_dict['PANDDA_site_ligand_id'] = ''

                db_pandda_dict['PANDDA_site_occupancy'] = ''
                db_pandda_dict['PANDDA_site_B_average'] = ''
                db_pandda_dict['PANDDA_site_B_ratio_residue_surroundings'] = ''
                db_pandda_dict['PANDDA_site_rmsd'] = ''
                db_pandda_dict['PANDDA_site_RSCC'] = ''
                db_pandda_dict['PANDDA_site_RSR'] = ''
                db_pandda_dict['PANDDA_site_RSZD'] = ''

                if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'refine.pdb')):
                    ligands_in_file=pdbtools(os.path.join(self.initial_model_directory,xtal,'refine.pdb')).find_xce_ligand_details()
                    if not ligands_in_file:
                        self.Logfile.warning('{0!s}: could not find any ligands in refine.pdb'.format(xtal))
                        continue
                    else:
                        self.Logfile.insert('{0!s}: found the following ligands in refine.pdb: {1!s}'.format(xtal, str(ligands_in_file)))

                    distanceList=[]
                    for ligand in ligands_in_file:
                        residue_name=   ligand[0]
                        residue_chain=  ligand[1]
                        residue_number= ligand[2]
                        residue_altLoc= ligand[3]
                        residue_xyz = pdbtools(os.path.join(self.initial_model_directory,xtal,'refine.pdb')).get_center_of_gravity_of_residue_ish(residue_chain, residue_number)
                        distance = misc().calculate_distance_between_coordinates(residue_xyz[0], residue_xyz[1],residue_xyz[2],event_x, event_y,event_z)
                        distanceList.append([distance,residue_name,residue_chain,residue_number,residue_altLoc])
                        self.Logfile.insert('{0!s}: calculating distance between event and ligand ({1!s} {2!s} {3!s}): {4!s}'.format(xtal, residue_name, residue_chain, residue_number, str(distance)))

                    # now take the ligand that is closest to the event
                    try:
                        smallestDistance = min(distanceList, key=lambda x: x[0])
                    except ValueError:
                        self.Logfile.error('could not determine smallest distance between current ligand pandda events')
                        continue
                    distance =          smallestDistance[0]
                    residue_name =      smallestDistance[1]
                    residue_chain =     smallestDistance[2]
                    residue_number =    smallestDistance[3]
                    residue_altLoc =    smallestDistance[4]
                    self.Logfile.insert('{0!s}: ligand with the shorted distance ({1!s}A) to the current event (id: {2!s}): {3!s} {4!s} {5!s} {6!s}'.format(xtal, str(distance), event_index, residue_name, residue_chain, residue_number, residue_altLoc))
                    db_pandda_dict['PANDDA_site_ligand_resname'] = residue_name
                    db_pandda_dict['PANDDA_site_ligand_chain'] = residue_chain
                    db_pandda_dict['PANDDA_site_ligand_sequence_number'] = residue_number
                    db_pandda_dict['PANDDA_site_ligand_altLoc'] = residue_altLoc
                    db_pandda_dict['PANDDA_site_ligand_placed'] = 'True'
                    db_pandda_dict['PANDDA_site_ligand_id']=residue_name+'-'+residue_chain+'-'+residue_number

                    if xtal+'/Refine_' in os.path.realpath(os.path.join(self.initial_model_directory,xtal,'refine.pdb')):
                        tmp=os.path.realpath(os.path.join(self.initial_model_directory,xtal,'refine.pdb'))
                        spider_plot=os.path.join(tmp[:tmp.rfind('/')],'residue_plots',residue_chain+'-'+residue_number+'.png').replace(' ','')
                        if os.path.isfile(spider_plot):
                            db_pandda_dict['PANDDA_site_spider_plot']=os.path.realpath(spider_plot)
                        if os.path.isfile(os.path.join(tmp[:tmp.rfind('/')],'residue_scores.csv')):
                            with open(os.path.join(tmp[:tmp.rfind('/')],'residue_scores.csv'), 'rb') as csv_import:
                                csv_dict = csv.DictReader(csv_import)
                                for i, line in enumerate(csv_dict):
                                    residueNameChainNumber = line['']
                                    if residueNameChainNumber == residue_chain+'-'+residue_number:
                                        db_pandda_dict['PANDDA_site_occupancy'] = line['Occupancy']
                                        db_pandda_dict['PANDDA_site_B_average'] = line['Average B-factor (Residue)']
                                        db_pandda_dict['PANDDA_site_B_ratio_residue_surroundings'] = line['Surroundings B-factor Ratio']
                                        db_pandda_dict['PANDDA_site_rmsd'] = line['Model RMSD']
                                        db_pandda_dict['PANDDA_site_RSCC'] = line['RSCC']
                                        db_pandda_dict['PANDDA_site_RSR'] = line['RSR']
                                        db_pandda_dict['PANDDA_site_RSZD'] = line['RSZD']
#                    break

                if db_pandda_dict != {}:
#            self.db.update_panddaTable(xtal, site_index, db_pandda_dict)
                    self.db.update_site_event_panddaTable(xtal, site_index, event_index, db_pandda_dict)
                    self.Logfile.insert('updating panddaTable for xtal: {0!s}, site: {1!s}'.format(xtal, site_index))

#                progress += progress_step
#                self.emit(QtCore.SIGNAL('update_progress_bar'), progress)





class create_png_and_cif_of_compound(QtCore.QThread):
    def __init__(self,external_software,initial_model_directory,compound_list,database_directory,data_source_file,todo,ccp4_scratch_directory,xce_logfile,max_queue_jobs,restraints_program):
        QtCore.QThread.__init__(self)
        self.external_software=external_software
        self.initial_model_directory=initial_model_directory
        self.compound_list=compound_list
        self.database_directory=database_directory
        self.data_source_file=data_source_file
        self.todo=todo
        self.ccp4_scratch_directory=ccp4_scratch_directory
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.max_queue_jobs=max_queue_jobs
        self.db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
        self.restraints_program=restraints_program

    def run(self):
        # first remove all ACEDRG input scripts in ccp4_scratch directory
        self.Logfile.insert('removing all xce_acedrg scripts from '+self.ccp4_scratch_directory)
        os.chdir(self.ccp4_scratch_directory)
        os.system('/bin/rm -f xce_acedrg*')

        progress_step=100/float(len(self.compound_list))
        progress=0
        counter=1
        for item in self.compound_list:
            sampleID=item[0]
            compoundID=item[1]
            sm=item[2]
            # if counterions are present, split and take the longest string in list
            smiles = max(sm.split('.'), key=len)
            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating ACEDRG shell script for '+sampleID)
            if compoundID=='' or compoundID is None:
                compoundID='compound'

            if not os.path.isdir(os.path.join(self.initial_model_directory,sampleID)):
                os.mkdir(os.path.join(self.initial_model_directory,sampleID))

            if self.todo=='ALL' or self.todo=='SELECTED':
                # remove symbolic links if present
                if os.path.isfile(os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.pdb')):
                    os.system('/bin/rm '+os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.pdb'))
#               commented this out since people found the presence of the old.cif file interfering with pandda.inspect
#                if os.path.isfile(os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.cif')):
#                    # copy existing CIF file to old.cif so that it can be used as input in
#                    # restraints generating program
#                    os.chdir(os.path.join(self.initial_model_directory,sampleID))
#                    os.system('/bin/cp %s old.cif' %(compoundID.replace(' ','')+'.cif'))
#                    os.system('/bin/rm '+os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.cif'))
                if os.path.isfile(os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.png')):
                    os.system('/bin/rm '+os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.png'))
                if os.path.isdir(os.path.join(self.initial_model_directory,sampleID,'compound')):
                    os.system('/bin/rm -fr '+os.path.join(self.initial_model_directory,sampleID,'compound'))
                    db_dict= {'RefinementCIFStatus': 'pending'}
                    self.Logfile.insert('{0!s}: removed compound directory and all its contents'.format(sampleID))
                    self.Logfile.insert('{0!s}: setting RefinementCIFStatus flag to started'.format(sampleID))
                    self.db.update_data_source(sampleID,db_dict)

            # create 'compound' directory if not present
            if not os.path.isdir(os.path.join(self.initial_model_directory,sampleID,'compound')):
                os.mkdir(os.path.join(self.initial_model_directory,sampleID,'compound'))

            # create text file which contains the smiles string
            if not os.path.isfile(os.path.join(self.initial_model_directory,sampleID,'compound',compoundID+'.smiles')):
                os.chdir(os.path.join(self.initial_model_directory,sampleID,'compound'))
                f=open(compoundID+'.smiles','w')
                f.write(smiles)
                f.close()

            if not os.path.isfile(os.path.join(self.initial_model_directory,sampleID,compoundID.replace(' ','')+'.cif')) \
                    or self.todo=='SELECTED':
                os.chdir(os.path.join(self.initial_model_directory,sampleID))
                os.system('/bin/rm -f %s*' %compoundID.replace(' ',''))
                os.chdir(os.path.join(self.initial_model_directory,sampleID,'compound'))


                helpers().make_png( self.initial_model_directory,
                                    sampleID,
                                    compoundID,
                                    smiles,
                                    self.external_software,
                                    self.database_directory,
                                    self.data_source_file,
                                    self.ccp4_scratch_directory,
                                    counter,
                                    self.xce_logfile,
                                    self.restraints_program )
                counter += 1

                db_dict= {'RefinementCIFprogram': self.restraints_program, 'RefinementCIFStatus': 'started'}
                self.Logfile.insert('{0!s}: setting RefinementCIFStatus flag to started'.format(sampleID))
                self.db.update_data_source(sampleID,db_dict)

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        # submit array job at Diamond
        self.Logfile.insert('created input scripts for '+str(counter)+' ACEDRG jobs in '+self.ccp4_scratch_directory)
        os.chdir(self.ccp4_scratch_directory)
        self.Logfile.insert('changing directory to '+self.ccp4_scratch_directory)
        if counter > 1:
#            if os.getcwd().startswith('/dls'):
            if os.path.isdir('/dls'):
                if self.external_software['qsub_array']:
                    Cmds = (
                            '#PBS -joe -N xce_acedrg_master\n'
                            './xce_acedrg_$SGE_TASK_ID.sh\n'
                            )
                    f = open('acedrg_master.sh','w')
                    f.write(Cmds)
                    f.close()
                    self.Logfile.insert('submitting array job with maximal 100 jobs running on cluster')
                    self.Logfile.insert('using the following command:')
                    self.Logfile.insert('         qsub -P labxchem -q medium.q -t 1:{0!s} -tc {1!s} acedrg_master.sh'.format(str(counter), self.max_queue_jobs))
                    os.system('qsub -P labxchem -q medium.q -t 1:{0!s} -tc {1!s} -N acedrg acedrg_master.sh'.format(str(counter), self.max_queue_jobs))
                else:
                    self.Logfile.insert("cannot start ARRAY job: make sure that 'module load global/cluster' is in your .bashrc or .cshrc file")
            elif self.external_software['qsub']:
                self.Logfile.insert('submitting {0!s} individual jobs to cluster'.format((str(counter))))
                self.Logfile.insert('WARNING: this could potentially lead to a crash...')
                for i in range(counter):
                    self.Logfile.insert('qsub -q medium.q -N acedrg xce_acedrg_{0!s}.sh'.format((str(i+1))))
                    os.system('qsub -N acedrg xce_acedrg_{0!s}.sh'.format((str(i+1))))
            else:
                self.Logfile.insert('running %s consecutive ACEDRG jobs on your local machine')
                for i in range(counter):
                    self.Logfile.insert('starting xce_acedrg_{0!s}.sh'.format((str(i+1))))
                    os.system('./xce_acedrg_{0!s}.sh'.format((str(i+1))))

#        self.emit(QtCore.SIGNAL("finished()"))
        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))

class run_dimple_on_all_autoprocessing_files(QtCore.QThread):
    def __init__(self,sample_list,initial_model_directory,external_software,ccp4_scratch_directory,database_directory,
                 data_source_file,max_queue_jobs,xce_logfile, remote_submission, remote_submission_string, dimple_twin_mode):
        QtCore.QThread.__init__(self)
        self.sample_list=sample_list
        self.initial_model_directory=initial_model_directory
        self.external_software=external_software
        self.queueing_system_available=external_software['qsub']
        self.ccp4_scratch_directory=ccp4_scratch_directory
        self.database_directory=database_directory
        self.data_source_file=data_source_file
        self.max_queue_jobs=max_queue_jobs
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.pipeline='dimple'
        self.using_remote_qsub_submission = remote_submission
        self.remote_qsub_submission = remote_submission_string
        self.dimple_twin_mode = dimple_twin_mode


    def run(self):
        progress_step=1
        if len(self.sample_list) != 0:
            progress_step=100/float(len(self.sample_list))
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        os.chdir(self.ccp4_scratch_directory)
        os.system('/bin/rm -f xce_dimple*sh')

        db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
        database=os.path.join(self.database_directory,self.data_source_file)

        print(self.sample_list)

        for n,item in enumerate(self.sample_list):

            xtal =                  item[0]
            visit_run_autoproc =    item[1]
            mtzin =                 item[2]
            ref_pdb =               item[3]
            ref_mtz =               item[4]
            ref_cif =               item[5]

            # check if reference mtzfile has an Rfree column; if not, then ignore
            # DIMPLE assumes an Rfree column and barfs if it is not present
            # note: ref_mtz looks like this: ref mtz  -R reference.mtz
#            if os.path.isfile(ref_mtz.split()[len(ref_mtz.split())-1]):
#                mtz_column_dict=mtztools(ref_mtz.split()[len(ref_mtz.split())-1]).get_all_columns_as_dict()
            if os.path.isfile(ref_mtz):
                if 'FreeR_flag' not in mtz.object(ref_mtz).column_labels():
#                mtz_column_dict=mtztools(ref_mtz).get_all_columns_as_dict()
#                if 'FreeR_flag' not in mtz_column_dict['RFREE']:
                    self.Logfile.insert('cannot find FreeR_flag in reference mtz file: {0!s} -> ignoring reference mtzfile!!!'.format(ref_mtz))
                    ref_mtz = ''
#                    if mtz_column_dict['RFREE'] != []:
#                        self.Logfile.insert('found Rfree set with other column name though: {0!s}'.format(str(mtz_column_dict['RFREE'])))
#                        self.Logfile.insert('try renaming Rfree column to FreeR_flag with CAD!')

#            uniqueify = ''
#            if 'FreeR_flag' in mtz.object(mtzin).column_labels():
#                uniqueify = 'uniqueify -f FreeR_flag ' + mtzin + ' ' + xtal + '-unique.mtz' + '\n'
#            else:
#                uniqueify = 'uniqueify ' + mtzin + ' ' + xtal + '-unique.mtz' + '\n'

            db_dict= {'DimpleReferencePDB': ref_pdb}
            db.update_data_source(xtal,db_dict)

            self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating input script for '+xtal+' in '+visit_run_autoproc)

            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
                os.mkdir(os.path.join(self.initial_model_directory,xtal))
            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple')):
                os.mkdir(os.path.join(self.initial_model_directory,xtal,'dimple'))
            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc)):
                os.mkdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
            os.chdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
            os.system('touch dimple_run_in_progress')
            os.system('/bin/rm final.mtz 2> /dev/null')
            os.system('/bin/rm final.pdb 2> /dev/null')

            if self.queueing_system_available:
                top_line='#PBS -joe -N XCE_dimple\n'
            else:
                top_line='#!'+os.getenv('SHELL')+'\n'

            if 'csh' in os.getenv('SHELL'):
                ccp4_scratch='setenv CCP4_SCR '+self.ccp4_scratch_directory+'\n'
            elif 'bash' in os.getenv('SHELL'):
                ccp4_scratch='export CCP4_SCR='+self.ccp4_scratch_directory+'\n'
            else:
                ccp4_scratch=''

            if 'dimple_rerun_on_selected_file' in visit_run_autoproc:
                additional_cmds = (
                            'cd {0!s}\n'.format(os.path.join(self.initial_model_directory,xtal)) +
                            '/bin/rm dimple.pdb\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/final.pdb dimple.pdb\n'
                            '/bin/rm dimple.mtz\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/final.mtz dimple.mtz\n'
                            '/bin/rm 2fofc.map\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/2fofc.map .\n'
                            '/bin/rm fofc.map\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/fofc.map .\n'
                            '\n'
                            '$CCP4/bin/ccp4-python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','update_data_source_for_new_dimple_pdb.py')+
                            ' {0!s} {1!s} {2!s}\n'.format(os.path.join(self.database_directory,self.data_source_file), xtal, self.initial_model_directory)  )

            else:
                additional_cmds=''

#            resolution_high = 0.1
#            o = iotbx.mtz.object(mtzin)
#            low, high = o.max_min_resolution()
            hkl = any_reflection_file(file_name=mtzin)
            miller_arrays = hkl.as_miller_arrays()
            mtzFile = miller_arrays[0]

            if mtzFile.space_group_info().symbol_and_number() ==  'R 3 :H (No. 146)':
                symNoAbsence = 'H3'
            else:
                symNoAbsence = str([x[0] for x in str(mtzFile.space_group_info().symbol_and_number().split('(')[0]).split()]).replace('[','').replace(']','').replace("'","").replace(',','').replace(' ','')

            dls_stuff = ''
            if os.path.isdir('/dls'):
                dls_stuff = (
                    'module unload ccp4\n'
                    'source /dls/science/groups/i04-1/software/pandda_0.2.12/ccp4/ccp4-7.0/bin/ccp4.setup-sh\n'
                )

            twin = ''
            if self.dimple_twin_mode:
                twin = "--refmac-key 'TWIN'"

            Cmds = (
                    '{0!s}\n'.format(top_line)+
                    '\n'
                    'export XChemExplorer_DIR="'+os.getenv('XChemExplorer_DIR')+'"\n'
                    '\n'
                    'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc) +
                    '\n'
                    'source $XChemExplorer_DIR/setup-scripts/xce.setup-sh\n'
                    '\n'
                    + dls_stuff +
                    '\n'
                    + ccp4_scratch +
                    '\n'
                    '$CCP4/bin/ccp4-python $XChemExplorer_DIR/helpers/update_status_flag.py %s %s %s %s\n' %(database,xtal,'DimpleStatus','running') +
                    '\n'
                    'unique hklout unique.mtz << eof\n'
                    ' cell %s\n' %str([round(float(i),2) for i in mtzFile.unit_cell().parameters()]).replace('[','').replace(']','')+
                    ' symmetry %s\n' %symNoAbsence+
                    ' resolution %s\n' %str(round(float(mtzFile.d_min()),3))+
                    'eof\n'
                    '\n'
#                    'freerflag hklin unique.mtz hklout free.mtz > freerflag.log\n'
                    '\n'
                    'sftools << eof > sftools.log\n'
                    ' read unique.mtz\n'
                    ' calc col F = 10.0\n'
                    ' calc col SIGF = 1.0\n'
                    ' write sftools.mtz\n'
                    'eof\n'
                    '\n'
                    'cad hklin1 sftools.mtz hklin2 %s hklout %s.999A.mtz << eof\n' %(mtzin,xtal) +
                    ' monitor BRIEF\n'
                    ' labin file 1 E1=F E2=SIGF\n'
                    ' labout file 1 E1=F_unique E2=SIGF_unique\n'
                    ' labin file 2 ALL\n'
                    ' resolution file 1 999.0 %s\n' %str(round(float(mtzFile.d_min()),2))+
                    'eof\n'
                    '\n'
                    "dimple --no-cleanup %s.999A.mtz %s %s %s %s dimple\n" %(xtal,ref_pdb,ref_mtz,ref_cif,twin) +
                    '\n'
                    'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc,'dimple') +
#                    + uniqueify +
#                    '\n'
#                    'cad hklin1 %s hklout %s <<eof\n' %(str(xtal + '-unique.mtz'), str(xtal + '-unique-2.mtz')) +
#                    'monitor BRIEF\n'
#                    'labin file 1 -\n' 
#                    '    ALL\n'
#                    'resolution file 1 999.0 %s\n' %(high) +
#                    'eof\n'
#                    'dimple --no-cleanup %s-unique-2.mtz %s %s %s dimple\n' %(xtal,ref_pdb,ref_mtz,ref_cif) +
#                    '\n'
#                    'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc,'dimple') +
                    '\n'
                    'fft hklin final.mtz mapout 2fofc.map << EOF\n'
                    ' labin F1=FWT PHI=PHWT\n'
                    'EOF\n'
                    '\n'
                    'fft hklin final.mtz mapout fofc.map << EOF\n'
                    ' labin F1=DELFWT PHI=PHDELWT\n'
                    'EOF\n'
                    '\n'
                    +additional_cmds+
                    '\n'
                    'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc) +
                    '\n'
                    '/bin/rm dimple_run_in_progress\n'
                    '\n'
                    'ln -s dimple/final.pdb .\n'
                    'ln -s dimple/final.mtz .\n'
                    )

            # print(Cmds)

            os.chdir(self.ccp4_scratch_directory)
            f = open('xce_dimple_{0!s}.sh'.format(str(n+1)),'w')
            f.write(Cmds)
            f.close()
            os.system('chmod +x xce_dimple_{0!s}.sh'.format(str(n+1)))
            db_dict= {'DimpleStatus': 'started'}
            self.Logfile.insert('{0!s}: setting DataProcessingStatus flag to started'.format(xtal))
            db.update_data_source(xtal,db_dict)


            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        # submit job
        self.Logfile.insert('created input scripts for '+str(n+1)+' in '+self.ccp4_scratch_directory)
        os.chdir(self.ccp4_scratch_directory)
        if os.path.isdir('/dls'):
            if self.external_software['qsub_array']:
                Cmds = (
                        '#PBS -joe -N xce_dimple_master\n'
                        './xce_dimple_$SGE_TASK_ID.sh\n'
                        )
                f = open('dimple_master.sh','w')
                f.write(Cmds)
                f.close()
                print(os.getcwd())
                self.Logfile.insert('submitting array job with maximal 100 jobs running on cluster')
                self.Logfile.insert('using the following command:')
                self.Logfile.insert('qsub -P labxchem -q medium.q -t 1:{0!s} -tc {1!s} dimple_master.sh'.format(str(n+1), self.max_queue_jobs))
                if self.using_remote_qsub_submission:
                    os.system(str(self.remote_qsub_submission).replace("qsub'", str('cd ' + str(os.getcwd()) + '; ' + 'qsub -P labxchem -q medium.q -t 1:{0!s} -tc {1!s} dimple_master.sh'
                                                                  .format(str(n+1), self.max_queue_jobs))) + "'")

                else:
                    os.system('qsub -P labxchem -t 1:{0!s} -tc {1!s} -N dimple-master dimple_master.sh'.format(str(n+1), self.max_queue_jobs))
            else:
                self.Logfile.insert("cannot start ARRAY job: make sure that 'module load global/cluster' is in your .bashrc or .cshrc file")
        elif self.external_software['qsub']:
            self.Logfile.insert('submitting {0!s} individual jobs to cluster'.format((str(n+1))))
            self.Logfile.insert('WARNING: this could potentially lead to a crash...')
            for i in range(n+1):
                self.Logfile.insert('qsub -q medium.q -N dimple xce_dimple_{0!s}.sh'.format((str(i+1))))
                os.system('qsub -N dimple xce_dimple_{0!s}.sh'.format((str(i+1))))
        else:
            self.Logfile.insert('running %s consecutive DIMPLE jobs on your local machine')
            for i in range(n+1):
                self.Logfile.insert('starting xce_dimple_{0!s}.sh'.format((str(i+1))))
                os.system('./xce_dimple_{0!s}.sh'.format((str(i+1))))

        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))

class run_dimple_on_all_autoprocessing_files_new(QtCore.QThread):
    def __init__(self,sample_list,initial_model_directory,external_software,ccp4_scratch_directory,database_directory,data_source_file,max_queue_jobs,xce_logfile):
        QtCore.QThread.__init__(self)
        self.sample_list=sample_list
        self.initial_model_directory=initial_model_directory
        self.external_software=external_software
        self.queueing_system_available=external_software['qsub']
        self.ccp4_scratch_directory=ccp4_scratch_directory
        self.database_directory=database_directory
        self.data_source_file=data_source_file

        self.db=XChemDB.data_source(os.path.join(self.database_directory,self.data_source_file))
        self.database=os.path.join(self.database_directory,self.data_source_file)

        self.max_queue_jobs=max_queue_jobs
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)

        self.n=0

        self.pipeline='dimple'


    def run(self):
        progress_step=1
        if len(self.sample_list) != 0:
            progress_step=100/float(len(self.sample_list))
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        os.chdir(self.ccp4_scratch_directory)
        os.system('/bin/rm -f xce_{0!s}*sh'.format(self.pipeline))


        for item in sorted(self.sample_list):

            xtal =                  item[0]
            visit_run_autoproc =    item[1]
            mtzin =                 item[2]
            ref_pdb =               item[3]
            ref_mtz =               item[4]
            ref_cif =               item[5]

            if 'dimple_rerun_on_selected_file' in visit_run_autoproc:
                if self.pipeline=='dimple':
                    self.prepare_dimple_shell_script(xtal,visit_run_autoproc,mtzin,ref_pdb,ref_mtz,ref_cif)
                elif self.pipeline=='pipedream':
                    self.prepare_pipedream_shell_script()
                elif self.pipeline=='phenix.ligand_pipeline':
                    self.prepare_phenix_ligand_pipeline_shell_script()
            else:
                self.prepare_dimple_shell_script(xtal,visit_run_autoproc,mtzin,ref_pdb,ref_mtz,ref_cif)

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.run_script()

        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))






    def prepare_phenix_ligand_pipeline_shell_script(self,xtal,visit_run_autoproc,mtzin,ref_pdb,ref_mtz,ref_cif):

        # check if reference mtzfile has an Rfree column; if not, then ignore
        # DIMPLE assumes an Rfree column and barfs if it is not present
        # note: ref_mtz looks like this: ref mtz  -R reference.mtz
#        if os.path.isfile(ref_mtz):
#            mtz_column_dict=mtztools(ref_mtz).get_all_columns_as_dict()
#            if 'FreeR_flag' not in mtz_column_dict['RFREE']:
#                self.Logfile.insert('cannot find FreeR_flag in reference mtz file: %s -> ignoring reference mtzfile!!!' %ref_mtz)
#                ref_mtz = ''
#                if mtz_column_dict['RFREE'] != []:
#                    self.Logfile.insert('found Rfree set with other column name though: %s' %str(mtz_column_dict['RFREE']))
#                    self.Logfile.insert('try renaming Rfree column to FreeR_flag with CAD!')
#
#        db_dict={}
#        db_dict['DimpleReferencePDB']=ref_pdb
#        self.db.update_data_source(xtal,db_dict)
#
#        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating input script for '+xtal+' in '+visit_run_autoproc)
#

        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
            os.mkdir(os.path.join(self.initial_model_directory,xtal))
        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'phenix.ligand_pipeline')):
            os.mkdir(os.path.join(self.initial_model_directory,xtal,'phenix.ligand_pipeline'))
        os.chdir(os.path.join(self.initial_model_directory,xtal,'phenix.ligand_pipeline'))
        os.system('touch dimple_run_in_progress')
        os.system('/bin/rm final.mtz 2> /dev/null')
        os.system('/bin/rm final.pdb 2> /dev/null')

        if self.queueing_system_available:
            top_line='#PBS -joe -N XCE_{0!s}\n'.format(self.pipeline)
        else:
            top_line='#!'+os.getenv('SHELL')+'\n'

        if 'csh' in os.getenv('SHELL'):
            ccp4_scratch='setenv CCP4_SCR '+self.ccp4_scratch_directory+'\n'
        elif 'bash' in os.getenv('SHELL'):
            ccp4_scratch='export CCP4_SCR='+self.ccp4_scratch_directory+'\n'
        else:
            ccp4_scratch=''

        if os.path.isdir('/dls'):
            ccp4_scratch+='module load phenix\n'

        Cmds = (
                '{0!s}\n'.format(top_line)+
                '\n'
                'export XChemExplorer_DIR="'+os.getenv('XChemExplorer_DIR')+'"\n'
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'phenix.ligand_pipeline') +
                '\n'
                'source $XChemExplorer_DIR/setup-scripts/xce.setup-sh\n'
                '\n'
                'module unload ccp4\n'
                'source /dls/science/groups/i04-1/software/pandda_0.2.12/ccp4/ccp4-7.0/bin/ccp4.setup-sh\n'
                +ccp4_scratch+
                '\n'
                '$CCP4/bin/ccp4-python $XChemExplorer_DIR/helpers/update_status_flag.py %s %s %s %s\n' %(database,xtal,'DimpleStatus','running') +
                '\n'
                'phenix.ligand_pipeline %s %s' %(ref_pdb,mtzin)+
                ' mr=False'
                ' ligand_copies=0'
                ' build=False'
                ' prune=False'
                ' remove_waters=False'
                ' stop_if_r_free_greater_than=0.4'
                ' update_waters=False'
                ' build_hydrogens=False\n'
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal) +
                '\n'
                'ln -s phenix.ligand_pipeline/pipeline_1/refine_final.pdb dimple.pdb'
                'ln -s phenix.ligand_pipeline/pipeline_1/refine_final.mtz dimple.mtz'
                '\n'
                'fft hklin dimple.mtz mapout 2fofc.map << EOF\n'
                ' labin F1=2FOFCWT PHI=PH2FOFCWT\n'
                'EOF\n'
                '\n'
                'fft hklin dimple.mtz mapout fofc.map << EOF\n'
                ' labin F1=FOFCWT PHI=PHFOFCWT\n'
                'EOF\n'
                '\n'
                '$CCP4/bin/ccp4-python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','update_data_source_for_new_dimple_pdb.py')+
                ' {0!s} {1!s} {2!s}\n'.format(os.path.join(self.database_directory,self.data_source_file), xtal, self.initial_model_directory)+
                '\n'
                '/bin/rm dimple_run_in_progress\n'
                )

        os.chdir(self.ccp4_scratch_directory)
        f = open('xce_{0!s}_{1!s}.sh'.format(self.pipeline, str(n+1)),'w')
        f.write(Cmds)
        f.close()
        self.n+=1
        os.system('chmod +x xce_{0!s}_{1!s}.sh'.format(self.pipeline, str(n+1)))
        db_dict= {'DimpleStatus': 'started'}
        self.Logfile.insert('{0!s}: setting DataProcessingStatus flag to started'.format(xtal))
        self.db.update_data_source(xtal,db_dict)


    def prepare_pipedream_shell_script(self,xtal,visit_run_autoproc,mtzin,ref_pdb,ref_mtz,ref_cif):

        # check if reference mtzfile has an Rfree column; if not, then ignore
        # DIMPLE assumes an Rfree column and barfs if it is not present
        # note: ref_mtz looks like this: ref mtz  -R reference.mtz
#        if os.path.isfile(ref_mtz):
#            mtz_column_dict=mtztools(ref_mtz).get_all_columns_as_dict()
#            if 'FreeR_flag' not in mtz_column_dict['RFREE']:
#                self.Logfile.insert('cannot find FreeR_flag in reference mtz file: %s -> ignoring reference mtzfile!!!' %ref_mtz)
#                ref_mtz = ''
#                if mtz_column_dict['RFREE'] != []:
#                    self.Logfile.insert('found Rfree set with other column name though: %s' %str(mtz_column_dict['RFREE']))
#                    self.Logfile.insert('try renaming Rfree column to FreeR_flag with CAD!')
#
#        db_dict={}
#        db_dict['DimpleReferencePDB']=ref_pdb
#        self.db.update_data_source(xtal,db_dict)
#
#        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating input script for '+xtal+' in '+visit_run_autoproc)
#

        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
            os.mkdir(os.path.join(self.initial_model_directory,xtal))
        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'pipedream')):
            os.mkdir(os.path.join(self.initial_model_directory,xtal,'pipedream'))
        os.chdir(os.path.join(self.initial_model_directory,xtal,'pipedream'))
        os.system('touch dimple_run_in_progress')
        os.system('/bin/rm final.mtz 2> /dev/null')
        os.system('/bin/rm final.pdb 2> /dev/null')

        if self.queueing_system_available:
            top_line='#PBS -joe -N XCE_{0!s}\n'.format(self.pipeline)
        else:
            top_line='#!'+os.getenv('SHELL')+'\n'

        if 'csh' in os.getenv('SHELL'):
            ccp4_scratch='setenv CCP4_SCR '+self.ccp4_scratch_directory+'\n'
        elif 'bash' in os.getenv('SHELL'):
            ccp4_scratch='export CCP4_SCR='+self.ccp4_scratch_directory+'\n'
        else:
            ccp4_scratch=''

        if os.path.isdir('/dls'):
            ccp4_scratch+='module load buster\n'

        if os.path.isfile(ref_mtz):
            hklref_line=' -hklref {0!s}'.format(ref_mtz)
        else:
            hklref_line=''

        Cmds = (
                '{0!s}\n'.format(top_line)+
                '\n'
                'export XChemExplorer_DIR="'+os.getenv('XChemExplorer_DIR')+'"\n'
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'pipedream') +
                '\n'
                'source $XChemExplorer_DIR/setup-scripts/xce.setup-sh\n'
                '\n'
                +ccp4_scratch+
                '\n'
                '$CCP4/bin/ccp4-python $XChemExplorer_DIR/helpers/update_status_flag.py %s %s %s %s\n' %(database,xtal,'DimpleStatus','running') +
                '\n'
                'pipedream '
                ' -d pipedreamDir'
                ' -xyzin %s' %ref_pdb+
                hklref_line+
                ' -hklin {0!s}'.format(mtzin)+
                ' -keepwater\n'
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal) +
                '\n'
                'ln -s phenix.ligand_pipeline/pipeline_1/refine_final.pdb dimple.pdb'
                'ln -s phenix.ligand_pipeline/pipeline_1/refine_final.mtz dimple.mtz'
                '\n'
                'fft hklin dimple.mtz mapout 2fofc.map << EOF\n'
                ' labin F1=FWT PHI=PHWT\n'
                'EOF\n'
                '\n'
                'fft hklin dimple.mtz mapout fofc.map << EOF\n'
                ' labin F1=DELFWT PHI=PHDELWT\n'
                'EOF\n'
                '\n'
                '$CCP4/libexec/python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','update_data_source_for_new_dimple_pdb.py')+
                ' {0!s} {1!s} {2!s}\n'.format(os.path.join(self.database_directory,self.data_source_file), xtal, self.initial_model_directory)+
                '\n'
                '/bin/rm dimple_run_in_progress\n'
                )

        os.chdir(self.ccp4_scratch_directory)
        f = open('xce_{0!s}_{1!s}.sh'.format(self.pipeline, str(n+1)),'w')
        f.write(Cmds)
        f.close()
        self.n+=1
        os.system('chmod +x xce_{0!s}_{1!s}.sh'.format(self.pipeline, str(n+1)))
        db_dict= {'DimpleStatus': 'started'}
        self.Logfile.insert('{0!s}: setting DataProcessingStatus flag to started'.format(xtal))
        self.db.update_data_source(xtal,db_dict)




#pipedream -d pipedreamDir -xyzin $reference.pdb -hklref $referenceDataset.mtz -hklin $mtzFileUnique.mtz -rhofit $inhib.cif  -target $reference.pdb -keepwater -nthreads 12 -postquick -allclusters |tee -a pipedream.log


    def prepare_dimple_shell_script(self,xtal,visit_run_autoproc,mtzin,ref_pdb,ref_mtz,ref_cif):

        # check if reference mtzfile has an Rfree column; if not, then ignore
        # DIMPLE assumes an Rfree column and barfs if it is not present
        # note: ref_mtz looks like this: ref mtz  -R reference.mtz
        if os.path.isfile(ref_mtz):
            mtz_column_dict=mtztools(ref_mtz).get_all_columns_as_dict()
            if 'FreeR_flag' not in mtz_column_dict['RFREE']:
                self.Logfile.insert('cannot find FreeR_flag in reference mtz file: {0!s} -> ignoring reference mtzfile!!!'.format(ref_mtz))
                ref_mtz = ''
                if mtz_column_dict['RFREE']:
                    self.Logfile.insert('found Rfree set with other column name though: {0!s}'.format(str(mtz_column_dict['RFREE'])))
                    self.Logfile.insert('try renaming Rfree column to FreeR_flag with CAD!')

        db_dict= {'DimpleReferencePDB': ref_pdb}
        self.db.update_data_source(xtal,db_dict)

        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'creating input script for '+xtal+' in '+visit_run_autoproc)

        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
            os.mkdir(os.path.join(self.initial_model_directory,xtal))
        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple')):
            os.mkdir(os.path.join(self.initial_model_directory,xtal,'dimple'))
        if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc)):
            os.mkdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
        os.chdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
        os.system('touch dimple_run_in_progress')
        os.system('/bin/rm final.mtz 2> /dev/null')
        os.system('/bin/rm final.pdb 2> /dev/null')

        if self.queueing_system_available:
            top_line='#PBS -joe -N XCE_dimple\n'
        else:
            top_line='#!'+os.getenv('SHELL')+'\n'

        if 'csh' in os.getenv('SHELL'):
            ccp4_scratch='setenv CCP4_SCR '+self.ccp4_scratch_directory+'\n'
        elif 'bash' in os.getenv('SHELL'):
            ccp4_scratch='export CCP4_SCR='+self.ccp4_scratch_directory+'\n'
        else:
            ccp4_scratch=''

        if 'dimple_rerun_on_selected_file' in visit_run_autoproc:
            additional_cmds = (
                            'cd {0!s}\n'.format(os.path.join(self.initial_model_directory,xtal)) +
                            '/bin/rm dimple.pdb\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/final.pdb dimple.pdb\n'
                            '/bin/rm dimple.mtz\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/final.mtz dimple.mtz\n'
                            '/bin/rm 2fofc.map\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/2fofc.map .\n'
                            '/bin/rm fofc.map\n'
                            'ln -s dimple/dimple_rerun_on_selected_file/dimple/fofc.map .\n'
                            '\n'
                            '$CCP4/libexec/python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','update_data_source_for_new_dimple_pdb.py')+
                            ' {0!s} {1!s} {2!s}\n'.format(os.path.join(self.database_directory,self.data_source_file), xtal, self.initial_model_directory)  )

        else:
            additional_cmds=''

        Cmds = (
                '{0!s}\n'.format(top_line)+
                '\n'
                'export XChemExplorer_DIR="'+os.getenv('XChemExplorer_DIR')+'"\n'
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc) +
                '\n'
                'source $XChemExplorer_DIR/setup-scripts/xce.setup-sh\n'
                '\n'
                +ccp4_scratch+
                '\n'
                '$CCP4/bin/ccp4-python $XChemExplorer_DIR/helpers/update_status_flag.py %s %s %s %s\n' %(database,xtal,'DimpleStatus','running') +
                '\n'
                'dimple --no-cleanup %s %s %s %s dimple\n' %(mtzin,ref_pdb,ref_mtz,ref_cif) +
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc,'dimple') +
                '\n'
                'fft hklin final.mtz mapout 2fofc.map << EOF\n'
                ' labin F1=FWT PHI=PHWT\n'
                'EOF\n'
                '\n'
                'fft hklin final.mtz mapout fofc.map << EOF\n'
                ' labin F1=DELFWT PHI=PHDELWT\n'
                'EOF\n'
                '\n'
                +additional_cmds+
                '\n'
                'cd %s\n' %os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc) +
                '\n'
                '/bin/rm dimple_run_in_progress\n'
                '\n'
                'ln -s dimple/final.pdb .\n'
                'ln -s dimple/final.mtz .\n'
                )

        os.chdir(self.ccp4_scratch_directory)
        f = open('xce_{0!s}_{1!s}.sh'.format(self.pipeline, str(n+1)),'w')
        f.write(Cmds)
        f.close()
        self.n+=1
        os.system('chmod +x xce_{0!s}_{1!s}.sh'.format(self.pipeline, str(n+1)))
        db_dict= {'DimpleStatus': 'started'}
        self.Logfile.insert('{0!s}: setting DataProcessingStatus flag to started'.format(xtal))
        self.db.update_data_source(xtal,db_dict)


    def run_script(self):
        # submit job
        self.Logfile.insert('created input scripts for '+str(self.n)+' in '+self.ccp4_scratch_directory)
        os.chdir(self.ccp4_scratch_directory)
        if os.path.isdir('/dls'):
            if self.external_software['qsub_array']:
                Cmds = (
                        '#PBS -joe -N xce_{0!s}_master\n'.format(self.pipeline)+
                        './xce_{0!s}_$SGE_TASK_ID.sh\n'.format(self.pipeline)
                        )
                f = open('{0!s}_master.sh'.format(self.pipeline),'w')
                f.write(Cmds)
                f.close()
                self.Logfile.insert('submitting array job with maximal 100 jobs running on cluster')
                self.Logfile.insert('using the following command:')
                self.Logfile.insert('qsub -P labxchem -q medium.q -t 1:{0!s} -tc {1!s} {2!s}_master.sh'.format(str(self.n), self.max_queue_jobs, self.pipeline))
                os.system('qsub -P labxchem -q medium.q -t 1:{0!s} -tc {1!s} {2!s}_master.sh'.format(str(self.n), self.max_queue_jobs, self.pipeline))
            else:
                self.Logfile.insert("cannot start ARRAY job: make sure that 'module load global/cluster' is in your .bashrc or .cshrc file")
        elif self.external_software['qsub']:
            self.Logfile.insert('submitting {0!s} individual jobs to cluster'.format((str(self.n))))
            self.Logfile.insert('WARNING: this could potentially lead to a crash...')
            for i in range(self.n):
                self.Logfile.insert('qsub -q medium.q xce_{0!s}_{1!s}.sh'.format(str(i+1), self.pipeline))
                os.system('qsub -q medium.q xce_{0!s}_{1!s}.sh'.format(str(i+1), self.pipeline))
        else:
            self.Logfile.insert('running {0!s} consecutive {1!s} jobs on your local machine'.format(*(self.pipeline)))
            for i in range(self.n):
                self.Logfile.insert('starting xce_{0!s}_{1!s}.sh'.format(str(i+1), self.pipeline))
                os.system('./xce_{0!s}_{1!s}.sh'.format(str(i+1), self.pipeline))



class remove_selected_dimple_files(QtCore.QThread):
    def __init__(self,sample_list,initial_model_directory,xce_logfile,database_directory,data_source_file,):
        QtCore.QThread.__init__(self)
        self.sample_list=sample_list
        self.initial_model_directory=initial_model_directory
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.db=XChemDB.data_source(os.path.join(database_directory,data_source_file))

    def run(self):
        progress_step=1
        if len(self.sample_list) != 0:
            progress_step=100/float(len(self.sample_list))
        progress=0
        self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        for n,xtal in enumerate(self.sample_list):
            db_dict={}
            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
                self.Logfile.insert('{0!s}: directory does not exist'.format(xtal))
                continue
            os.chdir(os.path.join(self.initial_model_directory,xtal))
            self.Logfile.insert('{0!s}: removing dimple.pdb/dimple.mtz'.format(xtal))
            os.system('/bin/rm dimple.pdb 2> /dev/null')
            os.system('/bin/rm dimple.mtz 2> /dev/null')
            if os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple','dimple_rerun_on_selected_file')):
                os.chdir('dimple')
                self.Logfile.insert('{0!s} removing directory dimple/dimple_rerun_on_selected_file'.format(xtal))
                os.system('/bin/rm -fr dimple_rerun_on_selected_file')

            db_dict['DimpleResolutionHigh']=''
            db_dict['DimpleRcryst']=''
            db_dict['DimpleRfree']=''
            db_dict['DimplePathToPDB']=''
            db_dict['DimplePathToMTZ']=''
            db_dict['DimpleReferencePDB']=''
            db_dict['DimplePANDDAwasRun']='False'
            db_dict['DimplePANDDAhit']='False'
            db_dict['DimplePANDDAreject']='False'
            db_dict['DimplePANDDApath']=''
            db_dict['DimpleStatus']='pending'

            self.Logfile.insert('{0!s}: updating database'.format(xtal))
            self.db.update_data_source(xtal,db_dict)


            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.emit(QtCore.SIGNAL('datasource_menu_reload_samples'))




class start_COOT(QtCore.QThread):

    def __init__(self,settings,interface):
        QtCore.QThread.__init__(self)
        self.settings=settings
        if interface=='old':
            self.pylib='XChemCoot.py'
        elif interface=='new':
            self.pylib='XChemCootNew.py'
        elif interface=='panddaV1':
            self.pylib='XChemCootOld.py'
        elif interface=='reference':
            self.pylib='XChemCootReference.py'

    def run(self):
        cwd=os.getcwd()
        # coot at Diamond always or sometimes at least open in home directory, so then it won't find the .pkl file
        pickle.dump(self.settings,open(os.path.join(os.getenv('HOME'),'.xce_settings.pkl'),'wb'))
        os.system('cd {0!s}\ncoot --no-guano --no-state-script --script {1!s}'.format(os.getenv('HOME'), os.path.join(os.getenv('XChemExplorer_DIR'),'lib',self.pylib)))



class start_ICM(QtCore.QThread):

    def __init__(self,html_export_directory):
        QtCore.QThread.__init__(self)
        self.html_export_directory=html_export_directory

    def run(self):
        cwd=os.getcwd()
        if cwd.startswith('/dls'):
#            if not os.path.isfile(os.path.join('/home',getpass.getuser(),'.flexlmrc')):
#                os.system('touch '+os.path.join('/home',getpass.getuser(),'.flexlmrc'))
#                f=open(os.path.join('/home',getpass.getuser(),'.flexlmrc'),'w')
#                f.write('MOLSOFTD_LICENSE_FILE=@diamvicmpro.diamond.ac.uk')
#                f.close()
            os.system('nautilus {0!s} &'.format(self.html_export_directory))
            os.system('/dls/science/groups/i04-1/software/icm-3.8-5/icm64 -g')

class start_pandda_inspect(QtCore.QThread):

    def __init__(self,settings,xce_logfile):
        QtCore.QThread.__init__(self)
#        self.settings=settings
        self.panddas_directory=settings['panddas_directory']
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)

    def run(self):
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
                'cd '+self.panddas_directory+'\n'
                'pandda.inspect\n'
            )

        self.Logfile.insert('starting pandda.inspect with the following command:\n'+Cmds)
        os.system(Cmds)

#        Cmds = (
#                'source '+os.path.join(os.getenv('XChemExplorer_DIR'),'setup-scripts','pandda.setup-sh')+'\n'
#                'cd '+self.panddas_directory+'\n'
#                'pandda.inspect\n'
#            )
#
#        self.Logfile.insert('starting pandda.inspect with the following command:\n/bin/bash\n'+Cmds)
#        # need to do this because we're having all sorts of csh bash issues at SGC
#        os.system('/bin/bash\n'+Cmds)





class start_dials_image_viewer(QtCore.QThread):

    def __init__(self,diffraction_image):
        QtCore.QThread.__init__(self)
        self.diffraction_image=diffraction_image

    def run(self):
        os.system('dials.image_viewer '+self.diffraction_image)




#
# --- new module from hell -------------------------------------------------------------------------------------------
#

class read_pinIDs_from_gda_logs(QtCore.QThread):
    def __init__(self,
                 beamline,
                 visit,
                 database,
                 gdaLogInstructions,
                 xce_logfile):
        QtCore.QThread.__init__(self)
        self.beamline = beamline
        self.visit = visit

        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)

        self.db = XChemDB.data_source(os.path.join(database))
        self.allSamples = self.db.collected_xtals_during_visit_for_scoring(visit)

        self.gdaLogInstructions = gdaLogInstructions
        self.gda_log_start_line = gdaLogInstructions[0]
        self.gzipped_logs_parsed = gdaLogInstructions[1]

    def run(self):
        if self.gzipped_logs_parsed:
            self.Logfile.insert("parsed gzipped gda logfiles before, won't do again...")
        self.Logfile.insert("will start parsing of current gda logfile at line {0!s}".format(str(self.gda_log_start_line)))
        self.emit(QtCore.SIGNAL('update_status_bar(QString)'), 'checking GDA logiles for pinID details')
        pinDict, self.gda_log_start_line = XChemMain.get_gda_barcodes(self.allSamples,
                                                                 self.gzipped_logs_parsed,
                                                                 self.gda_log_start_line,
                                                                 self.beamline,
                                                                 self.xce_logfile)

        self.update_database(pinDict)

        self.gdaLogInstructions = [self.gda_log_start_line,True]
        self.Logfile.insert('====== finished checking GDA logfiles ======')
        self.emit(QtCore.SIGNAL('update_gdaLog_parsing_instructions_and_score'), self.gdaLogInstructions)
        self.emit(QtCore.SIGNAL("finished()"))

    def update_database(self,pinDict):
        self.Logfile.insert('updating database with pinDIs from GDA logfiles')
        for sample in pinDict:
            dbDict = {'DataCollectionPinBarcode': pinDict[sample]}
            self.db.update_data_source(sample,dbDict)



class choose_autoprocessing_outcome(QtCore.QThread):
    def __init__(self,
                database,
                visit,
                reference_file_list,
                preferences,
                projectDir,
                rescore,
                xce_logfile):
        QtCore.QThread.__init__(self)
        self.visit = visit
        self.projectDir = projectDir
        self.reference_file_list = reference_file_list
        self.rescore = rescore
        self.selection_mechanism = preferences['dataset_selection_mechanism']
        self.acceptable_unitcell_volume_difference = preferences['allowed_unitcell_difference_percent']
        self.acceptable_low_resolution_limit_for_data = preferences['acceptable_low_resolution_limit_for_data']
        self.acceptable_low_resolution_Rmerge = preferences['acceptable_low_resolution_Rmerge']

        self.xce_logfile = xce_logfile
        self.Logfile = XChemLog.updateLog(xce_logfile)

        self.db = XChemDB.data_source(os.path.join(database))
        self.allSamples = self.db.collected_xtals_during_visit_for_scoring(visit)
#        self.allSamples = self.db.collected_xtals_during_visit_for_scoring(visit,rescore)


    def run(self):
        for sample in sorted(self.allSamples):
            if self.db.autoprocessing_result_user_assigned(sample) and not self.rescore:
                self.Logfile.warning('{0!s}: user has manually selected auto-processing result; will NOT auto-select!'.format(sample))
                continue
            elif self.rescore:
                self.Logfile.warning('{0!s}: rescore selected -> might overwrite user selection!'.format(sample))
            else:
                self.Logfile.insert('%s: selecting autoprocessing result' % sample)
            dbList = self.db.all_autoprocessing_results_for_xtal_as_dict(sample)
            self.Logfile.insert('%s: found %s different autoprocessing results' %(sample,str(len(dbList))))

            # 1.) if posssible, only carry forward samples with similar UCvolume and same point group
            dbList = self.selectResultsSimilarToReference(dbList)

            # 2.) if possible, only carry forward samples with low resolution Rmerge smaller than
            #     specified in perferences
            dbList = self.selectResultsWithAcceptableLowResoRmerge(dbList)

            # 3.) Make selection based on speified selection mechanism
            if self.selection_mechanism == 'IsigI*Comp*UniqueRefl':
                dbDict = self.selectHighestScore(dbList)
            elif self.selection_mechanism == 'highest_resolution':
                dbDict = self.selectHighestResolution(dbList)
            elif 'only' in self.selection_mechanism:
                dbDict = self.selectSpecificPipelineOnly(dbList)

            # 4.) Set new symbolic links in project directory
            XChemMain.linkAutoProcessingResult(sample,dbDict,self.projectDir,self.xce_logfile)

            # 5.) Determine DataProcessing Outcome
            dbDict['DataCollectionOutcome'] = self.determine_processing_outcome(dbDict)

            # 6.) update database
            dbDict['DataProcessingAutoAssigned'] = 'True'
            self.updateDB(sample,dbDict)

        self.emit(QtCore.SIGNAL('populate_datasets_summary_table_NEW'))
        self.emit(QtCore.SIGNAL("finished()"))

    def report_forward_carried_pipelines(self,dbListOut,dbList):
        if not dbListOut:
            dbListOut = dbList
            self.Logfile.warning('none of the MTZ files fulfilled criteria; will carry forward all results:')
        else:
            self.Logfile.insert('will carry forward the MTZ files from the following auto-processing pipelines:')
        self.Logfile.insert('{0:30} {1:10} {2:10} {3:10}'.format('pipeline','Rmerge(Low)','PG','Score'))
        self.Logfile.insert('----------------------------------------------------------------------')
        for resultDict in dbListOut:
            self.Logfile.insert('{0:30} {1:10} {2:10} {3:10}'.format(resultDict['DataProcessingProgram'],
                                                                     resultDict['DataProcessingRmergeLow'],
                                                                     resultDict['DataProcessingPointGroup'],
                                                                     resultDict['DataProcessingScore']))
        return dbListOut

    def selectResultsSimilarToReference(self,dbList):
        self.Logfile.insert('checking if MTZ files are similar to supplied reference files')
        dbListOut = []
        for resultDict in dbList:
            try:
                if isinstance(float(resultDict['DataProcessingUnitCellVolume']),float):
                    self.Logfile.insert('checking unit cell volume difference and point group:')
                    for reference_file in self.reference_file_list:
                        if not reference_file[4]==0:
                            unitcell_difference=round((math.fabs(reference_file[4]-float(resultDict['DataProcessingUnitCellVolume']))/reference_file[4])*100,1)
                            self.Logfile.insert(resultDict['DataProcessingProgram'] + ': ' +
                                                str(unitcell_difference) + '% difference -> pg(ref): ' +
                                                reference_file[5] + ' -> pg(mtz): ' + resultDict['DataProcessingPointGroup'])
                            if unitcell_difference < self.acceptable_unitcell_volume_difference and reference_file[5]==resultDict['DataProcessingPointGroup']:
                                self.Logfile.insert('=> passed -> mtz file has same point group as reference file and similar unit cell volume')
                                dbListOut.append(resultDict)
                            else:
                                self.Logfile.warning('mtz file has different point group/ unit cell volume as reference file')
            except ValueError:
                pass
        dbListOut = self.report_forward_carried_pipelines(dbListOut,dbList)
        return dbListOut


    def selectResultsWithAcceptableLowResoRmerge(self,dbList):
        self.Logfile.insert('checking if MTZ files have acceptable low resolution Rmerge values (currently set to %s)' %str(self.acceptable_low_resolution_Rmerge))
        dbListOut = []
        for resultDict in dbList:
            try:
                if float(resultDict['DataProcessingRmergeLow']) < self.acceptable_low_resolution_Rmerge:
                    self.Logfile.insert(resultDict['DataProcessingProgram'] +
                                        ': Rmerge(low) of MTZ file is below threshold: ' +
                                        str(resultDict['DataProcessingRmergeLow']))
                    dbListOut.append(resultDict)
                else:
                    self.Logfile.warning(resultDict['DataProcessingProgram'] +
                                         ': Rmerge(low) of MTZ file is ABOVE threshold: ' +
                                         str(resultDict['DataProcessingRmergeLow']))
            except ValueError:
                pass
        dbListOut = self.report_forward_carried_pipelines(dbListOut,dbList)
        return dbListOut


    def selectHighestScore(self,dbList):
        tmp = []
        for resultDict in dbList:
            try:
                tmp.append(float(resultDict['DataProcessingScore']))
            except ValueError:
                tmp.append(0.0)
        highestScoreDict = dbList[tmp.index(max(tmp))]
        return highestScoreDict

    def selectHighestResolution(self,dbList):
        tmp = []
        for resultDict in dbList:
            try:
                tmp.append(float(resultDict['DataProcessingResolutionHigh']))
            except ValueError:
                tmp.append(100.0)
        highestResoDict = dbList[tmp.index(min(tmp))]
        return highestResoDict

    def selectSpecificPipelineOnly(self,dbList):
        tmp = []
        for resultDict in dbList:
            if self.selection_mechanism == 'dials - only' and 'dials' in resultDict['DataProcessingProgram']:
                tmp.append(resultDict)
            if self.selection_mechanism == 'xia2 3d - only' and '3d-run' in resultDict['DataProcessingProgram']:
                tmp.append(resultDict)
            if self.selection_mechanism == 'xia2 3dii - only' and '3dii-run' in resultDict['DataProcessingProgram']:
                tmp.append(resultDict)
            if self.selection_mechanism == 'autoProc - only' and resultDict['DataProcessingProgram'] == 'autoPROC':
                tmp.append(resultDict)
            if self.selection_mechanism == 'autoProc_staraniso - only' and resultDict['DataProcessingProgram'] == 'aP_staraniso':
                tmp.append(resultDict)
        if not tmp:
            tmp = dbList
        pipelineDict = self.selectHighestScore(tmp)
        return pipelineDict


    def determine_processing_outcome(self,db_dict):
        outcome = 'Failed - unknown'
        try:
            if float(db_dict['DataProcessingResolutionHigh']) < self.acceptable_low_resolution_limit_for_data:
                outcome = 'success'
            else:
                outcome = 'Failed - low resolution'
        except ValueError:
            pass
        return outcome


    def updateDB(self,sample, dbDict):
        self.Logfile.insert('{0!s}: updating database'.format(sample))
        self.db.update_insert_data_source(sample, dbDict)


class read_write_autoprocessing_results_from_to_disc(QtCore.QThread):

    """
    major changes:
    - pkl file is obsolete
    - results for every autoprocessing result is recorded in new DB table
    - crystal centring images are copied into project directory
    - beamline directory in project directory will not be used anymore
    - users need to actively select the actual data collection visit as Data Collection Directory in the
      settings tab (e.g. /dls/i04-1/data/2017/mx15433-50
    - results from fast_dp are not copied over and included in analysis

    - DB mainTable gets flag if user updated autoprocessing selection
    - checking of reprocessed files needs to be explicit
    - all dictionaries used to store information are retired
    - at the moment one can only review/ rescore crystals collected during the selected visit
    - parsing of pinIDs in gda logfiles is still missing
    """
    def __init__(self,
                 processedDir,
                 database,
                 projectDir,
                 xce_logfile,
                 target):
        QtCore.QThread.__init__(self)
#        self.target = target
        self.processedDir =  processedDir
        self.visit,self.beamline = XChemMain.getVisitAndBeamline(self.processedDir)
        self.projectDir = projectDir
        self.Logfile = XChemLog.updateLog(xce_logfile)
        self.target = target

        self.db = XChemDB.data_source(os.path.join(database))
        self.exisitingSamples = self.getExistingSamples()

        self.toParse = [
                [   os.path.join('xia2', '*'),
                    os.path.join('LogFiles', '*aimless.log'),
                    os.path.join('DataFiles', '*free.mtz')],
                [   os.path.join('multi-xia2', '*'),
                    os.path.join('LogFiles', '*aimless.log'),
                    os.path.join('DataFiles', '*free.mtz')],
                [   os.path.join('autoPROC', '*'),
                    '*aimless.log',
                    '*truncate-unique.mtz'],
                [   os.path.join('autoPROC', '*'),
                    '*staraniso_alldata-unique.table1',
                    '*staraniso_alldata-unique.mtz'],
                [   os.path.join('autoPROC'),
                    '*aimless.log',
                    '*truncate-unique.mtz'],
                [   os.path.join('autoPROC'),
                    '*staraniso_alldata-unique.table1',
                    '*staraniso_alldata-unique.mtz'],
                [   os.path.join('*'),
                    os.path.join('LogFiles', '*aimless.log'),
                    os.path.join('DataFiles', '*free.mtz')]
                        ]

    def run(self):
        self.parse_file_system()

    def getExistingSamples(self):
        existingSamples={}
        self.Logfile.insert('reading existing samples from collectionTable')
        allEntries = self.db.execute_statement('select CrystalName,DataCollectionVisit,DataCollectionRun,DataProcessingProgram from collectionTable where DataCollectionOutcome = "success"')
        for item in allEntries:
            if str(item[0]) not in existingSamples:
                existingSamples[str(item[0])]=[]
            self.Logfile.insert('%s: adding %s' %(str(item[0]),str(item[1])))
            existingSamples[str(item[0])].append(str(item[1])+ '-' + str(item[2])+str(item[3]))
        return existingSamples

    def createSampleDir(self,xtal):
        if not os.path.isdir(os.path.join(self.projectDir,xtal)):
            os.mkdir(os.path.join(self.projectDir,xtal))

    def createAutoprocessingDir(self,xtal,run,autoproc):
        # create all the directories if necessary
        if not os.path.isdir(os.path.join(self.projectDir,xtal,'autoprocessing')):
            os.mkdir(os.path.join(self.projectDir,xtal,'autoprocessing'))
        if not os.path.isdir(
                os.path.join(self.projectDir,xtal,'autoprocessing', self.visit + '-' + run + autoproc)):
            os.mkdir(os.path.join(self.projectDir,xtal,'autoprocessing', self.visit + '-' + run + autoproc))

    def copyMTZandLOGfiles(self,xtal,run,autoproc,mtzfile,logfile):
        mtzNew = ''
        logNew = ''
        os.chdir(os.path.join(self.projectDir,xtal,'autoprocessing', self.visit + '-' + run + autoproc))
        # MTZ file
        if not os.path.isfile(mtzfile[mtzfile.rfind('/') + 1:]):
            self.Logfile.insert('%s: copying %s' % (xtal, mtzfile))
            os.system('/bin/cp ' + mtzfile + ' .')
        if os.path.isfile(mtzfile[mtzfile.rfind('/')+1:]) and not os.path.isfile(xtal+'.mtz'):
            os.symlink(mtzfile[mtzfile.rfind('/')+1:], xtal + '.mtz')
        if os.path.isfile(mtzfile[mtzfile.rfind('/') + 1:]):
            mtzNew=os.path.join(self.projectDir,xtal,'autoprocessing', self.visit + '-' + run + autoproc, mtzfile[mtzfile.rfind('/')+1:])
        # AIMLESS logfile
        if not os.path.isfile(logfile[logfile.rfind('/')+1:]):
            self.Logfile.insert('%s: copying %s' % (xtal, logfile))
            os.system('/bin/cp ' + logfile + ' .')
        if os.path.isfile(logfile[logfile.rfind('/')+1:]) and not os.path.isfile(xtal+'.log'):
            os.symlink(logfile[logfile.rfind('/')+1:], xtal + '.log')
        if os.path.isfile(logfile[logfile.rfind('/') + 1:]):
            logNew=os.path.join(self.projectDir,xtal,'autoprocessing', self.visit + '-' + run + autoproc, logfile[logfile.rfind('/') + 1:])
        return mtzNew,logNew

    def makeJPGdir(self,xtal,run):
        if not os.path.isdir(os.path.join(self.projectDir,xtal,'jpg')):
            self.Logfile.insert('making jpg directory in '+xtal)
            os.mkdir(os.path.join(self.projectDir,xtal,'jpg'))
        if not os.path.isdir(
                os.path.join(self.projectDir,xtal,'jpg', self.visit + '-' + run)):
            os.mkdir(os.path.join(self.projectDir,xtal,'jpg', self.visit + '-' + run))

    def copyJPGs(self,xtal,run):
        for img in glob.glob(os.path.join(self.processedDir.replace('processed','jpegs'),xtal,run+'*t.png')):
            if not os.path.isfile(os.path.join(self.projectDir,xtal,'jpg', self.visit +'-'+ run,img[img.rfind('/')+1:])):
                os.system('/bin/cp %s %s' %(img,os.path.join(self.projectDir,xtal,'jpg', self.visit + '-' + run)))

    def findJPGs(self,xtal,run):
        jpgDict={}
        for n,img in enumerate(glob.glob(os.path.join(self.projectDir,xtal,'jpg', self.visit +'-'+ run,'*t.png'))):
            if n <= 3: jpgDict['DataCollectionCrystalImage'+str(n+1)]=img
        return jpgDict


    def readProcessingUpdateResults(self,xtal,folder,log,mtz,timestamp,current_run,autoproc):
        db_dict = {}
        for mtzfile in glob.glob(os.path.join(folder,mtz)):
            for logfile in glob.glob(os.path.join(folder, log)):
                self.createAutoprocessingDir(xtal, current_run, autoproc)
                mtzNew,logNew = self.copyMTZandLOGfiles(xtal,current_run,autoproc,mtzfile,logfile)
                if self.target == '=== project directory ===':
                    target = 'unknown'
                else:
                    target = self.target
                db_dict = { 'DataCollectionDate':               timestamp,
                            'DataProcessingPathToLogfile':      logNew,
                            'DataProcessingPathToMTZfile':      mtzNew,
                            'DataProcessingDirectoryOriginal':  folder,
                            'DataCollectionOutcome':            'success',  # success in collection Table only means that a logfile was found
                            'ProteinName':                      target     }
                db_dict.update(parse().read_aimless_logfile(logNew))
                db_dict.update(self.findJPGs(xtal,current_run))     # image exist even if data processing failed
                db_dict['DataCollectionBeamline'] = self.beamline
                print '====>',db_dict['DataProcessingProgram'],db_dict['DataCollectionRun']
                self.update_data_collection_table(xtal,current_run,autoproc,db_dict)
#        return db_dict

    def getAutoProc(self,folder):
        autoproc='unkown'
        if 'ap-run' in folder:
            autoproc = 'autoPROC'
        else:
            autoproc = folder.split('/')[len(folder.split('/'))-1]
        return autoproc

    def update_data_collection_table(self,xtal,current_run,autoproc,db_dict):
        condition_dict = {  'CrystalName':              xtal,
                            'DataCollectionVisit':      self.visit,
                            'DataCollectionRun':        current_run,
                            'DataProcessingProgram':    autoproc    }
        self.db.update_insert_any_table('collectionTable', db_dict, condition_dict)

    def alreadyParsed(self,xtal,current_run,autoproc):
        parsed=False
        if xtal in self.exisitingSamples:
            if self.visit + '-' + current_run + autoproc in self.exisitingSamples[xtal]:
                self.Logfile.insert(
                    '%s: results from %s already parsed; skipping...' % (
                        xtal, self.visit + '-' + current_run + autoproc))
                parsed=True
        return parsed

    def empty_folder(self,xtal,folder):
        empty = True
        stuff = []
        for x in glob.glob(os.path.join(folder,'*')):
            stuff.append(x)
        if not stuff:
            self.Logfile.warning(
                '{0!s}: {1!s} is empty; probably waiting for autoprocessing to finish; try later!'.format(xtal, folder))
        else:
            empty = False
        return empty

    def junk(self,folder):
        do_not_parse = False
        if not os.path.isdir(folder):
            do_not_parse = True
        elif 'dimple' in folder:
            do_not_parse = True
        return do_not_parse

    def parse_file_system(self):
        self.Logfile.insert('checking for new data processing results in '+self.processedDir)
        progress = 0
        progress_step = XChemMain.getProgressSteps(len(glob.glob(os.path.join(self.processedDir,'*'))))

        for collected_xtals in sorted(glob.glob(os.path.join(self.processedDir,'*'))):
            if 'tmp' in collected_xtals or 'results' in collected_xtals or 'scre' in collected_xtals:
                continue
            if not os.path.isdir(collected_xtals):
                self.Logfile.warning(collected_xtals + ' is not a directory')
                continue

            xtal = collected_xtals[collected_xtals.rfind('/')+1:]
            self.Logfile.insert('%s: checking auto-processing results' %xtal)
            self.createSampleDir(xtal)

            if self.target == '=== project directory ===':
                runDir = os.path.join(collected_xtals,'processed','*')
            else:
                runDir = os.path.join(collected_xtals,'*')

            for run in sorted(glob.glob(runDir)):
                current_run=run[run.rfind('/')+1:]
                timestamp=datetime.fromtimestamp(os.path.getmtime(run)).strftime('%Y-%m-%d %H:%M:%S')

                # create directory for crystal aligment images in projectDir
                self.makeJPGdir(xtal,current_run)
                self.copyJPGs(xtal, current_run)

                for item in self.toParse:
                    procDir = os.path.join(run,item[0])
                    logfile = item[1]
                    mtzfile = item[2]

                    for folder in glob.glob(procDir):
                        if self.junk(folder):
                            continue
                        if self.empty_folder(xtal,folder):
                            continue
                        if 'staraniso' in logfile:
                            autoproc = 'aP_staraniso'
                        else:
                            autoproc = self.getAutoProc(folder)
                        if self.alreadyParsed(xtal,current_run,autoproc):
                            continue
                        self.readProcessingUpdateResults(xtal,folder,logfile,mtzfile,timestamp,current_run,autoproc)
#                        db_dict = self.readProcessingResults(xtal,folder,logfile,mtzfile,timestamp,current_run,autoproc)
#                        self.update_data_collection_table(xtal,current_run,autoproc,db_dict)

            progress += progress_step
            self.emit(QtCore.SIGNAL('update_progress_bar'), progress)

        self.Logfile.insert('====== finished parsing beamline directory ======')
        self.emit(QtCore.SIGNAL('read_pinIDs_from_gda_logs'))
        self.emit(QtCore.SIGNAL("finished()"))
