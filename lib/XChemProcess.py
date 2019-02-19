# last edited: 24/03/2017, 15:00

import os, sys, glob
from PyQt4 import QtGui, QtCore

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))
import XChemLog
import XChemDB

class run_xia2(QtCore.QThread):
    def __init__(self,initial_model_directory,run_dict,protocol,spg,ref,reso_limit,cc_half,xce_logfile,external_software,ccp4_scratch_directory,max_queue_jobs,database,overwrite):
        QtCore.QThread.__init__(self)
        self.initial_model_directory=initial_model_directory
        self.run_dict=run_dict
        self.protocol=protocol
        self.spg=spg
        self.ref=ref
        self.reso_limit=reso_limit
        self.cc_half=cc_half
        self.xce_logfile=xce_logfile
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.external_software=external_software
        self.ccp4_scratch_directory=ccp4_scratch_directory
        self.max_queue_jobs=max_queue_jobs
        self.database=database
        self.db=XChemDB.data_source(database)
        self.overwrite=overwrite

    def run(self):
        os.chdir(os.path.join(self.initial_model_directory))
        # first create directories if they do not exist

        if not self.protocol:
            self.Logfile.insert('please select data processing protocol first!')
            return None

        for i,xtal in enumerate(self.run_dict):

            script=''

            if self.external_software['qsub']:
                script+='#PBS -joe -N XCE_reprocess\n'
            else:
                script+='#!'+os.getenv('SHELL')+'\n'

            if 'bash' in os.getenv('SHELL'):
                script += 'export XChemExplorer_DIR="'+os.getenv('XChemExplorer_DIR')+'"\n'
            elif 'csh' in os.getenv('SHELL'):
                script += 'setenv XChemExplorer_DIR '+os.getenv('XChemExplorer_DIR')+'\n'

            if os.getcwd().startswith('/dls'):
                script += 'module load ccp4\n'
#                script += 'module load phenix\n'

# 2018-09-19 - removed this for Joe's test
#                script += 'module load ccp4\n'
#                script += 'module load XDS\n'
#                script += 'module load ccp4 xia2\n'
#                print 'hallo'

            if not self.spg:
                spg_option=''
            else:
                spg_option='space_group='+str(self.spg[0])

            if not self.ref:
                ref_option=''
            else:
                ref_option='reference_reflection_file='+str(self.ref[0])

            if not self.reso_limit:
                reso_limit_option=''
            else:
                reso_limit_option='misigma='+str(self.reso_limit[0])

            if not self.cc_half:
                cc_half_option=''
            else:
                cc_half_option='cc_half='+str(self.cc_half[0])


            # first link diffraction images into directory
            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
                os.mkdir(os.path.join(self.initial_model_directory,xtal))
            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'diffraction_images')):
                os.mkdir(os.path.join(self.initial_model_directory,xtal,'diffraction_images'))

            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'processed')):
                os.mkdir(os.path.join(self.initial_model_directory,xtal,'processed'))
# 2018-09-19 - disabled this for Joe's test
#            if os.path.isfile(os.path.join(self.initial_model_directory,xtal,'processed','run_in_progress')):
#                self.Logfile.insert('data processing is in progress; skipping...')
#                continue
            f = False
            if f:
                print 'hallo'
            else:
                if self.overwrite:
                    os.chdir(os.path.join(self.initial_model_directory,xtal,'diffraction_images'))
                    self.Logfile.insert('removing all links to diffraction images in '+os.path.join(self.initial_model_directory,xtal,'diffraction_images'))
                    os.system('/bin/rm -fr *')
                    os.chdir(os.path.join(self.initial_model_directory,xtal,'processed'))
                    self.Logfile.insert('removing all links to diffraction images in '+os.path.join(self.initial_model_directory,xtal,'processed'))
                    os.system('/bin/rm -fr *')
                os.chdir(os.path.join(self.initial_model_directory,xtal,'processed'))
# 2018-09-19 - disabled this for Joe's test
#                os.system('touch run_in_progress')

            for n,entry in enumerate(sorted(self.run_dict[xtal])):
                newRun=1
                for runDir in sorted(glob.glob(os.path.join(self.initial_model_directory,xtal,'diffraction_images','run_*'))):
                    newRun=int(runDir[runDir.rfind('_')+1:])+1

                os.mkdir(os.path.join(self.initial_model_directory,xtal,'diffraction_images','run_'+str(newRun)))
                image_dir=os.path.join(self.initial_model_directory,xtal,'diffraction_images','run_'+str(newRun))
                os.chdir(os.path.join(self.initial_model_directory,xtal,'diffraction_images','run_'+str(newRun)))
                os.system('ln -s '+os.path.join(entry[0],entry[1])+'* .')

                os.chdir(os.path.join(self.initial_model_directory,xtal,'processed'))
                os.mkdir(os.path.join(self.initial_model_directory,xtal,'processed','run_'+str(newRun)))

                for pipeline in self.protocol:
                    script+='cd '+os.path.join(self.initial_model_directory,xtal,'processed','run_'+str(newRun),pipeline)+'\n'
                    if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'processed','run_'+str(newRun),pipeline)):
                        os.mkdir(os.path.join(self.initial_model_directory,xtal,'processed','run_'+str(newRun),pipeline))

                    script+='$CCP4/bin/ccp4-python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','update_status_flag.py')+' {0!s} {1!s} {2!s} {3!s}\n'.format(self.database, xtal, 'DataProcessingStatus', 'running')

# 2018-09-19 - changed this for Joe's test
#                    script+='$CCP4/bin/xia2 pipeline='+pipeline+' '+ref_option+' '+spg_option+' '+reso_limit_option+' '+cc_half_option+' '+image_dir+'\n'
                    script+='xia2 pipeline='+pipeline+' '+ref_option+' '+spg_option+' '+reso_limit_option+' '+cc_half_option+' '+image_dir+'\n'

            script+='$CCP4/bin/ccp4-python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','update_status_flag.py')+' {0!s} {1!s} {2!s} {3!s}\n'.format(self.database, xtal, 'DataProcessingStatus', 'finished')
            script+='cd '+os.path.join(self.initial_model_directory,xtal,'processed')+'\n'
            script+='/bin/rm run_in_progress\n'

            os.chdir(self.ccp4_scratch_directory)
            f = open('xce_xia2_{0!s}.sh'.format(str(i+1)),'w')
            f.write(script)
            f.close()
            os.system('chmod +x xce_xia2_{0!s}.sh'.format(str(i+1)))
            db_dict= {'DataProcessingStatus': 'started'}
            self.Logfile.insert('{0!s}: setting DataProcessingStatus flag to started'.format(xtal))
            self.db.update_data_source(xtal,db_dict)

        # submit job
        self.Logfile.insert('created input scripts for '+str(n+1)+' in '+self.ccp4_scratch_directory)
        os.chdir(self.ccp4_scratch_directory)
        if os.getcwd().startswith('/dls'):
            if self.external_software['qsub_array']:
                Cmds = (
                        '#PBS -joe -N xce_xia2_master\n'
                        './xce_xia2_$SGE_TASK_ID.sh\n'
                        )
                f = open('xia2_master.sh','w')
                f.write(Cmds)
                f.close()
                self.Logfile.insert('submitting array job with maximal 100 jobs running on cluster')
                self.Logfile.insert('using the following command:')
                self.Logfile.insert('qsub -t 1:{0!s} -tc {1!s} xia2_master.sh'.format(str(i+1), self.max_queue_jobs))
                os.system('qsub -P labxchem -q medium.q -N xia2 -t 1:{0!s} -tc {1!s} xia2_master.sh'.format(str(i+1), self.max_queue_jobs))
            else:
                self.Logfile.insert("cannot start ARRAY job: make sure that 'module load global/cluster' is in your .bashrc or .cshrc file")
        elif self.external_software['qsub']:
            self.Logfile.insert('submitting {0!s} individual jobs to cluster'.format((str(i+1))))
            self.Logfile.insert('WARNING: this could potentially lead to a crash...')
            for n in range(i+1):
                self.Logfile.insert('qsub xce_xia2_{0!s}.sh'.format((str(n+1))))
                os.system('qsub -q medium.q -N xia2 xce_xia2_{0!s}.sh'.format((str(n+1))))
        else:
            self.Logfile.insert('running %s consecutive XIA2 jobs on your local machine')
            for n in range(i+1):
                self.Logfile.insert('starting xce_xia2_{0!s}.sh'.format((str(n+1))))
                os.system('./xce_xia2_{0!s}.sh'.format((str(n+1))))


#            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal)):
#                os.mkdir(os.path.join(self.initial_model_directory,xtal))
#            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'autoprocessing')):
#                os.mkdir(os.path.join(self.initial_model_directory,xtal,'autoprocessing'))
#
#
#
#
#            if not os.path.isdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc)):
#                os.mkdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
#            os.chdir(os.path.join(self.initial_model_directory,xtal,'dimple',visit_run_autoproc))
#            os.system('touch dimple_run_in_progress')




#        header='#!'+os.getenv('SHELL')+'\n'
#        if external_software['qsub']:
#            if not external_software['qsub_array']:
#                header='#PBS -joe -N xce_acedrg\n'
#
#        Cmds = (
#                    header+
#                    '\n'
#                    'export XChemExplorer_DIR="'+os.getenv('XChemExplorer_DIR')+'"\n'
#                    '\n'
#                    'source $XChemExplorer_DIR/setup-scripts/xce.setup-sh\n'
#                    '\n'
#                    '$CCP4/bin/ccp4-python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','create_png_of_compound.py')+
#                    ' "%s" %s %s %s\n' %(smiles,compoundID.replace(' ',''),sample,initial_model_directory)+
#                    '\n'
#                    'cd '+os.path.join(initial_model_directory,sample,'compound')+'\n'
#                    '\n'
#                    'acedrg --res LIG -i "%s" -o %s\n' %(smiles,compoundID.replace(' ',''))+
#                    '\n'
#                    'cd '+os.path.join(initial_model_directory,sample)+'\n'
#                    '\n'
#                    'ln -s compound/%s.cif .\n' %compoundID.replace(' ','')+
#                    'ln -s compound/%s.pdb .\n' %compoundID.replace(' ','')+
#                    'ln -s compound/%s.png .\n' %compoundID.replace(' ','')+
#                    '\n'
#                    '$CCP4/bin/ccp4-python '+os.path.join(os.getenv('XChemExplorer_DIR'),'helpers','update_data_source_for_new_cif_files.py')+
#                    ' %s %s %s %s\n' %(os.path.join(database_directory,data_source_file),sample,initial_model_directory,compoundID.replace(' ','') )+
#                    '\n'
#                    '/bin/rm compound/ACEDRG_IN_PROGRESS\n'
