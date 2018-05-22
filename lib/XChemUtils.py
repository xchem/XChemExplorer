# last edited: 07/07/2017, 17:00

import sys
import os
import glob
import math
import subprocess
import getpass
import shutil
import math
import platform


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

import iotbx.pdb

from iotbx.reflection_file_reader import any_reflection_file
from iotbx import mtz

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
import XChemDB
import XChemLog


class process:

    def __init__(self,dimple):
        self.project_directory=dimple['project_directory']
        self.delete_old=dimple['delete_old']
        self.xtalID=dimple['xtalID']
        self.reference=dimple['reference']
        self.queueing_system_available=dimple['queueing_system_available']
        self.ccp4_scratch_directory=dimple['ccp4_scratch']
        self.fileroot_in=dimple['fileroot_in']
        self.mtz_free=self.fileroot_in+'.free.mtz'
        self.mtz_in=self.fileroot_in+'.mtz'
        # need to set $CCP4_SCR, because the default directory in /home/zqr16691/tmp fills up easily
        # and then dimple will not run
#        if not os.path.isdir(self.project_directory[:self.project_directory.find('processing')+11]+'/tmp'):
#            os.mkdir(self.project_directory[:self.project_directory.find('processing')+11]+'/tmp')
#        self.ccp4_scratch=self.project_directory[:self.project_directory.find('processing')+11]+'/tmp'


    def get_Rfree(self):
        Cmds=''
        if not os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.mtz_free)):
            if os.path.isfile(os.path.join(self.project_directory,self.xtalID,self.reference+'.mtz')):
                Cmds = (
#                    '#!'+os.getenv('SHELL')+'\n'
                    '\n'
                    'cd %s/%s/Dimple\n' %(self.project_directory,self.xtalID) +
                    '\n'
                    'pointless hklin ../%s hklref %s.mtz hklout %s.reind.mtz << EOF > pointless.log\n' %(self.mtz_in,self.reference,self.xtalID)+
                    ' tolerance 5\n'
                    'EOF\n'
                    '\n'
                    'cad hklin1 %s.reind.mtz hklin2 %s.mtz hklout cad.mtz << EOF > cad.log\n' %(self.xtalID,self.reference) +
                    '   labin file_number 1 E1=F E2=SIGF E3=IMEAN E4=SIGIMEAN\n'
                    '   labout file_number 1 E1=F E2=SIGF E3=IMEAN E4=SIGIMEAN\n'
                    '   labin file_number 2 E1=FreeR_flag\n'
                    '   labout file_number 2 E1=FreeR_flag\n'
                    '   END\n'
                    'EOF\n'
                    '\n'
                    'freerflag hklin cad.mtz hklout ../%s << EOF > freerflag.log\n' %self.mtz_free +
                    '   COMPLETE FREE=FreeR_flag\n'
                    '   END\n'
                    'EOF\n'
                    '#uniqueify -f FreeR_flag cad.mtz ../%s\n' %self.mtz_free
                    )
            else:
                Cmds = (
#                    '#!'+os.getenv('SHELL')+'\n'
                    '\n'
                    'cd %s/%s/Dimple\n' %(self.project_directory,self.xtalID) +
                    '\n'
                    'pointless hklin ../%s xyzin %s.pdb hklout %s.reind.mtz << EOF > pointless.log\n' %(self.mtz_in,self.reference,self.xtalID)+
                    ' tolerance 5\n'
                    'EOF\n'
                    '\n'
                    'cad hklin1 %s.reind.mtz hklout cad.mtz << EOF > cad.log\n' %self.xtalID +
                    '   labin file_number 1 E1=F E2=SIGF E3=IMEAN E4=SIGIMEAN\n'
                    '   labout file_number 1 E1=F E2=SIGF E3=IMEAN E4=SIGIMEAN\n'
                    '   END\n'
                    'EOF\n'
                    '\n'
                    'freerflag hklin cad.mtz hklout ../%s > freerflag.log\n' %self.mtz_free +
                    '#uniqueify cad.mtz ../{0!s}\n'.format(self.mtz_free)
                    )
#            os.chdir(os.path.join(self.project_directory,self.xtalID))
#            os.system(Cmds)
##            os.symlink(os.path.join('Dimple',self.xtalID+'.free.mtz'),self.xtalID+'.free.mtz')
            return Cmds

    def dimple(self):
        os.chdir('{0!s}/{1!s}'.format(self.project_directory, self.xtalID))
        generate_Rfree_Cmds=''
        if os.path.isdir('{0!s}/{1!s}/Dimple'.format(self.project_directory, self.xtalID)):
            # this 'if not' step is usually not necessary; only if process was stopped in an odd way
            if not os.path.isfile(self.project_directory+'/'+self.xtalID+'/'+self.mtz_free):
                generate_Rfree_Cmds=self.get_Rfree()
            if self.delete_old==True:
                os.system('rm -fr Dimple')
                os.mkdir('{0!s}/{1!s}/Dimple'.format(self.project_directory, self.xtalID))
            else:
                return None
        else:
            os.mkdir('{0!s}/{1!s}/Dimple'.format(self.project_directory, self.xtalID))
            generate_Rfree_Cmds=self.get_Rfree()
            os.chdir('{0!s}/{1!s}/Dimple'.format(self.project_directory, self.xtalID))

        if self.queueing_system_available:
            top_line='#PBS -joe -N XCE_dimple'
        else:
            top_line='#!'+os.getenv('SHELL')

        if 'csh' in os.getenv('SHELL'):
            ccp4_scratch='setenv CCP4_SCR '+self.ccp4_scratch_directory+'\n'
        elif 'bash' in os.getenv('SHELL'):
            ccp4_scratch='export CCP4_SCR='+self.ccp4_scratch_directory+'\n'
        else:
            ccp4_scratch=''

        ################################################################################
        # in case additional ligands are present in the reference.pdb file,
        # the user needs to provide the following files:
        # - <my_file_1>.pdb
        # - <my_file_1>.cif     # this file contains restraints for all the ligands in the respective pdb file

        if os.path.isfile(self.reference+'.cif'):
            ref_lib=' --libin '+self.reference+'.cif'
        else:
            ref_lib=''

        Cmds = (
                '{0!s}\n'.format(top_line)+
                generate_Rfree_Cmds+
                '\n'
                'cd %s/%s/Dimple\n' %(self.project_directory,self.xtalID) +
                '\n'
                +ccp4_scratch+
                '\n'
                'dimple %s ../%s %s.pdb dimple' %(ref_lib,self.mtz_free,self.reference) +
                '\n'
                'cd %s/%s\n' %(self.project_directory,self.xtalID) +
                '\n'
                '/bin/rm refine.pdb\n'
                '/bin/rm refine.mtz\n'
                'ln -s Dimple/dimple/final.pdb refine.pdb\n'
                'ln -s Dimple/dimple/final.mtz refine.mtz\n'
                '\n'
                'fft hklin refine.mtz mapout 2fofc.map << EOF\n'
                ' labin F1=2FOFCWT PHI=PH2FOFCWT\n'
                'EOF\n'
                '\n'
                'fft hklin refine.mtz mapout fofc.map << EOF\n'
                ' labin F1=FOFCWT PHI=PHFOFCWT\n'
                'EOF\n'
                '\n'
                'fft hklin refine.mtz mapout 2fofc.map << EOF\n'
                ' labin F1=FWT PHI=PHWT\n'
                'EOF\n'
                '\n'
                'fft hklin refine.mtz mapout fofc.map << EOF\n'
                ' labin F1=DELFWT PHI=PHDELWT\n'
                'EOF\n'
                '\n'
                '/bin/rm dimple_run_in_progress\n'
                )

        self.Run(Cmds)



    def Run(self,Cmds):
        os.chdir(self.project_directory+'/'+self.xtalID)
        os.system('touch dimple_run_in_progress')
        os.chdir(self.project_directory+'/'+self.xtalID+'/Dimple')
        f = open('xce_dimple.sh','w')
        f.write(Cmds)
        f.close()
        if self.queueing_system_available:
            os.system('qsub xce_dimple.sh')
        else:
            os.system('chmod +x xce_dimple.sh')
            os.system('./xce_dimple.sh')

class helpers:


    def __init__(self):
        try:
            # this sys path append is only meaningful if non-CCP4 python is used to launch XCE
#            sys.path.append(os.path.join(os.getenv('CCP4'),'lib','python2.7','site-packages'))
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from rdkit.Chem import Draw
            self.pil_rdkit_present=True
        except ImportError:
            self.pil_rdkit_present=False
            pass

    def pil_rdkit_exist(self):
        return self.pil_rdkit_present

    def make_png(self,initial_model_directory,sample,compoundID,smiles,external_software,database_directory,data_source_file,ccp4_scratch_directory,counter,xce_logfile,restraints_program):
        Logfile=XChemLog.updateLog(xce_logfile)

#        if not os.path.isfile(os.path.join(initial_model_directory,sample,'compound','ACEDRG_IN_PROGRESS')):
        os.system('touch RESTRAINTS_IN_PROGRESS')

        header='#!'+os.getenv('SHELL')+'\n'
        if external_software['qsub']:
            if not external_software['qsub_array']:
                header='#PBS -joe -N xce_acedrg\n'

        software=''
        if restraints_program=='acedrg':
            if os.path.isfile(os.path.join(initial_model_directory,sample,'old.cif')):
                software='acedrg --res LIG -c ../old.cif -o {0!s}\n'.format((compoundID.replace(' ','')))
            else:
                software='acedrg --res LIG -i "{0!s}" -o {1!s}\n'.format(smiles, compoundID.replace(' ',''))
        elif restraints_program=='phenix.elbow':
            if os.path.isfile(os.path.join(initial_model_directory,sample,'old.cif')):
                software='phenix.elbow --file=../old.cif --id LIG --output {0!s}\n'.format((compoundID.replace(' ','')))
            else:
                software='phenix.elbow --smiles="{0!s}" --id LIG --output {1!s}\n'\
                    .format(smiles, compoundID.replace(' ',''))
        elif restraints_program=='grade':
            if os.getcwd().startswith('/dls'):
                software+='module load buster\n'
            software+="export BDG_TOOL_OBABEL='none'\n"
            if os.path.isfile(os.path.join(initial_model_directory,sample,'old.cif')):
                software+='grade -resname LIG -nomogul -in ../old.cif -ocif {0!s}.cif -opdb {1!s}.pdb\n'\
                    .format(compoundID.replace(' ',''), compoundID.replace(' ',''))
            else:
                software+='grade -resname LIG -nomogul "{0!s}" -ocif {1!s}.cif -opdb {2!s}.pdb\n'\
                    .format(smiles, compoundID.replace(' ',''), compoundID.replace(' ',''))
        # Removal of the hydrogen atoms in PDB files is required for REFMAC 5 run. With hydrogens some ligands fail to
        # pass the external restraints in pandda.giant.make_restraints.
        # Copy the file with hydrogens to retain in case needed

        copy_with_hydrogens = 'cp {0!s}.pdb {0!s}_with_H.pdb'.format(compoundID.replace(' ', ''))

        strip_hydrogens = 'phenix.reduce {0!s}.pdb -trim > {0!s}_tmp.pdb'.format(compoundID.replace(' ', ''))

        Cmds = (

            header +
            '\n'
            'export XChemExplorer_DIR="' + os.getenv('XChemExplorer_DIR') + '"' +
            '\n'
            'source $XChemExplorer_DIR/setup-scripts/xce.setup-sh'
            '\n'
            '$CCP4/bin/ccp4-python $XChemExplorer_DIR/helpers/update_status_flag.py {0!s} {1!s} {2!s} {3!s}'
            .format(os.path.join(database_directory, data_source_file), sample, 'RefinementCIFStatus', 'running') +
            '\n'
            '$CCP4/bin/ccp4-python {0!s} "{1!s}" {2!s} {3!s} {4!s}\n'
            .format(os.path.join(os.getenv('XChemExplorer_DIR'), 'helpers','create_png_of_compound.py'),
                    smiles, compoundID.replace(' ', ''), sample, initial_model_directory) +
            '\n'
            'cd ' + os.path.join(initial_model_directory, sample, 'compound') +
            '\n'
            + software +
            '\n'
            + copy_with_hydrogens +
            '\n'
            + strip_hydrogens +
            '\n'
            'mv {0!s}_tmp.pdb {0!s}.pdb'.format(compoundID.replace(' ', '')) +
            '\n'
            'rm -f {0!s}_tmp.pdb'.format(compoundID.replace(' ', '')) +
            '\n'
            'cd ' + os.path.join(initial_model_directory, sample) +
            '\n'
            'ln -s compound/%s.cif .\n' % compoundID.replace(' ', '') +
            'ln -s compound/{0!s}.pdb .\n'.format(compoundID.replace(' ', '')) +
            'ln -s compound/{0!s}.png .\n'.format(compoundID.replace(' ', '')) +
            '\n'
            '$CCP4/bin/ccp4-python ' + os.path.join(os.getenv('XChemExplorer_DIR'), 'helpers',
                                                    'update_data_source_for_new_cif_files.py') +
            ' {0!s} {1!s} {2!s} {3!s}\n'.format(os.path.join(database_directory, data_source_file), sample,
                                                initial_model_directory, compoundID.replace(' ', '')) +
            '\n'
            '/bin/rm -f compound/RESTRAINTS_IN_PROGRESS\n'
        )


        os.chdir(ccp4_scratch_directory)
        Logfile.insert('creating ACEDRG shell script for {0!s},{1!s} in {2!s}'.format(sample, compoundID, ccp4_scratch_directory))
        print Cmds
        print 'ccp4_scratch',ccp4_scratch_directory
        f = open('xce_acedrg_{0!s}.sh'.format(str(counter)),'w')
        f.write(Cmds)
        f.close()
        os.system('chmod +x xce_acedrg_{0!s}.sh'.format(str(counter)))



#                if queueing_system_available:
#                    os.system('qsub acedrg.sh')
#                else:
#                    os.system('chmod +x acedrg.sh')
#                    os.system('./acedrg.sh')

#        os.chdir(os.path.join(initial_model_directory,sample))
#        if os.path.isfile(os.path.join(initial_model_directory,sample,'compound',compoundID.replace(' ','')+'.pdb'))\
#                and not os.path.isfile(os.path.join(initial_model_directory,sample,compoundID.replace(' ','')+'.pdb')):
#            os.symlink(os.path.join(initial_model_directory,sample,'compound',compoundID.replace(' ','')+'.pdb'),compoundID.replace(' ','')+'.pdb')
#        if os.path.isfile(os.path.join(initial_model_directory,sample,'compound',compoundID.replace(' ','')+'.cif'))\
#                and not os.path.isfile(os.path.join(initial_model_directory,sample,compoundID.replace(' ','')+'.cif')):
#            os.symlink(os.path.join(initial_model_directory,sample,'compound',compoundID.replace(' ','')+'.cif'),compoundID.replace(' ','')+'.cif')
#        if os.path.isfile(os.path.join(initial_model_directory,sample,'compound',compoundID.replace(' ','')+'.png'))\
#                and not os.path.isfile(os.path.join(initial_model_directory,sample,compoundID.replace(' ','')+'.png')):
#            os.symlink(os.path.join(initial_model_directory,sample,'compound',compoundID.replace(' ','')+'.png'),compoundID.replace(' ','')+'.png')


class parse:

    def __init__(self):
        self.space_group_dict=   {  'triclinic':    ['P1'],
#                                    'monoclinic':   ['P2','P21','C121','P1211','P121'],
                                    'monoclinic_P':   ['P2','P21','P1211','P121'],
                                    'monoclinic_C':   ['C2','C121'],
                                    'orthorhombic': ['P222','P2122','P2212','P2221',
                                                     'P21212','P21221','P22121','P212121',
                                                     'C222','C2221',
                                                     'F222', 'I222', 'I212121'],
                                    'tetragonal':   ['P4','P41','P42','P43','I4','I41',
                                                     'P422','P4212','P4122','P41212','P4222',
                                                     'P42212','P4322','P43212','I422','I4122' ],
                                    'hexagonal':    ['P3','P31','P32',
                                                     'P312','P321','P3112','P3121','P3212','P3221',
                                                     'P6','P61','P65','P62','P64','P63',
                                                     'P622','P6122','P6522','P6222','P6422','P6322'    ],
                                    'rhombohedral': ['H3','H32'],
                                    'cubic':        ['P23','F23','I23','P213','I213',
                                                     'P432','P4232','F432','F4132','I432','P4332','P4132','I4132' ] }

        self.point_group_dict=   {  '1':    ['P1'],
                                    '2':    ['P2','P21','C121','P1211','P121'],
                                    '222':  ['P222','P2122','P2212','P2221',
                                             'P21212','P21221','P22121','P212121',
                                             'C222','C2221',
                                             'F222', 'I222', 'I212121'],
                                    '4':    ['P4','P41','P42','P43','I4','I41'],
                                    '422':  ['P422','P4212','P4122','P41212','P4222',
                                             'P42212','P4322','P43212','I422','I4122' ],
                                    '3':    ['P3','P31','P32','H3'],
                                    '32':   ['P312','P321','P3112','P3121','P3212','P3221','H32'],
                                    '6':    ['P6','P61','P65','P62','P64','P63'],
                                    '622':  ['P622','P6122','P6522','P6222','P6422','P6322'    ],
                                    '23':   ['P23','F23','I23','P213','I213'],
                                    '432':  ['P432','P4232','F432','F4132','I432','P4332','P4132','I4132' ] }

        self.nr_asu_in_unitcell_for_point_group =   {   '1':            1,
                                                        '2':            2,
                                                        '222':          4,
                                                        '4':            4,
                                                        '422':          8,
                                                        '3':            3,
                                                        '32':           6,
                                                        '6':            6,
                                                        '622':          12,
                                                        '23':           12,
                                                        '432':          24  }

        self.aimless = {    'DataProcessingProgram':                        'n/a',
                            'DataCollectionRun':                            'n/a',
                            'DataProcessingSpaceGroup':                     'n/a',
                            'DataProcessingUnitCell':                       'n/a',
                            'DataProcessingA':                              'n/a',
                            'DataProcessingB':                              'n/a',
                            'DataProcessingC':                              'n/a',
                            'DataProcessingAlpha':                          'n/a',
                            'DataProcessingBeta':                           'n/a',
                            'DataProcessingGamma':                          'n/a',
                            'DataProcessingResolutionLow':                  'n/a',
                            'DataProcessingResolutionLowInnerShell':        'n/a',
                            'DataProcessingResolutionHigh':                 'n/a',
                            'DataProcessingResolutionHighOuterShell':       'n/a',
                            'DataProcessingResolutionOverall':              'n/a',
                            'DataProcessingRmergeOverall':                  'n/a',
                            'DataProcessingRmergeLow':                      'n/a',
                            'DataProcessingRmergeHigh':                     'n/a',
                            'DataProcessingIsigOverall':                    'n/a',
                            'DataProcessingIsigLow':                        'n/a',
                            'DataProcessingIsigHigh':                       'n/a',
                            'DataProcessingCompletenessOverall':            'n/a',
                            'DataProcessingCompletenessLow':                'n/a',
                            'DataProcessingCompletenessHigh':               'n/a',
                            'DataProcessingMultiplicityOverall':            'n/a',
                            'DataProcessingMultiplicityLow':                'n/a',
                            'DataProcessingMultiplicityHigh':               'n/a',
                            'DataProcessingCChalfOverall':                  'n/a',
                            'DataProcessingCChalfLow':                      'n/a',
                            'DataProcessingCChalfHigh':                     'n/a',
                            'DataProcessingResolutionHigh15sigma':         'n/a',
                            'DataProcessingUniqueReflectionsOverall':       'n/a',
                            'DataProcessingLattice':                        'n/a',
                            'DataProcessingPointGroup':                     'n/a',
                            'DataProcessingUnitCellVolume':                 0,
                            'DataProcessingAlert':                          '#FF0000',
                            'DataCollectionWavelength':                     'n/a',
                            'DataProcessingScore':                          'n/a'             }



    def GetAimlessLog(self,Logfile):
        self.Logfile=Logfile
        Aimless = { 'AutoProc': 'n/a',
                    'Run': 'n/a',
                    'SpaceGroup': 'n/a',
                    'UnitCell': 'n/a',
                    'ResolutionLow': 'n/a',
                    'ResolutionLowInnerShell': 'n/a',
                    'ResolutionHigh': 'n/a',
                    'ResolutionHighOuterShell': 'n/a',
                    'RmergeOverall': 'n/a',
                    'RmergeLow': 'n/a',
                    'RmergeHigh': 'n/a',
                    'IsigOverall': 'n/a',
                    'IsigLow': 'n/a',
                    'IsigHigh': 'n/a',
                    'CompletenessOverall': 'n/a',
                    'CompletenessLow': 'n/a',
                    'CompletenessHigh': 'n/a',
                    'MultiplicityOverall': 'n/a',
                    'MultiplicityLow': 'n/a',
                    'MultiplicityHigh': 'n/a',
                    'UniqueReflectionsOverall': 'n/a',
                    'Lattice': 'n/a',
                    'UnitCellVolume': 0,
                    'Alert': '#FF0000' }

        spg='n/a'
        a='n/a'
        b='n/a'
        c='n/a'
        alpha='n/a'
        beta='n/a'
        gamma='n/a'

        if 'fast_dp' in self.Logfile:
            Aimless['AutoProc']='fast_dp'
        if '3d-run' in self.Logfile:
            Aimless['AutoProc']='xia2 3d'
        if '3dii-run' in self.Logfile:
            Aimless['AutoProc']='xia2 3dii'
        if 'dials-run' in self.Logfile:
            Aimless['AutoProc']='dials'
        if 'autoPROC' in self.Logfile:
            Aimless['AutoProc']='autoPROC'

        # get run number from logfile
        # Note: only works if file is in original directory, but not once it lifes in 'inital_model'
#        print self.Logfile.split('/')[9].split('_')[1]
#        if len(self.Logfile.split('/'))>8 and len(self.Logfile.split('/')[9].split('_'))==1:
        try:
            Aimless['Run']=self.Logfile.split('/')[9].split('_')[1]
        except IndexError:
            pass

        for line in open(self.Logfile):
            if line.startswith('Low resolution limit') and len(line.split())==6:
                Aimless['ResolutionLow'] = line.split()[3]
                Aimless['ResolutionHighOuterShell'] = line.split()[5]
            if line.startswith('High resolution limit') and len(line.split())==6:
                Aimless['ResolutionHigh'] = line.split()[3]
                Aimless['ResolutionLowInnerShell'] = line.split()[4]
            if line.startswith('Rmerge  (all I+ and I-)') and len(line.split())==8:
                Aimless['RmergeOverall'] = line.split()[5]
                Aimless['RmergeLow'] = line.split()[6]
                Aimless['RmergeHigh']  = line.split()[7]
            if line.startswith('Mean((I)/sd(I))') and len(line.split())==4:
                Aimless['IsigOverall'] = line.split()[1]
                Aimless['IsigHigh'] = line.split()[3]
                Aimless['IsigLow'] = line.split()[2]
            if line.startswith('Completeness') and len(line.split())==4:
                Aimless['CompletenessOverall'] = line.split()[1]
                Aimless['CompletenessHigh'] = line.split()[3]
                Aimless['CompletenessLow'] = line.split()[2]
            if line.startswith('Multiplicity') and len(line.split())==4:
                Aimless['MultiplicityOverall'] = line.split()[1]
                Aimless['MultiplicityHigh'] = line.split()[3]
                Aimless['MultiplicityLow'] = line.split()[3]
            if line.startswith('Average unit cell:') and len(line.split())==9:
                tmp = []
                tmp.append(line.split())
                a = int(float(tmp[0][3]))
                b = int(float(tmp[0][4]))
                c = int(float(tmp[0][5]))
                alpha = int(float(tmp[0][6]))
                beta = int(float(tmp[0][7]))
                gamma = int(float(tmp[0][8]))
            if line.startswith('Total number unique') and len(line.split())==6:
                Aimless['UniqueReflectionsOverall']=line.split()[3]
            if line.startswith('Space group:'):
                Aimless['SpaceGroup']=line.replace('Space group: ','')[:-1]
                Aimless['Lattice']=self.get_lattice_from_space_group(Aimless['SpaceGroup'])
                if a != 'n/a' and b != 'n/a' and c != 'n/a' and \
                   alpha != 'n/a' and beta != 'n/a' and gamma != 'n/a' and Aimless['Lattice'] != 'n/a':
                    Aimless['UnitCellVolume']=self.calc_unitcell_volume_from_logfile(float(a),float(b),float(c),
                                                                                 math.radians(float(alpha)),
                                                                                 math.radians(float(beta)),
                                                                                 math.radians(float(gamma)),
                                                                                 Aimless['Lattice'])
        Aimless['UnitCell']=str(a)+' '+str(b)+' '+str(c)+' '+str(alpha)+' '+str(beta)+' '+str(gamma)

        # Hex Color code:
        # red:      #FF0000
        # orange:   #FF9900
        # green:    #00FF00
        # gray:     #E0E0E0

        if Aimless['ResolutionHigh']=='n/a' or Aimless['RmergeLow'] =='n/a':
            Aimless['Alert'] = '#FF0000'
        else:
            if float(Aimless['ResolutionHigh']) > 3.5 or float(Aimless['RmergeLow']) > 0.1:
                Aimless['Alert'] = '#FF0000'
            if (float(Aimless['ResolutionHigh']) <= 3.5 and float(Aimless['ResolutionHigh']) > 2.5) or \
               (float(Aimless['RmergeLow']) <= 0.1 and float(Aimless['RmergeLow']) > 0.05):
                Aimless['Alert'] = '#FF9900'
            if float(Aimless['ResolutionHigh']) <= 2.5 and float(Aimless['RmergeLow']) <= 0.05:
                Aimless['Alert'] = '#00FF00'

        return Aimless

    def return_empty_aimless_dict_for_db(self):
        return self.aimless

    def read_aimless_logfile(self,logfile):
        # essentially same as above, but compatible with datasource
        # will hopefully supersede function above

        spg='n/a'
        a='n/a'
        b='n/a'
        c='n/a'
        alpha='n/a'
        beta='n/a'
        gamma='n/a'

        if 'fast_dp' in logfile:
            self.aimless['DataProcessingProgram']='fast_dp'
        elif '3d-run' in logfile:
            self.aimless['DataProcessingProgram']='xia2 3d'
        elif '3dii-run' in logfile:
            self.aimless['DataProcessingProgram']='xia2 3dii'
        elif 'dials-run' in logfile:
            self.aimless['DataProcessingProgram']='dials'
        elif 'autoPROC' in logfile:
            self.aimless['DataProcessingProgram']='autoPROC'

        # get run number from logfile
        # Note: only works if file is in original directory, but not once it moved to 'inital_model' folder
#        print self.Logfile.split('/')[9].split('_')[1]
#        if len(self.Logfile.split('/'))>8 and len(self.Logfile.split('/')[9].split('_'))==1:
        try:
            self.aimless['DataCollectionRun']=logfile.split('/')[9].split('_')[1]
        except IndexError:
            pass

        resolution_at_sigma_line_overall_found=False
        for line_number,line in enumerate(open(logfile)):
#            if 'Wavelength' in line:
#                print 'here'
#                print line.split()
            if 'Wavelength' in line and len(line.split())==3:
                self.aimless['DataCollectionWavelength']=line.split()[1]
            if line.startswith('Low resolution limit') and len(line.split())==6:
                self.aimless['DataProcessingResolutionLow'] = line.split()[3]
                self.aimless['DataProcessingResolutionHighOuterShell'] = line.split()[5]
            if line.startswith('High resolution limit') and len(line.split())==6:
                self.aimless['DataProcessingResolutionHigh'] = line.split()[3]
                self.aimless['DataProcessingResolutionLowInnerShell'] = line.split()[4]
            if line.startswith('Rmerge  (all I+ and I-)') and len(line.split())==8:
                self.aimless['DataProcessingRmergeOverall'] = line.split()[5]
                self.aimless['DataProcessingRmergeLow'] = line.split()[6]
                self.aimless['DataProcessingRmergeHigh']  = line.split()[7]
            if line.startswith('Mean((I)/sd(I))') and len(line.split())==4:
                self.aimless['DataProcessingIsigOverall'] = line.split()[1]
                self.aimless['DataProcessingIsigHigh'] = line.split()[3]
                self.aimless['DataProcessingIsigLow'] = line.split()[2]
            if line.startswith('Completeness') and len(line.split())==4:
                self.aimless['DataProcessingCompletenessOverall'] = line.split()[1]
                self.aimless['DataProcessingCompletenessHigh'] = line.split()[3]
                self.aimless['DataProcessingCompletenessLow'] = line.split()[2]
            if line.startswith('Multiplicity') and len(line.split())==4:
                self.aimless['DataProcessingMultiplicityOverall'] = line.split()[1]
                self.aimless['DataProcessingMultiplicityHigh'] = line.split()[3]
                self.aimless['DataProcessingMultiplicityLow'] = line.split()[3]
            if line.startswith('Mn(I) half-set correlation CC(1/2)') and len(line.split())==7:
                self.aimless['DataProcessingCChalfOverall'] = line.split()[4]
                self.aimless['DataProcessingCChalfLow'] = line.split()[5]
                self.aimless['DataProcessingCChalfHigh'] = line.split()[6]
            if line.startswith('Estimates of resolution limits: overall'):
                resolution_at_sigma_line_overall_found=True
            if resolution_at_sigma_line_overall_found:
                if 'from Mn(I/sd)' in line and len(line.split()) >= 7:
                    if '1.5' in line.split()[3]:
                        self.aimless['DataProcessingResolutionHigh15sigma']=line.split()[6][:-1]
                        resolution_at_sigma_line_overall_found=False
            if line.startswith('Average unit cell:') and len(line.split())==9:
                tmp = []
                tmp.append(line.split())
                a = int(float(tmp[0][3]))
                b = int(float(tmp[0][4]))
                c = int(float(tmp[0][5]))
                alpha = int(float(tmp[0][6]))
                beta = int(float(tmp[0][7]))
                gamma = int(float(tmp[0][8]))
                self.aimless['DataProcessingA']=str(a)
                self.aimless['DataProcessingB']=str(b)
                self.aimless['DataProcessingC']=str(c)
                self.aimless['DataProcessingAlpha']=str(alpha)
                self.aimless['DataProcessingBeta']=str(beta)
                self.aimless['DataProcessingGamma']=str(gamma)
            if line.startswith('Total number unique') and len(line.split())==6:
                self.aimless['DataProcessingUniqueReflectionsOverall']=line.split()[3]
            if line.startswith('Space group:'):
                self.aimless['DataProcessingSpaceGroup']=line.replace('Space group: ','')[:-1]
                self.aimless['DataProcessingLattice']=self.get_lattice_from_space_group(self.aimless['DataProcessingSpaceGroup'])
                self.aimless['DataProcessingPointGroup']=self.get_pointgroup_from_space_group(self.aimless['DataProcessingSpaceGroup'])
                if a != 'n/a' and b != 'n/a' and c != 'n/a' and \
                   alpha != 'n/a' and beta != 'n/a' and gamma != 'n/a' and self.aimless['DataProcessingLattice'] != 'n/a':
                    self.aimless['DataProcessingUnitCellVolume']=str(self.calc_unitcell_volume_from_logfile(float(a),float(b),float(c),
                                                                                 math.radians(float(alpha)),
                                                                                 math.radians(float(beta)),
                                                                                 math.radians(float(gamma)),
                                                                                 self.aimless['DataProcessingLattice']))
                    try:
                        high_symmetry_boost=self.nr_asu_in_unitcell_for_point_group[self.aimless['DataProcessingPointGroup']]
                        self.aimless['DataProcessingScore'] = (float(self.aimless['DataProcessingUniqueReflectionsOverall'])*\
                                                               float(self.aimless['DataProcessingCompletenessOverall'])*\
                                                               high_symmetry_boost*\
                                                               float(self.aimless['DataProcessingIsigOverall']))/float(self.aimless['DataProcessingUnitCellVolume'])
                    except ValueError:
                        self.aimless['DataProcessingScore']=0.0
        self.aimless['DataProcessingUnitCell']=str(a)+' '+str(b)+' '+str(c)+' '+str(alpha)+' '+str(beta)+' '+str(gamma)
        self.aimless['DataProcessingResolutionOverall']=str(self.aimless['DataProcessingResolutionLow'])+' - '+str(self.aimless['DataProcessingResolutionHigh'])

        # Hex Color code:
        # red:      #FF0000
        # orange:   #FF9900
        # green:    #00FF00
        # gray:     #E0E0E0

        if self.aimless['DataProcessingResolutionHigh']=='n/a' or self.aimless['DataProcessingRmergeLow'] =='n/a':
            self.aimless['DataProcessingAlert'] = '#FF0000'
        else:
            if float(self.aimless['DataProcessingResolutionHigh']) > 3.5 or float(self.aimless['DataProcessingRmergeLow']) > 0.1:
                self.aimless['DataProcessingAlert'] = '#FF0000'
            if (float(self.aimless['DataProcessingResolutionHigh']) <= 3.5 and float(self.aimless['DataProcessingResolutionHigh']) > 2.5) or \
               (float(self.aimless['DataProcessingRmergeLow']) <= 0.1 and float(self.aimless['DataProcessingRmergeLow']) > 0.05):
                self.aimless['DataProcessingAlert'] = '#FF9900'
            if float(self.aimless['DataProcessingResolutionHigh']) <= 2.5 and float(self.aimless['DataProcessingRmergeLow']) <= 0.05:
                self.aimless['DataProcessingAlert'] = '#00FF00'

        return self.aimless




    def get_lattice_from_space_group(self,logfile_spacegroup):
        lattice_type='n/a'
        for lattice in self.space_group_dict:
            for spacegroup in self.space_group_dict[lattice]:
                if logfile_spacegroup.replace(' ','')==spacegroup:
                    lattice_type=lattice
                    break
        return lattice_type

    def get_pointgroup_from_space_group(self,logfile_spacegroup):
        pointgroup='n/a'
        for pg in self.point_group_dict:
            for spacegroup in self.point_group_dict[pg]:
                if logfile_spacegroup.replace(' ','')==spacegroup:
                    pointgroup=pg
                    break
        return pointgroup



    def calc_unitcell_volume_from_logfile(self,a,b,c,alpha,beta,gamma,lattice):
        unitcell_volume=0
        if lattice=='triclinic':
            unitcell_volume=a*b*c* \
                            math.sqrt((1-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2) \
                                      +2*(math.cos(alpha)*math.cos(beta)*math.cos(gamma)))
#        if lattice=='monoclinic':
        if 'monoclinic' in lattice:
            unitcell_volume=round(a*b*c*math.sin(beta),1)
        if lattice=='orthorhombic' or lattice=='tetragonal' or lattice=='cubic':
            unitcell_volume=round(a*b*c,1)
        if lattice=='hexagonal' or lattice=='rhombohedral':
            unitcell_volume=round(a*b*c*(math.sin(math.radians(60))),1)
        return unitcell_volume

#    def get_all_values_as_dict(self):
#        info = {    'unitcell':         'n/a',
#                    'spacegroup':       'n/a',
#                    'unitcell_volume':  'n/a',
#                    'bravais_lattice':  'n/a'   }





    def PDBheader(self,pdbfile):
        PDBinfo = { 'Rcryst':           'n/a',
                    'RcrystTL':         'gray',
                    'Rfree':            'n/a',
                    'RfreeTL':          'gray',
                    'SpaceGroup':       'n/a',
                    'PointGroup':       'n/a',
                    'UnitCell':         'n/a',
                    'ResolutionHigh':   'n/a',
                    'ResolutionColor':  'gray',
                    'Lattice':          'n/a',
                    'UnitCellVolume':       0,
                    'Alert':            '#E0E0E0',
                    'rmsdBonds':        'n/a',
                    'rmsdBondsTL':      'gray',
                    'rmsdAngles':       'n/a',
                    'rmsdAnglesTL':     'gray'}

        a='n/a'
        b='n/a'
        c='n/a'
        alpha='n/a'
        beta='n/a'
        gamma='n/a'

        if os.path.isfile(pdbfile):
            for line in open(pdbfile):
                try:
                    if line.startswith('REMARK   3   R VALUE     (WORKING + TEST SET) :'):
                        PDBinfo['Rcryst']=line.split()[9]
                        if float(PDBinfo['Rcryst']) > 0.4:
                            PDBinfo['Alert']='#FF0000'
                            PDBinfo['RcrystTL'] = 'red'
                        if float(PDBinfo['Rcryst']) <= 0.4 and float(PDBinfo['Rcryst']) >= 0.3:
                            PDBinfo['Alert']='#FF9900'
                            PDBinfo['RcrystTL'] = 'orange'
                        if float(PDBinfo['Rcryst']) < 0.3:
                            PDBinfo['Alert']='#00FF00'
                            PDBinfo['RcrystTL'] = 'green'
                    if line.startswith('REMARK   3   FREE R VALUE                     :'):
                        PDBinfo['Rfree']=line.split()[6]
                        if float(PDBinfo['Rfree']) > 0.4:
                            PDBinfo['Alert']='#FF0000'
                            PDBinfo['RfreeTL'] = 'red'
                        if float(PDBinfo['Rfree']) <= 0.4 and float(PDBinfo['Rfree']) >= 0.3:
                            PDBinfo['Alert']='#FF9900'
                            PDBinfo['RfreeTL'] = 'orange'
                        if float(PDBinfo['Rfree']) < 0.3:
                            PDBinfo['Alert']='#00FF00'
                            PDBinfo['RfreeTL'] = 'green'
                except ValueError:
                    pass

                if line.startswith('REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) :'):
                    PDBinfo['ResolutionHigh']=line.split()[7]
                    try:
                        if float(line.split()[7]) < 2.4:                                   PDBinfo['ResolutionColor'] = 'green'
                        if float(line.split()[7]) >= 2.4 and float(line.split()[7]) < 2.8: PDBinfo['ResolutionColor'] = 'orange'
                        if float(line.split()[7]) >= 2.8:                                  PDBinfo['ResolutionColor'] = 'red'
                    except ValueError:
                        pass
                if line.startswith('REMARK   3   BOND LENGTHS REFINED ATOMS        (A):'):
                    PDBinfo['rmsdBonds'] = line.split()[9]
                    try:
                        if float(line.split()[9]) < 0.02:                                       PDBinfo['rmsdBondsTL'] = 'green'
                        if float(line.split()[9]) >= 0.02 and float(line.split()[9]) < 0.03:    PDBinfo['rmsdBondsTL'] = 'orange'
                        if float(line.split()[9]) >= 0.03:                                      PDBinfo['rmsdBondsTL'] = 'red'
                    except ValueError:
                        pass
                if line.startswith('REMARK   3   BOND ANGLES REFINED ATOMS   (DEGREES):'):
                    PDBinfo['rmsdAngles'] = line.split()[9]
                    try:
                        if float(line.split()[9]) < 2.0:                                        PDBinfo['rmsdAnglesTL'] = 'green'
                        if float(line.split()[9]) >= 2.0 and float(line.split()[9]) < 3.0:      PDBinfo['rmsdAnglesTL'] = 'orange'
                        if float(line.split()[9]) >= 3.0:                                       PDBinfo['rmsdAnglesTL'] = 'red'
                    except ValueError:
                        pass

                if line.startswith('CRYST1'):
                    a=int(float(line.split()[1]))
                    b=int(float(line.split()[2]))
                    c=int(float(line.split()[3]))
                    alpha=int(float(line.split()[4]))
                    beta=int(float(line.split()[5]))
                    gamma=int(float(line.split()[6]))

                    PDBinfo['UnitCell']=line.split()[1]+' '+line.split()[2]+' '+line.split()[3]+' '+ \
                                        line.split()[4]+' '+line.split()[5]+' '+line.split()[6]
#                    PDBinfo['SpaceGroup']=line[55:len(line)-1].replace(' ','').rstrip('\r')
                    PDBinfo['SpaceGroup']=str(line[55:65]).rstrip()

                    PDBinfo['Lattice']=self.get_lattice_from_space_group(PDBinfo['SpaceGroup'])
                    PDBinfo['PointGroup']=self.get_pointgroup_from_space_group(PDBinfo['SpaceGroup'])
                    if a != 'n/a' and b != 'n/a' and c != 'n/a' and \
                     alpha != 'n/a' and beta != 'n/a' and gamma != 'n/a' and PDBinfo['Lattice'] != 'n/a':
                        PDBinfo['UnitCellVolume']=self.calc_unitcell_volume_from_logfile(float(a),float(b),float(c),
                                                                                 math.radians(float(alpha)),
                                                                                 math.radians(float(beta)),
                                                                                 math.radians(float(gamma)),
                                                                                 PDBinfo['Lattice'])



        return PDBinfo


    def dict_for_datasource_update(self,pdbfile):
        db_dict={}
        db_dict['RefinementPDB_latest']=os.path.realpath(pdbfile)
        pdb=self.PDBheader(pdbfile)
        db_dict['RefinementRcryst'] =               pdb['Rcryst']
        db_dict['RefinementRcrystTraficLight'] =    pdb['RcrystTL']
        db_dict['RefinementRfree']=                 pdb['Rfree']
        db_dict['RefinementRfreeTraficLight'] =     pdb['RfreeTL']
        db_dict['RefinementRmsdBonds']  =           pdb['rmsdBonds']
        db_dict['RefinementRmsdBondsTL'] =          pdb['rmsdBondsTL']
        db_dict['RefinementRmsdAngles'] =           pdb['rmsdAngles']
        db_dict['RefinementRmsdAnglesTL'] =         pdb['rmsdAnglesTL']
        db_dict['RefinementSpaceGroup'] =           pdb['SpaceGroup']
        db_dict['RefinementResolution'] =           pdb['ResolutionHigh']
        db_dict['RefinementResolutionTL'] =         pdb['ResolutionColor']
        return db_dict



    def update_datasource_with_PDBheader(self,xtal,datasource,pdbfile):
        db_dict={}
        db_dict['RefinementPDB_latest']=os.path.realpath(pdbfile)
        pdb=self.PDBheader(pdbfile)
        db_dict['RefinementRcryst'] =               pdb['Rcryst']
        db_dict['RefinementRcrystTraficLight'] =    pdb['RcrystTL']
        db_dict['RefinementRfree']=                 pdb['Rfree']
        db_dict['RefinementRfreeTraficLight'] =     pdb['RfreeTL']
        db_dict['RefinementRmsdBonds']  =           pdb['rmsdBonds']
        db_dict['RefinementRmsdBondsTL'] =          pdb['rmsdBondsTL']
        db_dict['RefinementRmsdAngles'] =           pdb['rmsdAngles']
        db_dict['RefinementRmsdAnglesTL'] =         pdb['rmsdAnglesTL']
        db_dict['RefinementSpaceGroup'] =           pdb['SpaceGroup']
        db_dict['RefinementResolution'] =           pdb['ResolutionHigh']
        db_dict['RefinementResolutionTL'] =         pdb['ResolutionColor']
        print db_dict
        db=XChemDB.data_source(datasource)
        db.update_data_source(xtal,db_dict)

    def update_datasource_with_phenix_validation_summary(self,xtal,datasource,validation_summary):
        db_dict={}
        if os.path.isfile(validation_summary):
            for line in open(validation_summary):
                if 'molprobity score' in line.lower():
                    if len(line.split()) >= 4:
                        db_dict['RefinementMolProbityScore'] = line.split()[3]
                        if float(line.split()[3]) < 2:
                            db_dict['RefinementMolProbityScoreTL'] = 'green'
                        if float(line.split()[3]) >= 2 and float(line.split()[3]) < 3:
                            db_dict['RefinementMolProbityScoreTL'] = 'orange'
                        if float(line.split()[3]) >= 3:
                            db_dict['RefinementMolProbityScoreTL'] = 'red'

                if 'ramachandran outliers' in line.lower():
                    if len(line.split()) >= 4:
                        db_dict['RefinementRamachandranOutliers'] = line.split()[3]
                        if float(line.split()[3]) < 0.3:
                            db_dict['RefinementRamachandranOutliersTL'] = 'green'
                        if float(line.split()[3]) >= 0.3 and float(line.split()[3]) < 1:
                            db_dict['RefinementRamachandranOutliersTL'] = 'orange'
                        if float(line.split()[3]) >= 1:
                            db_dict['RefinementRamachandranOutliersTL'] = 'red'

                if 'favored' in line.lower():
                    if len(line.split()) >= 3:
                        db_dict['RefinementRamachandranFavored'] = line.split()[2]
                        if float(line.split()[2]) < 90:
                            db_dict['RefinementRamachandranFavoredTL'] = 'red'
                        if float(line.split()[2]) >= 90 and float(line.split()[2]) < 98:
                            db_dict['RefinementRamachandranFavoredTL'] = 'orange'
                        if float(line.split()[2]) >= 98:
                            db_dict['RefinementRamachandranFavoredTL'] = 'green'
        else:
            db_dict['RefinementMolProbityScore']        = '-'
            db_dict['RefinementMolProbityScoreTL']      = 'gray'
            db_dict['RefinementRamachandranOutliers']   = '-'
            db_dict['RefinementRamachandranOutliersTL'] = 'gray'
            db_dict['RefinementRamachandranFavored']    = '-'
            db_dict['RefinementRamachandranFavoredTL']  = 'gray'
        db=XChemDB.data_source(datasource)
        db.update_data_source(xtal,db_dict)

class mtztools:

    def __init__(self,mtzfile):
        self.mtzfile=mtzfile
        self.hkl = any_reflection_file(file_name=self.mtzfile)
        self.miller_arrays = self.hkl.as_miller_arrays()
        self.mtz = self.miller_arrays[0]
        self.iotbxMTZ = mtz.object(self.mtzfile)

        self.space_group_dict=   {  'triclinic':    [1],
#                                    'monoclinic':   [3,4,5],
                                    'monoclinic_P': [3,4],
                                    'monoclinic_C': [5],
                                    'orthorhombic': [16,17,18,19,20,21,22,23,24],
                                    'tetragonal':   [75,76,77,78,79,80,89,90,91,92,93,94,95,96,97,98],
                                    'hexagonal':    [143,144,145,149,150,151,152,153,154,168,169,170,
                                                     171,172,173,177,178,179,180,181,182],
                                    'rhombohedral': [146,155],
                                    'cubic':        [195,196,197,198,199,
                                                     207,208,209,210,211,212,213,214]  }

        self.translate_spg_to_number_dict = {
            'p1': 1, 'p121': 3, 'p1211': 4, 'c2': 5, 'c121': 5, 'i2': 5, 'i121': 5, 'p222': 16,
            'p2122': 17, 'p2212': 17, 'p2221': 17, 'p21212': 18, 'p21221': 18, 'p22121': 18,
            'p212121': 19, 'c2221': 20, 'c222': 21, 'f222': 22, 'i222': 23, 'i212121': 24,
            'p4': 75, 'p41': 76, 'p42': 77, 'p43': 78, 'i4': 79, 'i41': 80,
            'p422': 89, 'p4212': 90, 'p4122': 91, 'p41212': 92, 'p4222': 93, 'p42212': 94,
            'p4322': 95, 'p43212': 96, 'i422': 97, 'i4122': 98,
            'p3': 143, 'p31': 144, 'p32': 145, 'p312': 149, 'p321': 150, 'p3112': 151, 'p3121': 152,
            'p3212': 153, 'p3221': 154, 'p6': 168, 'p61': 169, 'p65': 170, 'p62': 171, 'p64': 172, 'p63': 173,
            'p622': 177, 'p6122': 178, 'p6522': 179, 'p6222': 180, 'p6422': 181, 'p6322': 182,
            'r3': 146, 'h3': 146, 'r32': 155, 'h32': 155,
            'p23': 195, 'f23': 196, 'i23': 197, 'p213': 198, 'i213': 199,
            'p432': 207, 'p4232': 208, 'f432': 209, 'f4132': 210, 'i432': 211, 'p4332': 212,
            'p4132': 213, 'i4132': 214
        }


        self.point_group_dict=   {  '1':            [1],
                                    '2':            [3,4,5],
                                    '222':          [16,17,18,19,20,21,22,23,24],
                                    '4':            [75,76,77,78,79,80],
                                    '422':          [89,90,91,92,93,94,95,96,97,98],
                                    '3':            [143,144,145,146],
                                    '32':           [149,150,151,152,153,154,155],
                                    '6':            [168,169,170,171,172,173],
                                    '622':          [177,178,179,180,181,182],
                                    '23':           [195,196,197,198,199],
                                    '432':          [207,208,209,210,211,212,213,214]  }

        self.nr_asu_in_unitcell_for_point_group =   {   '1':            1,
                                                        '2':            2,
                                                        '222':          4,
                                                        '4':            4,
                                                        '422':          8,
                                                        '3':            3,
                                                        '32':           6,
                                                        '6':            6,
                                                        '622':          12,
                                                        '23':           12,
                                                        '432':          24  }

        self.aimless = {    'DataProcessingProgram':                        'n/a',
                            'DataCollectionRun':                            'n/a',
                            'DataProcessingSpaceGroup':                     'n/a',
                            'DataProcessingUnitCell':                       'n/a',
                            'DataProcessingA':                              'n/a',
                            'DataProcessingB':                              'n/a',
                            'DataProcessingC':                              'n/a',
                            'DataProcessingAlpha':                          'n/a',
                            'DataProcessingBeta':                           'n/a',
                            'DataProcessingGamma':                          'n/a',
                            'DataProcessingResolutionLow':                  'n/a',
                            'DataProcessingResolutionLowInnerShell':        'n/a',
                            'DataProcessingResolutionHigh':                 'n/a',
                            'DataProcessingResolutionHighOuterShell':       'n/a',
                            'DataProcessingResolutionOverall':              'n/a',
                            'DataProcessingRmergeOverall':                  'n/a',
                            'DataProcessingRmergeLow':                      'n/a',
                            'DataProcessingRmergeHigh':                     'n/a',
                            'DataProcessingIsigOverall':                    'n/a',
                            'DataProcessingIsigLow':                        'n/a',
                            'DataProcessingIsigHigh':                       'n/a',
                            'DataProcessingCompletenessOverall':            'n/a',
                            'DataProcessingCompletenessLow':                'n/a',
                            'DataProcessingCompletenessHigh':               'n/a',
                            'DataProcessingMultiplicityOverall':            'n/a',
                            'DataProcessingMultiplicityLow':                'n/a',
                            'DataProcessingMultiplicityHigh':               'n/a',
                            'DataProcessingCChalfOverall':                  'n/a',
                            'DataProcessingCChalfLow':                      'n/a',
                            'DataProcessingCChalfHigh':                     'n/a',
                            'DataProcessingResolutionHigh15sigma':         'n/a',
                            'DataProcessingUniqueReflectionsOverall':       'n/a',
                            'DataProcessingLattice':                        'n/a',
                            'DataProcessingPointGroup':                     'n/a',
                            'DataProcessingUnitCellVolume':                 0,
                            'DataProcessingAlert':                          '#FF0000',
                            'DataCollectionWavelength':                     'n/a',
                            'DataProcessingScore':                          'n/a'             }

    def get_dmin(self):
        return str(round(float(self.mtz.d_min()), 2))

    def get_wavelength(self):
        wavelength = 0.0
        for crystal in self.iotbxMTZ.crystals():
            for dataset in crystal.datasets():
                if not dataset.wavelength() == 0.0:
                    wavelength = str(round(dataset.wavelength(),5))
                    break
        return wavelength

    def get_information_for_datasource(self):
        db_dict={}
        mtz_dict=self.get_all_values_as_dict()
        pg=self.get_pointgroup_from_mtz()
        if mtz_dict != {}:
            db_dict['DataProcessingResolutionHigh']=mtz_dict['resolution_high']
            db_dict['DataProcessingUnitCell']=mtz_dict['unitcell']
            db_dict['DataProcessingSpaceGroup']=mtz_dict['spacegroup']
            db_dict['DataProcessingUnitCellVolume']=mtz_dict['unitcell_volume']
            db_dict['DataProcessingLattice']=mtz_dict['bravais_lattice']
        if pg != '':
            db_dict['DataProcessingPointGroup']=pg
        return db_dict

    def get_bravais_lattice_from_spg_number(self,number):
        lattice=''
        for bravaislattice in self.space_group_dict:
            for spg_number in self.space_group_dict[bravaislattice]:
                if spg_number==number:
                    lattice=bravaislattice
        return lattice

    def get_point_group_from_spg_number(self,number):
        pointgroup=''
        for pg in self.point_group_dict:
            for spg_number in self.point_group_dict[pg]:
                if spg_number==number:
                    pointgroup=pg
        return pointgroup


    def get_unit_cell_from_mtz(self):
        unitcell=[]
        cell_line=100000
        a=0
        b=0
        c=0
        alpha=0
        beta=0
        gamma=0
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' * Cell Dimensions :'):
                cell_line=n+2
            if n==cell_line and len(line.split())==6:
                a=line.split()[0]
                b=line.split()[1]
                c=line.split()[2]
                alpha=line.split()[3]
                beta=line.split()[4]
                gamma=line.split()[5]
        return [a,b,c,alpha,beta,gamma]

    def get_spg_number_from_mtz(self):
        spg_number=0
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' * Space group ='):
                spg_number=int(line[line.rfind(')')-3:line.rfind(')')])
        return spg_number

    def get_spg_from_mtz(self):
        spg=''
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' * Space group ='):
                spg=line[line.find("'")+1:line.rfind("'")]
        return spg

    def get_bravais_lattice_from_mtz(self):
        bravais_lattice=''
        spg_number=self.get_spg_number_from_mtz()
        bravais_lattice=self.get_bravais_lattice_from_spg_number(spg_number)
        return bravais_lattice

    def get_pointgroup_from_mtz(self):
        pointgroup=''
        spg_number=self.get_spg_number_from_mtz()
        pointgroup=self.get_point_group_from_spg_number(spg_number)
        return pointgroup

    def calc_unitcell_volume_from_mtz(self):
        spg_number=self.get_spg_number_from_mtz()
        lattice=self.get_bravais_lattice_from_spg_number(spg_number)
        unitcell=self.get_unit_cell_from_mtz()
        a=float(unitcell[0])
        b=float(unitcell[1])
        c=float(unitcell[2])
        alpha=math.radians(float(unitcell[3]))
        beta= math.radians(float(unitcell[4]))
        gamma=math.radians(float(unitcell[5]))
        unitcell_volume=0
        if lattice=='triclinic':
            unitcell_volume=a*b*c* \
                            math.sqrt((1-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2) \
                                      +2*(math.cos(alpha)*math.cos(beta)*math.cos(gamma)))
        if 'monoclinic' in lattice:
            unitcell_volume=round(a*b*c*math.sin(beta),1)
        if lattice=='orthorhombic' or lattice=='tetragonal' or lattice=='cubic':
            unitcell_volume=round(a*b*c,1)
        if lattice=='hexagonal' or lattice=='rhombohedral':
            unitcell_volume=round(a*b*c*(math.sin(math.radians(60))),1)
        return unitcell_volume

    def get_high_resolution_from_mtz(self):
        resolution_high='n/a'
        resolution_line=1000000
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' *  Resolution Range :'):
                resolution_line=n+2
            if n==resolution_line and len(line.split())==8:
                resolution_high=line.split()[5]
        return resolution_high

    def get_low_resolution_from_mtz(self):
        resolution_low='n/a'
        resolution_line=1000000
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' *  Resolution Range :'):
                resolution_line=n+2
            if n==resolution_line and len(line.split())==8:
                resolution_low=line.split()[3]
        return resolution_low

    def get_number_measured_reflections(self):
        missing_reflections='0'
        all_reflections='0'
        meassured_reflections='0'
        resolution_line=1000000
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        foundTable=False
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' Col Sort    Min    Max    Num      %     Mean     Mean   Resolution   Type Column'):
                foundTable=True
            if foundTable and len(line.split())==12:
                if line.split()[11]=='F':
                    missing_reflections=line.split()[4]
                    foundTable=False
            if line.startswith(' No. of reflections used in FILE STATISTICS'):
                all_reflections=line.split()[7]
                break
        try:
            meassured_reflections=int(all_reflections)-int(missing_reflections)
        except ValueError:
            pass
        return meassured_reflections

    def calculate_correlaton_between_intensities_in_mtzfiles(self,mtzin):
        CC = '0.0'
        errorMessage=''
#        cmd = ( 'pointless hklin %s hklref %s << eof\n' %(mtzin,self.mtzfile)+
#                'labref F=F\n'
#                'labin F=F\n'
#                'eof\n' )
        cmd = ( 'pointless hklin %s hklref %s << eof\n' %(mtzin,self.mtzfile)+
                'labref I=IMEAN\n'
                'labin I=IMEAN\n'
                'eof\n' )

        pointless=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        foundLine=False
        for line in iter(pointless.stdout.readline,''):
            if foundLine:
                CC=line.split()[3]
                break
            if 'Alternative reindexing        Lklhd      CC' in line:
                foundLine=True
            if '**** Incompatible symmetries ****' in line:
                errorMessage='**** Incompatible symmetries ****'
                break
            if 'Merged test dataset (HKLIN) has different Laue symmetry to reference set' in line:
                errorMessage='%s has different Laue symmetry to %s' %(mtzin,self.mtzfile)
                break

        return CC,errorMessage

    def get_all_values_as_dict(self):
        mtz = { 'resolution_high':  'n/a',
                'unitcell':         'n/a',
                'spacegroup':       'n/a',
                'unitcell_volume':  'n/a',
                'bravais_lattice':  'n/a'   }
        a=0.0
        b=0.0
        c=0.0
        alpha_rad=0.0
        beta_rad=0.0
        gamma_rad=0.0
        resolution_line=1000000
        cell_line=1000000
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' *  Resolution Range :'):
                resolution_line=n+2
            if n==resolution_line and len(line.split())==8:
                mtz['resolution_high']=round(float(line.split()[5]),2)
            if line.startswith(' * Cell Dimensions :'):
                cell_line=n+2
            if n==cell_line and len(line.split())==6:
                a=      round(float(line.split()[0]),1)
                b=      round(float(line.split()[1]),1)
                c=      round(float(line.split()[2]),1)
                alpha=  round(float(line.split()[3]),1)
                beta=   round(float(line.split()[4]),1)
                gamma=  round(float(line.split()[5]),1)
                mtz['unitcell']=str(a)+' '+str(b)+' '+str(c)+' '+ \
                                str(alpha)+' '+str(beta)+' '+str(gamma)
                alpha_rad=math.radians(alpha)
                beta_rad= math.radians(beta)
                gamma_rad=math.radians(gamma)
            if line.startswith(' * Space group ='):
                spg_number=int(line[line.rfind(')')-3:line.rfind(')')])
                mtz['bravais_lattice']=self.get_bravais_lattice_from_spg_number(spg_number)
                mtz['spacegroup']=line[line.find("'")+1:line.rfind("'")]
        if mtz['bravais_lattice']=='triclinic':
            mtz['unitcell_volume']=a*b*c* \
                            math.sqrt((1-math.cos(alpha_rad)**2-math.cos(beta_rad)**2-math.cos(gamma_rad)**2) \
                                      +2*(math.cos(alpha_rad)*math.cos(beta_rad)*math.cos(gamma_rad)))
#        elif mtz['bravais_lattice']=='monoclinic':
        elif 'monoclinic' in mtz['bravais_lattice']:
            mtz['unitcell_volume']=round(a*b*c*math.sin(beta_rad),1)
        elif mtz['bravais_lattice']=='orthorhombic' or mtz['bravais_lattice']=='tetragonal' or mtz['bravais_lattice']=='cubic':
            mtz['unitcell_volume']=round(a*b*c,1)
        elif mtz['bravais_lattice']=='hexagonal' or mtz['bravais_lattice']=='rhombohedral':
            mtz['unitcell_volume']=round(a*b*c*(math.sin(math.radians(60))),1)

        return mtz

    def get_all_columns_as_dict(self):
        column_dict = { 'F':        [],
                        'I':        [],
                        'SIG':      [],
                        'PHS':      [],
                        'FOM':      [],
                        'RFREE':    []  }
        startline=1000000
        mtzdmp=subprocess.Popen(['mtzdmp',self.mtzfile],stdout=subprocess.PIPE)
        for n,line in enumerate(iter(mtzdmp.stdout.readline,'')):
            if line.startswith(' Col Sort    Min    Max    Num'):
                startline=n+2
            if n >= startline and len(line.split()) > 10:
                if line.split()[10] == 'F':
                    column_dict['F'].append(line.split()[11])
                if line.split()[10] == 'J':
                    column_dict['I'].append(line.split()[11])
                if line.split()[10] == 'Q':
                    column_dict['SIG'].append(line.split()[11])
                if line.split()[10] == 'I':
                    column_dict['RFREE'].append(line.split()[11])
                if line.split()[10] == 'P':
                    column_dict['PHS'].append(line.split()[11])
                if line.split()[10] == 'W':
                    column_dict['FOM'].append(line.split()[11])

        return column_dict


class queue:

    def jobs_in_queue(self):
            jobs=[]
            qstat=subprocess.Popen(['qstat -r'],stdout=subprocess.PIPE)
            counter=0
            for n,line in enumerate(iter(qstat.stdout.readline,'')):
                jobID=''
                job_name=''
                user_name=''
                if getpass.getuser() in line:
                    counter +=1
                    jobID=line.split()[0]
                    user_name=line.split()[3]
                if "Full jobname:" in line:
                    job_name=line.split()[2]
                jobs.append([counter,jobID,True,job_name,user_name])
            return jobs

    def remove_all_jobs_from_queue(self):
            jobs_in_queue=self.jobs_in_queue()
            if jobs_in_queue != []:
                for job in jobs_in_queue:
                    os.system('qdel '+job[0])



class external_software:

    def __init__(self,xce_logfile):

        self.available_programs = {}
        self.Logfile=XChemLog.updateLog(xce_logfile)

    def check(self):

        self.Logfile.insert('Searching for external software...')

        # default is False; user needs to explicitely set this
        self.available_programs['qsub_remote']=''

        FNULL = open(os.devnull, 'w')

        try:
#            subprocess.call(['qstat'], stdout=FNULL, stderr=subprocess.STDOUT)
            p = subprocess.Popen('qstat', shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
            status='found'
            try:
                for line in p.stdout:
                    if 'symbol lookup error:' in line or 'command not found' in line or 'error' in line:
                        self.available_programs['qsub']=False
                        self.available_programs['qsub_array']=False
                        status='not found'
                        array_status='not found'
            except IOError:
                self.available_programs['qsub']=False
                self.available_programs['qsub_array']=False
                status='not found'
                array_status='not found'

            if status == 'found':
                self.available_programs['qsub']=True
#                if os.getcwd().startswith('/dls'):
                if os.path.isdir('/dls'):
                    self.available_programs['qsub_array']=True
                    array_status='found'
                else:
                    self.available_programs['qsub_array']=False
                    array_status='not found'
        except OSError:
            self.available_programs['qsub']=False
            self.available_programs['qsub_array']=False
            status='not found'
            array_status='not found'
        self.Logfile.insert('{0:50} {1:10}'.format('checking for qsub:', status))
#        print '{0:50} {1:10}'.format('-checking for qsub:', status)
#        if os.getcwd().startswith('/dls'):
        if os.path.isdir('/dls'):
            self.Logfile.insert('{0:50} {1:10}'.format('checking for array qsub:', array_status))

        try:
            subprocess.call(['refmac5','end'], stdout=FNULL, stderr=subprocess.STDOUT)
            self.available_programs['refmac5']=True
            status='found'
        except OSError:
            self.available_programs['refmac5']=False
            status='not found'
        self.Logfile.insert('{0:50} {1:10}'.format('checking for refmac5:', status))

        try:
            subprocess.call(['phenix.molprobity'], stdout=FNULL, stderr=subprocess.STDOUT)
            self.available_programs['phenix.molprobity']=True
            status='found'
        except OSError:
            self.available_programs['phenix.molprobity']=False
            status='not found'
        self.Logfile.insert('{0:50} {1:10}'.format('checking for phenix.molprobity:', status))

        try:
            subprocess.call(['phenix.find_tls_groups'], stdout=FNULL, stderr=subprocess.STDOUT)
            self.available_programs['phenix.find_tls_groups']=True
            status='found'
        except OSError:
            self.available_programs['phenix.find_tls_groups']=False
            status='not found'
        self.Logfile.insert('{0:50} {1:10}'.format('checking for phenix.find_tls_groups:', status))

        try:
            subprocess.call(['mmtbx.validate_ligands'], stdout=FNULL, stderr=subprocess.STDOUT)
            self.available_programs['mmtbx.validate_ligands']=True
            status='found'
        except OSError:
            self.available_programs['mmtbx.validate_ligands']=False
            status='not found'
        self.Logfile.insert('{0:50} {1:10}'.format('checking for mmtbx.validate_ligands:', status))

        try:
            subprocess.call(['acedrg'], stdout=FNULL, stderr=subprocess.STDOUT)
            self.available_programs['acedrg']=True
            shutil.rmtree('AcedrgOut_TMP')
            status='found'
        except OSError:
            self.available_programs['acedrg']=False
            status='not found'
        self.Logfile.insert('{0:50} {1:10}'.format('checking for acedrg:', status))

        try:
            subprocess.call(['phenix.elbow'], stdout=FNULL, stderr=subprocess.STDOUT)
            self.available_programs['phenix.elbow']=True
            status='found'
        except OSError:
            self.available_programs['phenix.elbow']=False
            status='not found'
        self.Logfile.insert('{0:50} {1:10}'.format('checking for phenix.elbow:', status))

        try:
            subprocess.call(['grade'], stdout=FNULL, stderr=subprocess.STDOUT)
            self.available_programs['grade']=True
            status='found'
        except OSError:
            self.available_programs['grade']=False
            status='not found'
        self.Logfile.insert('{0:50} {1:10}'.format('checking for grade:', status))

        try:
            subprocess.call(['giant.create_occupancy_params'], stdout=FNULL, stderr=subprocess.STDOUT)
            self.available_programs['giant.create_occupancy_params']=True
            status='found'
        except OSError:
            self.available_programs['giant.create_occupancy_params']=False
            status='not found'
        self.Logfile.insert('{0:50} {1:10}'.format('checking for giant.create_occupancy_params:', status))


        return self.available_programs


class ParseFiles:

    def __init__(self,DataPath,xtalID):
        # probably need to read in compoundID, because for custom projects, need to take newest pdb file
        # that has not same root as compoundID
        self.DataPath = DataPath
        self.xtalID = xtalID

        if '<samplename>' in DataPath:
            self.newestPDB = max(glob.iglob(self.DataPath.replace('<samplename>',xtalID)+'/*.pdb'), key=os.path.getctime)
        else:
            self.newestPDB = self.DataPath+'/'+self.xtalID+'/refine.pdb'

    def UpdateQualityIndicators(self):

        QualityIndicators = {   'R':                            '-',
                                'RRfree':                       '-',
                                'RRfreeColor':                  'gray',
                                'Resolution':                   'n/a',
                                'ResolutionColor':              'gray',
                                'MolprobityScore':              'n/a',
                                'MolprobityScoreColor':         'gray',
                                'RamachandranOutliers':         'n/a',
                                'RamachandranOutliersColor':    'gray',
                                'RamachandranFavored':          'n/a',
                                'RamachandranFavoredColor':     'gray',
                                'LigandCCcolor':                'gray',
                                'LigandCC':                     'n/a',
                                'rmsdBonds':                    'n/a',
                                'rmsdBondsColor':               'gray',
                                'rmsdAngles':                   'n/a',
                                'rmsdAnglesColor':              'gray',
                                'MatrixWeight':                 'n/a'   }

        # R, Rfree, Resolution
        found = 0
        if os.path.isfile(self.newestPDB):
            for line in open(self.newestPDB):
                if line.startswith('REMARK   3   R VALUE     (WORKING + TEST SET) :'):
                    QualityIndicators['R'] = line.split()[9]
                if line.startswith('REMARK   3   FREE R VALUE                     :'):
                    QualityIndicators['RRfree'] = line.split()[6]
                    if float(line.split()[6]) < 0.3:                             QualityIndicators['RRfreeColor'] = 'green'
                    if float(line.split()[6]) >= 0.3 and float(line.split()[6]): QualityIndicators['RRfreeColor'] = 'orange'
                    if float(line.split()[6]) >= 0.4:                            QualityIndicators['RRfreeColor'] = 'red'
                if line.startswith('REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) :'):
                    QualityIndicators['Resolution'] = line.split()[7]
                    if float(line.split()[7]) < 2.4:                                   QualityIndicators['ResolutionColor'] = 'green'
                    if float(line.split()[7]) >= 2.4 and float(line.split()[7]) < 2.8: QualityIndicators['ResolutionColor'] = 'orange'
                    if float(line.split()[7]) >= 2.8:                                  QualityIndicators['ResolutionColor'] = 'red'
                if line.startswith('REMARK   3   BOND LENGTHS REFINED ATOMS        (A):'):
                    QualityIndicators['rmsdBonds'] = line.split()[9]
                    if float(line.split()[9]) < 0.02:                                   QualityIndicators['rmsdBondsColor'] = 'green'
                    if float(line.split()[9]) >= 0.02 and float(line.split()[9]) < 0.03: QualityIndicators['rmsdBondsColor'] = 'orange'
                    if float(line.split()[9]) >= 0.03:                                  QualityIndicators['rmsdBondsColor'] = 'red'
                if line.startswith('REMARK   3   BOND ANGLES REFINED ATOMS   (DEGREES):'):
                    QualityIndicators['rmsdAngles'] = line.split()[9]
                    if float(line.split()[9]) < 2.0:                                   QualityIndicators['rmsdAnglesColor'] = 'green'
                    if float(line.split()[9]) >= 2.0 and float(line.split()[9]) < 3.0: QualityIndicators['rmsdAnglesColor'] = 'orange'
                    if float(line.split()[9]) >= 3.0:                                  QualityIndicators['rmsdAnglesColor'] = 'red'


        # Molprobity
        if os.path.isfile(self.DataPath+'/'+self.xtalID+'/validation_summary.txt'):
            for line in open(self.DataPath+'/'+self.xtalID+'/validation_summary.txt'):
#                if line.startswith('  Molprobity score      ='):
#                if line.lower().startswith('  molprobity score'):
                if 'molprobity score' in line.lower():
                    if len(line.split()) >= 4:
                        QualityIndicators['MolprobityScore'] = line.split()[3]
                        try:
                            if float(line.split()[3]) < 2:                                 QualityIndicators['MolprobityScoreColor'] = 'green'
                            if float(line.split()[3]) >= 2 and float(line.split()[3]) < 3: QualityIndicators['MolprobityScoreColor'] = 'orange'
                            if float(line.split()[3]) >= 3:                                QualityIndicators['MolprobityScoreColor'] = 'red'
                        except ValueError:
                            pass
#                if line.lower().startswith('  ramachandran outliers ='):
                if 'ramachandran outliers' in line.lower():
                    if len(line.split()) >= 4:
                        QualityIndicators['RamachandranOutliers'] = line.split()[3]
                        try:
                            if float(line.split()[3]) < 0.3:                                 QualityIndicators['RamachandranOutliersColor'] = 'green'
                            if float(line.split()[3]) >= 0.3 and float(line.split()[3]) < 1: QualityIndicators['RamachandranOutliersColor'] = 'orange'
                            if float(line.split()[3]) >= 1:                                  QualityIndicators['RamachandranOutliersColor'] = 'red'
                        except ValueError:
                            pass
#                if line.startswith('               favored  ='):
                if 'favored' in line.lower():
                    if len(line.split()) >= 3:
                        QualityIndicators['RamachandranFavored'] = line.split()[2]
                        try:
                            if float(line.split()[2]) < 90:                                  QualityIndicators['RamachandranFavoredColor'] = 'red'
                            if float(line.split()[2]) >= 90 and float(line.split()[2]) < 98: QualityIndicators['RamachandranFavoredColor'] = 'orange'
                            if float(line.split()[2]) >= 98:                                 QualityIndicators['RamachandranFavoredColor'] = 'green'
                        except ValueError:
                            pass

        # LigandCC
        if os.path.isfile(self.DataPath+'/'+self.xtalID+'/validate_ligands.txt'):
            for line in open(self.DataPath+'/'+self.xtalID+'/validate_ligands.txt'):
#                if line.startswith('|  LIG'):
                if 'LIG' in line:
                    QualityIndicators['LigandCC'] = line.split()[6]
                    if float(line.split()[6]) < 0.8:                                   QualityIndicators['LigandCCcolor'] = 'red'
                    if float(line.split()[6]) >= 0.8 and float(line.split()[6]) < 0.9: QualityIndicators['LigandCCcolor'] = 'orange'
                    if float(line.split()[6]) >= 0.9:                                  QualityIndicators['LigandCCcolor'] = 'green'

        # Matrix Weight
        temp = []
        found = 0
        if os.path.isdir(os.path.join(self.DataPath,self.xtalID)):
            for item in glob.glob(os.path.join(self.DataPath,self.xtalID,'*')):
                if item.startswith(os.path.join(self.DataPath,self.xtalID,'Refine_')):
                        temp.append(int(item[item.rfind('_')+1:]))
                        found = 1
        if found:
            Serial = max(temp)
            if os.path.isfile(os.path.join(self.DataPath,self.xtalID,'Refine_'+str(Serial),'refmac.log') ):
                for line in open(os.path.join(self.DataPath,self.xtalID,'Refine_'+str(Serial),'refmac.log')):
                    if line.startswith(' Weight matrix') and len(line.split())==3:
                        QualityIndicators['MatrixWeight']=line.split()[2]

        return QualityIndicators

class pdbtools(object):

    def __init__(self,pdb):
        self.pdb = pdb
        self.pdb_inp = iotbx.pdb.input(file_name=self.pdb)
        self.hierarchy = self.pdb_inp.construct_hierarchy()


        self.AminoAcids = ['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY','HIS',
                           'ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
        self.Solvents = ['DMS','EDO','GOL','HOH']
        self.Ions = ['NA','MG','CL','K','SO4','PO4','CA']
        self.AAdict = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q',
                       'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
                       'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
        self.xce_ligands = ['LIG','DRG','FRS']

        self.space_group_dict=   {  'triclinic':    [1],
#                                    'monoclinic':   [3,4,5],
                                    'monoclinic_P': [3,4],
                                    'monoclinic_C': [5],
                                    'orthorhombic': [16,17,18,19,20,21,22,23,24],
                                    'tetragonal':   [75,76,77,78,79,80,89,90,91,92,93,94,95,96,97,98],
                                    'hexagonal':    [143,144,145,149,150,151,152,153,154,168,169,170,
                                                     171,172,173,177,178,179,180,181,182],
                                    'rhombohedral': [146,155],
                                    'cubic':        [195,196,197,198,199,
                                                     207,208,209,210,211,212,213,214]  }

        self.translate_spg_to_number_dict = {
            'p1': 1, 'p121': 3, 'p1211': 4, 'c2': 5, 'c121': 5, 'i2': 5, 'i121': 5, 'p222': 16,
            'p2122': 17, 'p2212': 17, 'p2221': 17, 'p21212': 18, 'p21221': 18, 'p22121': 18,
            'p212121': 19, 'c2221': 20, 'c222': 21, 'f222': 22, 'i222': 23, 'i212121': 24,
            'p4': 75, 'p41': 76, 'p42': 77, 'p43': 78, 'i4': 79, 'i41': 80,
            'p422': 89, 'p4212': 90, 'p4122': 91, 'p41212': 92, 'p4222': 93, 'p42212': 94,
            'p4322': 95, 'p43212': 96, 'i422': 97, 'i4122': 98,
            'p3': 143, 'p31': 144, 'p32': 145, 'p312': 149, 'p321': 150, 'p3112': 151, 'p3121': 152,
            'p3212': 153, 'p3221': 154, 'p6': 168, 'p61': 169, 'p65': 170, 'p62': 171, 'p64': 172, 'p63': 173,
            'p622': 177, 'p6122': 178, 'p6522': 179, 'p6222': 180, 'p6422': 181, 'p6322': 182,
            'r3': 146, 'h3': 146, 'r32': 155, 'h32': 155,
            'p23': 195, 'f23': 196, 'i23': 197, 'p213': 198, 'i213': 199,
            'p432': 207, 'p4232': 208, 'f432': 209, 'f4132': 210, 'i432': 211, 'p4332': 212,
            'p4132': 213, 'i4132': 214
        }


        self.point_group_dict=   {  '1':            [1],
                                    '2':            [3,4,5],
                                    '222':          [16,17,18,19,20,21,22,23,24],
                                    '4':            [75,76,77,78,79,80],
                                    '422':          [89,90,91,92,93,94,95,96,97,98],
                                    '3':            [143,144,145,146],
                                    '32':           [149,150,151,152,153,154,155],
                                    '6':            [168,169,170,171,172,173],
                                    '622':          [177,178,179,180,181,182],
                                    '23':           [195,196,197,198,199],
                                    '432':          [207,208,209,210,211,212,213,214]  }

        self.nr_asu_in_unitcell_for_point_group =   {   '1':            1,
                                                        '2':            2,
                                                        '222':          4,
                                                        '4':            4,
                                                        '422':          8,
                                                        '3':            3,
                                                        '32':           6,
                                                        '6':            6,
                                                        '622':          12,
                                                        '23':           12,
                                                        '432':          24  }

    def get_refinement_program(self):
        program = 'unknown'
        for remark in self.pdb_inp.remark_section():
            if 'PROGRAM' in remark:
                if 'refmac' in remark.lower():
                    program = 'REFMAC'
                elif 'phenix' in remark.lower():
                    program = 'PHENIX'
        return  program

    def get_residues_with_resname(self,resname):
        ligands = []
        for model in self.hierarchy.models():
            for chain in model.chains():
                for conformer in chain.conformers():
                    for residue in conformer.residues():
                        if residue.resname == resname:
                            if [residue.resname, residue.resseq, chain.id] not in ligands:
                                ligands.append([residue.resname, residue.resseq, chain.id])
        return ligands

    def get_centre_of_gravity_of_residue(self,resname_resseq_chain):
        resname_x = resname_resseq_chain.split('-')[0]
        resseq_x = resname_resseq_chain.split('-')[1]
        chain_x = resname_resseq_chain.split('-')[2]
        print '===>',resname_x,resseq_x,chain_x
        x = []
        y = []
        z = []
        for model in self.hierarchy.models():
            for chain in model.chains():
                for conformer in chain.conformers():
                    for residue in conformer.residues():
                        print residue.resname,residue.resseq,chain.id
                        if residue.resname == resname_x and residue.resseq == resseq_x and chain.id == chain_x:
                            for atom in residue.atoms():
                                x.append(atom.xyz[0])
                                y.append(atom.xyz[1])
                                z.append(atom.xyz[2])

        if x != [] and y != [] and z != []:
            x = ((max(x) - min(x)) / 2) + min(x)
            y = ((max(y) - min(y)) / 2) + min(y)
            z = ((max(z) - min(z)) / 2) + min(z)

        return x,y,z

    def save_residues_with_resname(self,outDir,resname):
        ligands = self.get_residues_with_resname(resname)
        ligList = []
        for l in ligands:
            sel_cache = self.hierarchy.atom_selection_cache()
            lig_sel = sel_cache.selection("(resname %s and resseq %s and chain %s)" % (l[0], l[1], l[2]))
            hierarchy_lig = self.hierarchy.select(lig_sel)

            ligName = (l[0] + '-' + l[1] + '-' + l[2] + '.pdb').replace(' ', '')
            ligList.append(ligName)

            f = open(os.path.join(outDir,ligName), "w")
            f.write(hierarchy_lig.as_pdb_string(crystal_symmetry=self.pdb_inp.crystal_symmetry()))
            f.close()

        return ligList

    def GetRefinementProgram(self):
        program=''
        for line in open(self.pdb):
            if line.startswith('REMARK') and 'REFMAC' in line:
                program='REFMAC'
                break
            if line.startswith('REMARK') and 'PHENIX' in line:
                program='PHENIX'
                break
            if line.startswith('REMARK') and 'BUSTER' in line:
                program='BUSTER'
                break

        return program


    def GetSequence(self):
        chain = []
        Sequence=''
        # need to count residue numbers in case of alternative conformations
        ResiNum=[]
        for line in open(self.pdb):
            if line.startswith('ATOM') and line[17:20] in self.AminoAcids:
                if line[21:22] not in chain:
                    chain.append(line[21:22])
                    Sequence=Sequence+'\n>chain{0!s}.\n'.format(line[21:22])
                    ResiNum=[]
                if line[13:15]=='CA' and line[22:27] not in ResiNum:
                    Sequence=Sequence+self.AAdict[line[17:20]]
                    ResiNum.append(line[22:27])
        return Sequence


    def GetProteinChains(self):
        chain = []
        for line in open(self.pdb):
            if line.startswith('ATOM'):
                if line[17:20] in self.AminoAcids:
                    if line[21:22] not in chain: chain.append(line[21:22])
        return chain


    def GetSymm(self):
        unitcell = []
        a=''
        b=''
        c=''
        alpha=''
        beta=''
        gamma=''
        spg=''
        for line in open(self.pdb):
            if line.startswith('CRYST1'):
                a=line.split()[1]
                b=line.split()[2]
                c=line.split()[3]
                alpha=line.split()[4]
                beta=line.split()[5]
                gamma=line.split()[6]
#                spg=line[55:len(line)-1]
                spg=line[55:65]
        return [a,b,c,alpha,beta,gamma,spg]

    def get_spg_from_pdb(self):
        spg=self.GetSymm()[6]
        return spg

    def get_spg_number_from_pdb(self):
        spg=self.get_spg_from_pdb().replace(' ','').lower().replace('\n','').replace('\r','')
        spg_number='0'
        for key in self.translate_spg_to_number_dict:
            if key==spg:
                spg_number=str(self.translate_spg_to_number_dict[key])
                break
        return spg_number

    def get_bravais_lattice_from_pdb(self):
        bravais_lattice=''
        spg_number=self.get_spg_number_from_pdb()
        bravais_lattice=self.get_bravais_lattice_from_spg_number(spg_number)
        return bravais_lattice

    def get_bravais_lattice_from_spg_number(self,number):
        lattice=''
        for bravaislattice in self.space_group_dict:
            for spg_number in self.space_group_dict[bravaislattice]:
                if str(spg_number)==str(number):
                    lattice=bravaislattice
        return lattice

    def get_pointgroup_from_pdb(self):
        pointgroup=''
        spg_number=self.get_spg_number_from_pdb()
        pointgroup=self.get_point_group_from_spg_number(spg_number)
        return pointgroup

    def get_unit_cell_from_pdb(self):
        symm=self.GetSymm()
        unit_cell=symm[0:6]
        return unit_cell

    def calc_unitcell_volume_from_pdb(self):
        spg_number=self.get_spg_number_from_pdb()
        lattice=self.get_bravais_lattice_from_spg_number(spg_number)
        unitcell=self.get_unit_cell_from_pdb()
        a=float(unitcell[0])
        b=float(unitcell[1])
        c=float(unitcell[2])
        alpha=math.radians(float(unitcell[3]))
        beta= math.radians(float(unitcell[4]))
        gamma=math.radians(float(unitcell[5]))
        unitcell_volume=0
        if lattice=='triclinic':
            unitcell_volume=a*b*c* \
                            math.sqrt((1-math.cos(alpha)**2-math.cos(beta)**2-math.cos(gamma)**2) \
                                      +2*(math.cos(alpha)*math.cos(beta)*math.cos(gamma)))
#        if lattice=='monoclinic':
        if 'monoclinic' in lattice:
            unitcell_volume=round(a*b*c*math.sin(beta),1)
        if lattice=='orthorhombic' or lattice=='tetragonal' or lattice=='cubic':
            unitcell_volume=round(a*b*c,1)
        if lattice=='hexagonal' or lattice=='rhombohedral':
            unitcell_volume=round(a*b*c*(math.sin(math.radians(60))),1)
        return unitcell_volume


    def MatthewsCoefficient(self,sequence):
        chain=self.GetProteinChains()
        symm=self.GetSymm()
        nres=0
        for char in sequence:
            if char != ' ': nres+=1

        cmd = (
            '#!/bin/csh -f\n'
            'matthews_coef << eof\n'
            ' cell %s\n' %(symm[0]+' '+symm[1]+' '+symm[2]+' '+symm[3]+' '+symm[4]+' '+symm[5])+
            ' symm {0!s}\n'.format(symm[6].replace(' ',''))+
            ' nres {0!s}\n'.format(nres)+
            ' auto\n'
            ' end\n'
            'eof\n'
            )
        Log = subprocess.check_output(cmd, shell=True)
        MatthewsCoeff='n/a'
        Solvent='n/a'
        found=0
        for line in Log.split(os.linesep):
            if found:
                if int(line.split()[0]) == len(chain):
                    MatthewsCoeff=line.split()[1]
                    Solvent=line.split()[2]
                    break
            if line.startswith('_______________'):
                found=1

        return [MatthewsCoeff,Solvent]


    def find_ligands(self):
        Ligands = []
        # need to count residue numbers in case of alternative conformations
        ResiNum=[]
        for line in open(self.pdb):
            if (line.startswith('ATOM') or line.startswith('HETATM')) \
                    and line[17:20].replace(' ','') not in self.AminoAcids+self.Solvents+self.Ions:
                if [line[17:20],line[21:22],line[23:26]] not in Ligands:
                    Ligands.append([line[17:20],line[21:22],line[23:26]])
        return Ligands

    def save_ligands_to_pdb(self):
        Ligands=self.find_ligands()
        if not Ligands == []:
            for n,item in enumerate(Ligands):
                pdb=''
                for line in open(self.pdb):
                    if line.startswith('CRYST'):
                        pdb+=line
                    if (line.startswith('ATOM') or line.startswith('HETATM')) and line[17:20]==item[0] and line[21:22]==item[1] and line[23:26]==item[2]:
                        pdb=pdb+line
                f=open('ligand_{0!s}.pdb'.format(n),'w')
                f.write(pdb)
                f.close()
        return Ligands

    def save_ligands_to_pdb_to_directory(self,outDir):
        Ligands=self.find_ligands()
        if not Ligands == []:
            for n,item in enumerate(Ligands):
                pdb=''
                for line in open(self.pdb):
                    if line.startswith('CRYST'):
                        pdb+=line
                    if (line.startswith('ATOM') or line.startswith('HETATM')) and line[17:20]==item[0] and line[21:22]==item[1] and line[23:26]==item[2]:
                        pdb=pdb+line
                f=open(os.path.join(outDir,'ligand_{0!s}.pdb'.format(n)),'w')
                f.write(pdb)
                f.close()
        return Ligands


    def save_all_ligands_to_pdb(self,outDir):
        Ligands=self.find_ligands()
        if not Ligands == []:
            pdb=''
            for line in open(self.pdb):
                if line.startswith('CRYST'):
                    pdb+=line
            for n,item in enumerate(Ligands):
                for line in open(self.pdb):
                    if (line.startswith('ATOM') or line.startswith('HETATM')) and line[17:20]==item[0] and line[21:22]==item[1] and line[23:26]==item[2]:
                        pdb+=line
            f=open(os.path.join(outDir,'all_ligands.pdb'),'w')
            f.write(pdb)
            f.close()

    def save_specific_ligands_to_pdb(self,resname,chainID,resseq,altLoc):
        pdb=''
        outDir=self.pdb[:self.pdb.rfind('/')]
        for line in open(self.pdb):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resname_line=str(line[17:20]).replace(' ','')
                chainID_line=str(line[21:23]).replace(' ','')
                resseq_line=str(line[23:26]).replace(' ','')
                altLoc_line=str(line[16:17]).replace(' ','')
                if resname_line==str(resname) and chainID_line==str(chainID) and resseq_line==str(resseq) and altLoc_line==str(altLoc):
                    pdb=pdb+line

        if pdb != '':
            f=open('{0!s}/ligand_{1!s}_{2!s}_{3!s}_{4!s}.pdb'.format(outDir, str(resname), str(chainID), str(resseq), str(altLoc)),'w')
            f.write(pdb)
            f.close()

    def find_xce_ligand_details(self):
        Ligands = []
        for line in open(self.pdb):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resname=str(line[17:20]).replace(' ','')
                if resname in self.xce_ligands:
                    chainID=str(line[21:23]).replace(' ','')
                    resseq=str(line[23:26]).replace(' ','')
                    altLoc=str(line[16:17]).replace(' ','')
                    if [resname,chainID,resseq,altLoc] not in Ligands:
                        Ligands.append([resname,chainID,resseq,altLoc])
        return Ligands

    def ligand_details_as_list(self):
        Ligands = []
        for line in open(self.pdb):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resname=str(line[17:20]).replace(' ','')
                if resname in self.xce_ligands:
                    chainID=str(line[21:23]).replace(' ','')
                    resseq=str(line[23:26]).replace(' ','')
                    altLoc=str(line[16:17]).replace(' ','')
                    occupancy=str(line[56:60]).replace(' ','')
                    if [resname,chainID,resseq,altLoc,occupancy] not in Ligands:
                        Ligands.append([resname,chainID,resseq,altLoc,occupancy])
        return Ligands

    def residue_details_as_list(self,pdbin):
        residueList = []
        for line in open(pdbin):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                resname=str(line[17:20]).replace(' ','')
                chainID=str(line[21:23]).replace(' ','')
                resseq=str(line[23:26]).replace(' ','')
                altLoc=str(line[16:17]).replace(' ','')
                occupancy=str(line[56:60]).replace(' ','')
                if [resname,chainID,resseq,altLoc,occupancy] not in residueList:
                    residueList.append([resname,chainID,resseq,altLoc,occupancy])
        return residueList


    def check_occupancies(self):
        errorText=''
        warningText=''
        residueDict = {}
        for line in open(self.pdb):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atomname=str(line[12:16]).replace(' ','')
                resname=str(line[17:20]).replace(' ','')
                chainID=str(line[21:22]).replace(' ','')
                resseq=str(line[22:26]).replace(' ','')
                altLoc=str(line[16:17]).replace(' ','')
                occupancy=str(line[56:60]).replace(' ','')
                if resname+'-'+chainID+'-'+resseq not in residueDict:
                    residueDict[resname+'-'+chainID+'-'+resseq] = []
                if altLoc == '':
                    altLoc='0'
                residueDict[resname+'-'+chainID+'-'+resseq].append([altLoc,occupancy,atomname,resname,chainID,resseq])
        for item in residueDict:
            altLocDict={}
            for atom in residueDict[item]:
                altLoc=atom[0]
                occupancy=atom[1]
                atomname=atom[2]
                resname=atom[3]
                chainID=atom[4]
                resseq=atom[5]
                if altLoc not in altLocDict:
                    altLocDict[altLoc]=[]
                altLocDict[altLoc].append([occupancy,atomname,resname,chainID,resseq])
                if float(occupancy) > 1:
                    errorText+='ERROR: %s %s %s (%s) %s: occupancy is %s\n' %(chainID,resname,resseq,altLoc,atomname,occupancy)
            occupancySumList=[]
            for altLoc in altLocDict:
                occupancySum=0.0
                nAtom=float(len(altLocDict[altLoc]))
                for n,atom in enumerate(altLocDict[altLoc]):
                    occupancy=atom[0]
                    occupancySum+=float(occupancy)
                    atomname=atom[1]
                    resname=atom[2]
                    chainID=atom[3]
                    resseq=atom[4]
                    if n==0:
                        occStart=occupancy
                    else:
                        if occupancy != occStart:
                            warningText+='%s %s %s (%s) %s: occupancy differs for altLoc -> %s\n' %(chainID,resname,resseq,altLoc,atomname,occupancy)
                occupancySumList.append(occupancySum/nAtom)
            occAdd=0.0
            for entry in occupancySumList:
                occAdd+=entry
            if occAdd > 1:
                errorText+='ERROR: '+item+' -> summarised occupanies of alternative conformations are > 1.0 ('+str(occupancySumList)+')\n'
            
        return errorText,warningText

    def get_xyz_coordinated_of_residue(self,chain,number):
        X=0.0
        Y=0.0
        Z=0.0
        # pdb definition see: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        for line in open(self.pdb):
            if (line.startswith('ATOM') or line.startswith('HETATM')) and line[21:22]==chain and line[22:26].replace(' ','')==str(number):
                X=float(line[30:38])
                Y=float(line[38:46])
                Z=float(line[46:54])
                break
        return X,Y,Z

    def get_center_of_gravity_of_residue_ish(self,chain,number):
        X=0.0
        Y=0.0
        Z=0.0
        x_list=[]
        y_list=[]
        z_list=[]
        # pdb definition see: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        for line in open(self.pdb):
            if (line.startswith('ATOM') or line.startswith('HETATM')) and line[21:22]==chain and line[22:26].replace(' ','')==str(number):
                X=float(line[30:38])
                x_list.append(X)
                Y=float(line[38:46])
                y_list.append(Y)
                Z=float(line[46:54])
                z_list.append(Z)
        # 'ish' because it's not really the centre of gravity, but the the middle of the min/max of each x,y,z
        X=((max(x_list)-min(x_list))/2)+min(x_list)
        Y=((max(y_list)-min(y_list))/2)+min(y_list)
        Z=((max(z_list)-min(z_list))/2)+min(z_list)
        return X,Y,Z

    def get_center_of_gravity_of_molecule_ish(self):
        X=0.0
        Y=0.0
        Z=0.0
        x_list=[]
        y_list=[]
        z_list=[]
        # pdb definition see: http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        for line in open(self.pdb):
            if (line.startswith('ATOM') or line.startswith('HETATM')):
                X=float(line[30:38])
                x_list.append(X)
                Y=float(line[38:46])
                y_list.append(Y)
                Z=float(line[46:54])
                z_list.append(Z)
        # 'ish' because it's not really the centre of gravity, but the the middle of the min/max of each x,y,z
        X=((max(x_list)-min(x_list))/2)+min(x_list)
        Y=((max(y_list)-min(y_list))/2)+min(y_list)
        Z=((max(z_list)-min(z_list))/2)+min(z_list)
        return X,Y,Z

    def get_init_pdb_as_list(self):
        return self.get_pdb_as_list(self.pdb)

    def get_pdb_as_list(self,pdbin):
        pdb_list=[]
        for line in open(pdbin):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_line=str(line[12:16]).replace(' ','')
                resname_line=str(line[17:20]).replace(' ','')
                chainID_line=str(line[21:23]).replace(' ','')
                resseq_line=str(line[23:26]).replace(' ','')
                altLoc_line=str(line[16:17]).replace(' ','')
                pdb_list.append([atom_line,altLoc_line,resname_line,chainID_line,resseq_line])
        return pdb_list

    def ElementDict(self,resname,chainID,resseq,altLoc):

        ElementDict = {
            'C':    0,
            'N':    0,
            'O':    0,
            'P':    0,
            'S':    0,
            'BR':   0,
            'CL':   0,
            'I':    0,
            'F':    0   }

        for line in open(self.pdb):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_line=str(line[12:16]).replace(' ','')
                resname_line=str(line[17:20]).replace(' ','')
                chainID_line=str(line[21:23]).replace(' ','')
                resseq_line=str(line[23:26]).replace(' ','')
                altLoc_line=str(line[16:17]).replace(' ','')
                element_line=   str(line[66:78]).replace(' ','')
                if resname_line==resname and chainID_line==chainID and resseq_line==resseq and altLoc_line==altLoc:
                    if element_line.upper() in ElementDict:
                        ElementDict[element_line.upper()] += 1

        return ElementDict


    def update_residue(self,resname_old,chainID_old,resseq_old,altLoc_old,occupancy_old,resname_new,chainID_new,resseq_new,altLoc_new,occupancy_new):
        outPDB=''
        for line in open(self.pdb):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_line=      str(line[0:6])
                serial_line=    str(line[6:11])
                atomName_line=  str(line[11:16])
                altLoc_line=    str(line[16:17])
                resname_line=   str(line[17:20])
                chainID_line=   str(line[20:23])
                resseq_line=    str(line[23:26])
                insert_line=    str(line[26:30])
                x_line=         str(line[30:38])
                y_line=         str(line[38:46])
                z_line=         str(line[46:54])
                occupancy_line= str(line[54:60])
                Bfac_line=      str(line[60:66])
                element_line=   str(line[66:78])
                charge_line=    str(line[78:80])
                if resname_line.replace(' ','')==resname_old and chainID_line.replace(' ','')==chainID_old and resseq_line.replace(' ','')==resseq_old and altLoc_line.replace(' ','')==altLoc_old:

                    NewResname=resname_line.replace(resname_old,resname_new)

                    NewChain    = chainID_line.replace(chainID_old,chainID_new)

                    if len(resseq_old) < len(resseq_new):
                        NewResseq=resseq_line.replace(resseq_old,' '*(len(resseq_old)-len(resseq_new))+resseq_new)
                    elif len(resseq_old) == len(resseq_new):
                        NewResseq=resseq_line.replace(resseq_old,resseq_new)
                    else:
                        NewResseq=resseq_line[len(resseq_new)-len(resseq_old):].replace(resseq_old,resseq_new)

                    if altLoc_line == ' ':
                        NewAltloc=altLoc_line.replace(' ',altLoc_new)
                    else:
                        NewAltloc=altLoc_line.replace(altLoc_old,altLoc_new)

                    NewOccupany=occupancy_line.replace(occupancy_old,occupancy_new)

                    outPDB+=atom_line+serial_line+atomName_line+NewAltloc+NewResname+NewChain+NewResseq+insert_line+x_line+y_line+z_line+NewOccupany+Bfac_line+element_line+charge_line
                else:
                    outPDB+=line
            else:
                outPDB+=line

        f=open(self.pdb,'w')
        f.write(outPDB)
        f.close()

    def merge_pdb_file(self,pdbin):
        outPDB=''

        # only merge if residue with same resname, chain and resnumber not already in self.pdb!!!
        residueListReference=   self.residue_details_as_list(self.pdb)
        residueListPDBin=       self.residue_details_as_list(pdbin)
        DuplicateResidue=False
        for entry in residueListPDBin:
            if entry in residueListReference:
                print 'residue already exisits in '+self.pdb+':'
                print str(entry)
                DuplicateResidue=True

        if not DuplicateResidue:
            for line in open(self.pdb):
                if not line.startswith('END'):
                    outPDB+=line

            for line in open(pdbin):
                if not line.startswith('END'):
                    outPDB+=line

            outPDB+='END\n'

            f=open(self.pdb,'w')
            f.write(outPDB)
            f.close()
        else:
            print 'cannot merge pdb files!'


    def get_symmetry_operators(self):
        symop = []
        spg_number=str(self.get_spg_number_from_pdb())
#        print 'spg',spg_number
        foundSPG=False
        if os.path.isfile(os.path.join(os.getenv('CCP4'),'lib','data','symop.lib')):
#            print 'a',os.path.join(os.getenv('CCP4'),'lib','data','symop.lib')
            for line in open(os.path.join(os.getenv('CCP4'),'lib','data','symop.lib')):
#                print line.split()[0]
                if foundSPG:
                    if line.startswith(' '):
                        tmp=line.replace('\n','')
                        symop.append(tmp.split(','))
                    else:
                        break
                if line.split()[0] == spg_number:
                    foundSPG=True
        else:
            print "CCP4 environmental variable is not set"

        return symop


    def save_sym_equivalents(self):
        unit_cell=self.get_unit_cell_from_pdb()
        spg=self.get_spg_from_pdb()
        symop=self.get_symmetry_operators()
        outDir=self.pdb[:self.pdb.rfind('/')]
        root=self.pdb[self.pdb.rfind('/')+1:self.pdb.rfind('.')]
        outPDB=os.path.join(outDir,root+'_sym.pdb')
        pdbset = (  '#!'+os.getenv('SHELL')+'\n'
                    'pdbset xyzin %s xyzout %s << eof\n' %(self.pdb,outPDB)+
                    'cell {0!s}\n'.format((str(','.join(unit_cell))))+
                    'spacegroup {0!s}\n'.format(spg)  )
        for op in symop:
            pdbset+='SYMGEN '+','.join(op)+'\n'
        pdbset+='eof\n'
        os.system(pdbset)

    def save_sym_equivalents_of_ligands_in_pdb(self,pdbIN):
        unit_cell=self.get_unit_cell_from_pdb()
        spg=self.get_spg_from_pdb()
        symop=self.get_symmetry_operators()
        outDir=pdbIN[:pdbIN.rfind('/')]
        root=pdbIN[pdbIN.rfind('/')+1:pdbIN.rfind('.')]
        pdbset = (  '#!'+os.getenv('SHELL')+'\n'
                    'pdbset xyzin %s xyzout %s/%s_0.pdb << eof\n' %(pdbIN,outDir,root)+
                    'cell {0!s}\n'.format((str(','.join(unit_cell))))+
                    'spacegroup {0!s}\n'.format(spg)  )
        for op in symop:
            pdbset+='SYMGEN '+','.join(op)+'\n'
        pdbset+='eof\n'
        os.system(pdbset)
        return os.path.join(outDir,root+'_0.pdb')

    def save_sym_equivalents_of_ligands_in_pdb_as_one_file_per_ligand(self,pdbIN):
        nres=0
        for line in open(pdbIN):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                nres+=1
        print 'NRES',nres
        unit_cell=self.get_unit_cell_from_pdb()
        spg=self.get_spg_from_pdb()
        symop=self.get_symmetry_operators()
        outDir=pdbIN[:pdbIN.rfind('/')]
        root=pdbIN[pdbIN.rfind('/')+1:pdbIN.rfind('.')]
        pdbset = (  '#!'+os.getenv('SHELL')+'\n'
                    'pdbset xyzin %s xyzout %s/out.pdb << eof  2> /dev/null\n' %(pdbIN,outDir)+
                    'cell {0!s}\n'.format((str(','.join(unit_cell))))+
                    'spacegroup {0!s}\n'.format(spg)  )
        for op in symop:
            pdbset+='SYMGEN '+','.join(op)+'\n'
        pdbset+='eof\n'
        os.system(pdbset)

        pdbList=[]
        counter=0
        out=''
        Cryst=''
        n=0
        for line in open(os.path.join(outDir,'out.pdb')):
            if line.startswith('CRYST'):
                Cryst=line
            if line.startswith('ATOM') or line.startswith('HETATM'):
                if counter < nres-1:
                    if counter == 0:
                        out=Cryst
                    out+=line
                    counter+=1
                else:
                    out+=line
                    pdbList.append(os.path.join(outDir,root+'_'+str(n)+'.pdb'))
                    f=open(os.path.join(outDir,root+'_'+str(n)+'.pdb'),'w')
                    f.write(out)
                    f.close()
                    out=''
                    counter=0
                    n+=1

        return pdbList

    def save_surounding_unit_cells(self,pdbIN):
        translations = [

                    [   0,  1,  0   ],
                    [   0,  1,  1   ],
                    [   0,  0,  1   ],
                    [   0,  -1, 1   ],
                    [   0,  -1, 0   ],
                    [   0,  -1, -1  ],
                    [   0,  0,  -1  ],
                    [   0,  1,  -1  ],

                    [   1,  0,  0   ],
                    [   1,  1,  0   ],
                    [   1,  1,  1   ],
                    [   1,  0,  1   ],
                    [   1,  -1, 1   ],
                    [   1,  -1, 0   ],
                    [   1,  -1, -1  ],
                    [   1,  0,  -1  ],
                    [   1,  1,  -1  ],

                    [   -1, 0,  0   ],
                    [   -1, 1,  0   ],
                    [   -1, 1,  1   ],
                    [   -1, 0,  1   ],
                    [   -1, -1, 1   ],
                    [   -1, -1, 0   ],
                    [   -1, -1, -1  ],
                    [   -1, 0,  -1  ],
                    [   -1, 1,  -1  ],

            ]

        unit_cell=self.get_unit_cell_from_pdb()
        spg=self.get_spg_from_pdb()
        outDir=pdbIN[:pdbIN.rfind('/')]
        root=pdbIN[pdbIN.rfind('/')+1:pdbIN.rfind('.')]

        for n,shift in enumerate(translations):
            pdbset = (  '#!'+os.getenv('SHELL')+'\n'
                        'pdbset xyzin %s/%s.pdb xyzout %s/%s_%s.pdb << eof  2> /dev/null\n' %(outDir,root,outDir,root,str(n+1))+
                        'cell {0!s}\n'.format((str(','.join(unit_cell))))+
                        'spacegroup {0!s}\n'.format(spg)+
                        'shift fractional {0!s}\n'.format(str(shift).replace('[','').replace(']',''))+
                        'eof\n'            )
            print 'pdinin',pdbIN
            print 'outdir',outDir
            print 'root',root
            print pdbset
            os.system(pdbset)



class logtools:

    def __init__(self,logfile):
        self.logfile=logfile

    def phenix_molprobity(self):
        QualityIndicators = {   'MolprobityScore':              'n/a',
                                'MolprobityScoreColor':         'gray',
                                'RamachandranOutliers':         'n/a',
                                'RamachandranOutliersColor':    'gray',
                                'RamachandranFavored':          'n/a',
                                'RamachandranFavoredColor':     'gray'   }


        # Molprobity = validation_summary.txt
        if os.path.isfile(self.logfile):
            for line in open(self.logfile):
                if 'molprobity score' in line.lower():
                    if len(line.split()) >= 4:
                        QualityIndicators['MolprobityScore'] = line.split()[3]
                        try:
                            if float(line.split()[3]) < 2:                                 QualityIndicators['MolprobityScoreColor'] = 'green'
                            if float(line.split()[3]) >= 2 and float(line.split()[3]) < 3: QualityIndicators['MolprobityScoreColor'] = 'orange'
                            if float(line.split()[3]) >= 3:                                QualityIndicators['MolprobityScoreColor'] = 'red'
                        except ValueError:
                            pass
#                if line.lower().startswith('  ramachandran outliers ='):
                if 'ramachandran outliers' in line.lower():
                    if len(line.split()) >= 4:
                        QualityIndicators['RamachandranOutliers'] = line.split()[3]
                        try:
                            if float(line.split()[3]) < 0.3:                                 QualityIndicators['RamachandranOutliersColor'] = 'green'
                            if float(line.split()[3]) >= 0.3 and float(line.split()[3]) < 1: QualityIndicators['RamachandranOutliersColor'] = 'orange'
                            if float(line.split()[3]) >= 1:                                  QualityIndicators['RamachandranOutliersColor'] = 'red'
                        except ValueError:
                            pass
#                if line.startswith('               favored  ='):
                if 'favored' in line.lower():
                    if len(line.split()) == 4:
                        QualityIndicators['RamachandranFavored'] = line.split()[2]
                        try:
                            if float(line.split()[2]) < 90:                                  QualityIndicators['RamachandranFavoredColor'] = 'red'
                            if float(line.split()[2]) >= 90 and float(line.split()[2]) < 98: QualityIndicators['RamachandranFavoredColor'] = 'orange'
                            if float(line.split()[2]) >= 98:                                 QualityIndicators['RamachandranFavoredColor'] = 'green'
                        except ValueError:
                            pass

        return QualityIndicators


    def refmac_log(self):
        QualityIndicators = {   'MatrixWeight':              'n/a'   }

        # Matrix Weight
        if os.path.isfile(self.logfile):
            for line in open(self.logfile):
                if line.startswith(' Weight matrix') and len(line.split())==3:
                    QualityIndicators['MatrixWeight']=line.split()[2]

        return QualityIndicators



class reference:

    def __init__(self,sample_mtz,reference_file_list):
        self.sample_mtz=sample_mtz
        self.reference_file_list=reference_file_list
#        print reference_file_list

    def find_suitable_reference(self,allowed_unitcell_difference_percent):
        found_suitable_reference=False
        unitcell_reference='n/a'
        reference=''
        spg_reference='n/a'
        unitcell_difference='n/a'
        resolution_high='n/a'
        spg_autoproc='n/a'
        unitcell_autoproc='n/a'
        if os.path.isfile(self.sample_mtz):
            mtz_autoproc=mtztools(self.sample_mtz).get_all_values_as_dict()
            resolution_high=mtz_autoproc['resolution_high']
            spg_autoproc=mtz_autoproc['spacegroup']
            unitcell_autoproc=mtz_autoproc['unitcell']
            lattice_autoproc=mtz_autoproc['bravais_lattice']
            unitcell_volume_autoproc=mtz_autoproc['unitcell_volume']
            # check which reference file is most similar
            for o,reference_file in enumerate(self.reference_file_list):
#                print reference_file
                try:
                    if not reference_file[4]==0:
                        unitcell_difference=round((math.fabs(reference_file[4]-unitcell_volume_autoproc)/reference_file[4])*100,1)
                        # reference file is accepted when different in unitcell volume < 5%
                        # and both files have the same lattice type
                        if unitcell_difference < allowed_unitcell_difference_percent and lattice_autoproc==reference_file[3]:
                            spg_reference=reference_file[1]
                            unitcell_reference=reference_file[2]
                            reference=reference_file[0]
                            found_suitable_reference=True
                            break
                except IndexError:
                    pass
        return (spg_reference,unitcell_reference,reference,found_suitable_reference,
                resolution_high,spg_autoproc,unitcell_autoproc,unitcell_difference)


class misc:

    def calculate_distance_between_coordinates(self,x1,y1,z1,x2,y2,z2):
#        print '==> XCE: calculating distance between two coordinates'
        distance=0.0
        distance=math.sqrt(math.pow(float(x1)-float(x2),2)+math.pow(float(y1)-float(y2),2)+math.pow(float(z1)-float(z2),2))
        return distance

class smilestools(object):
    def __init__(self,smiles):
        self.smiles=smiles

    def ElementDict(self):
        ElementDict = {
            'C':    0,
            'N':    0,
            'O':    0,
            'P':    0,
            'S':    0,
            'BR':   0,
            'CL':   0,
            'I':    0,
            'F':    0   }

        m = Chem.MolFromSmiles(self.smiles)
        for atom in m.GetAtoms():
            if str(atom.GetSymbol()).upper() in ElementDict:
                ElementDict[str(atom.GetSymbol()).upper()] += 1

        return ElementDict

class maptools(object):
    def __init__(self,map):

        self.translate_spg_to_number_dict = {
            'p1': 1, 'p2': 3, 'p121': 4, 'c2': 5, 'c121': 5, 'p222': 16,
            'p2122': 17, 'p2212': 17, 'p2221': 17, 'p21212': 18, 'p21221': 18, 'p22121': 18,
            'p212121': 19, 'c2221': 20, 'c222': 21, 'f222': 22, 'i222': 23, 'i212121': 24,
            'p4': 75, 'p41': 76, 'p42': 77, 'p43': 78, 'i4': 79, 'i41': 80,
            'p422': 89, 'p4212': 90, 'p4122': 91, 'p41212': 92, 'p4222': 93, 'p42212': 94,
            'p4322': 95, 'p43212': 96, 'i422': 97, 'i4122': 98,
            'p3': 143, 'p31': 144, 'p32': 145, 'p312': 149, 'p321': 150, 'p3112': 151, 'p3121': 152,
            'p3212': 153, 'p3221': 154, 'p6': 168, 'p61': 169, 'p65': 170, 'p62': 171, 'p64': 172, 'p63': 173,
            'p622': 177, 'p6122': 178, 'p6522': 179, 'p6222': 180, 'p6422': 181, 'p6322': 182,
            'r3': 146, 'h3': 146, 'r32': 155, 'h32': 155,
            'p23': 195, 'f23': 196, 'i23': 197, 'p213': 198, 'i213': 199,
            'p432': 207, 'p4232': 208, 'f432': 209, 'f4132': 210, 'i432': 211, 'p4332': 212,
            'p4132': 213, 'i4132': 214        }

        self.map=map
        cmd = ( 'mapdump mapin {0!s} << eof\n'.format(self.map)+
                'end\n'
                'eof'   )
        mapdump=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        self.grid_sampling=[0,0,0]
        self.cell_dimensions=[]
        self.space_group_number=0
        for line in iter(mapdump.stdout.readline,''):
            if 'Grid sampling on x, y, z' in line:
                if len(line.split()) == 10:
                    self.grid_sampling=[line.split()[7],line.split()[8],line.split()[9]]
            if 'Cell dimensions' in line:
                if len(line.split()) == 9:
                    self.cell_dimensions=[line.split()[3],line.split()[4],line.split()[5],line.split()[6],line.split()[7],line.split()[8]]
            if 'Space-group' in line:
                if len(line.split()) == 3:
                    self.space_group_number=line.split()[2]

    def grid_sampling(self):
        return self.grid_sampling

    def cell_dimensions(self):
        return self.cell_dimensions

    def space_group_number(self):
        return self.space_group_number

    def space_group(self):
        space_group=''
        space_group_number=self.space_group_number
        for spg in self.translate_spg_to_number_dict:
            if str(self.translate_spg_to_number_dict[spg]) == str(space_group_number):
                space_group=str(spg)
        return space_group

