# last edited: 09/02/2017, 15:00

from __future__ import print_function
import os,glob
import sys
from datetime import datetime

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'),'lib'))

class startLog:

    def __init__(self,logfile):
        self.logfile=logfile

    def create_logfile(self,version):

        pasteVersion=version
        for i in range(0,20-len(version)):
            pasteVersion+=' '

        message = (
            '\n\n'
            '     #######################################################################\n'
            '     #                                                                     #\n'
            '     # XCHEMEXPLORER - multi dataset analysis                              #\n'
            '     #                                                                     #\n'
            '     # Version: %s                                       #\n' %pasteVersion+
            '     #                                                                     #\n'
            '     # Date: 01/03/2019                                                    #\n'
            '     #                                                                     #\n'
            '     # Authors: Tobias Krojer, Structural Genomics Consortium, Oxford, UK  #\n'
            '     #          tobias.krojer@sgc.ox.ac.uk                                 #\n'
            '     #                                                                     #\n'
            '     #          Rachael Skyner, Diamond Light Source, UK                   #\n'
            '     #          rachael.skyner@diamond.ac.uk                               #\n'
            '     #                                                                     #\n'
            '     #######################################################################\n'
            '\n'
        )

        if not os.path.isfile(self.logfile):
            os.system('touch '+self.logfile)
            message+='creating new logfile for the current XChemExplorer ('+version+') session:\n'+self.logfile+'\n'

        else:
            message+='writing into existing logfile for current XChemExplorer ('+version+') session:\n'+self.logfile+'\n'
        updateLog(self.logfile).insert(message)
        #print(message)

class depositLog:

    def __init__(self,logfile):
        self.logfile=logfile
        self.suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']

    def humansize(self,nbytes):
        if nbytes == 0: return '0 B'
        i = 0
        while nbytes >= 1024 and i < len(self.suffixes)-1:
            nbytes /= 1024.
            i += 1
        f = ('{0:.2f}'.format(nbytes)).rstrip('0').rstrip('.')
        return '{0!s} {1!s}'.format(f, self.suffixes[i])

    def modelInfo(self,xtal,structureType):
        message=(xtal+' @ '+structureType)
        updateLog(self.logfile).insert(message)

    def nEvents(self,xtal,n_eventMtz):
        updateLog(self.logfile).insert(xtal+' contains '+str(n_eventMtz)+ ' sites which are ready for deposition')

    def site_xml(self,xtal,xml):
        updateLog(self.logfile).insert('site description for '+xtal)
        updateLog(self.logfile).insert(xml)

    def text(self,message):
        updateLog(self.logfile).insert(message)

    def summary(self,n_toDeposit,success,failureDict,structureType,successDict):
        message=    (   '\n'
                        '==========================================================================\n'
                        'SUMMARY:\n'
                        '--------------------------------------------------------------------------\n'
                        ' - structure type:             %s\n'   %structureType+
                        ' - n(structures to deposit):   {0!s}\n'.format(str(n_toDeposit))+
                        ' - n(successful mmcif):        {0!s}\n'.format(str(success))+
                        '--------------------------------------------------------------------------\n'  )

        if successDict != {}:
            message+='\n => the following structures were successfully created:\n'
            message+='--------------------------------------------------------------------------\n'
        for key in sorted(successDict):
            message+='{0:20} {1:5} {2:25} {3:10} {4:25} {5:10}\n'.format(key,'->',successDict[key][0],str(self.humansize(successDict[key][1])),successDict[key][2],str(self.humansize(successDict[key][3])))


        if failureDict != {}:
            message+='--------------------------------------------------------------------------\n'
            message+='\n => the following structures had a problem:\n'
            message+='--------------------------------------------------------------------------\n'
        for key in failureDict:
            for entry in failureDict[key]:
                message+='{0!s} -> {1!s}\n'.format(key, entry)

        message+='==========================================================================\n\n'

        updateLog(self.logfile).insert(message)

class updateLog:

    def __init__(self,logfile):
        self.logfile=open(logfile, "a")

    def insert(self,message):
        present_time=datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')
        print( str(present_time)+' ==> XCE: '+message, file = self.logfile)
        print('==> XCE: '+message)

    def warning(self,message):
        present_time=datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')
        print( str(present_time)+' ==> XCE: WARNING! '+message, file = self.logfile)
        print('==> XCE: WARNING! '+message)

    def error(self,message):
        present_time=datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S')
        print( str(present_time)+' ==> XCE: ERROR!!! '+message, file = self.logfile)
        print('==> XCE: ERROR!!! '+message)


