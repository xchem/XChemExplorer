# last edited: 05/07/2017, 15:00

import os,sys
import glob

sys.path.append(os.getenv('XChemExplorer_DIR')+'/lib')
import XChemLog
import XChemDB
import XChemMain

from XChemUtils import pdbtools
from XChemUtils import mtztools
import XChemUtils

class export_to_html:

    def __init__(self,htmlDir,projectDir,database,xce_logfile):
        self.htmlDir = htmlDir
        self.projectDir = projectDir
        self.Logfile=XChemLog.updateLog(xce_logfile)
        self.db=XChemDB.data_source(database)
        self.db_dict = None
        self.pdb = None
        self.protein_name = None

    def prepare(self,whichSamples):
        self.Logfile.insert('======== preparing HTML summary ========')
        self.makeFolders()
        self.copy_jscss()
        html = XChemMain.html_header()
        firstFile = True
        for xtal in self.db.samples_for_html_summary(whichSamples):
            self.db_dict = self.db.get_db_dict_for_sample(xtal)
            if firstFile:
                if self.db_dict['ProteinName'] == 'None':
                    self.Logfile.warning('could not determine protein name')
                    try:
                        self.protein_name = xtal.split('-')[0]
                        self.Logfile.warning('xtal name = %s => setting protein name to %s' %(xtal,self.protein_name))
                    except IndexError:
                        self.Logfile.warning('could not determine protein name from cystal name; setting to None')
                        self.protein_name = ''
                else:
                    self.protein_name = self.db_dict['ProteinName']
                    self.Logfile.insert('protein name is: ' + self.protein_name)
            self.copy_pdb(xtal)
            self.copy_mtz(xtal)
#            self.copy_electron_density(xtal)
            self.copy_ligand_files(xtal)
            os.chdir(os.path.join(self.projectDir,xtal))
            ligandDict = XChemUtils.pdbtools_gemmi('refine.pdb').center_of_mass_ligand_dict('LIG')
            for ligand in ligandDict:
#            for ligand in self.ligands_in_pdbFile(xtal):
                ligName = ligand.split('-')[0]
                ligChain = ligand.split('-')[1]
                ligNumber = ligand.split('-')[2]
#                eventMap = self.find_matching_event_map_from_database(xtal, ligand)
#                if eventMap:
#                    self.Logfile.insert('%s: using the following event map -> %s' %(xtal,eventMap))
#                    self.cut_and_copy_map(xtal, ligand+'.pdb', eventMap, xtal + '_' + ligand + '_event.ccp4','F','PHIF')
#                    eventMap = xtal + '_' + ligand + '_event.ccp4'
#                else:
#                    self.Logfile.error('%s: value of event map -> %s' %(xtal,eventMap))
                eventMap = ''
                if os.path.isfile(xtal + '_' + ligand + '_event.ccp4'):
                    eventMap = xtal + '_' + ligand + '_event.ccp4'
#                x,y,z = self.pdb.get_centre_of_gravity_of_residue(ligand)
                x = ligandDict[ligand][0]
                y = ligandDict[ligand][1]
                z = ligandDict[ligand][2]
                self.copy_spider_plot(xtal,ligand)
                pdbID = self.db_dict['Deposition_PDB_ID']
                compoundImage = xtal + '_' + self.db_dict['CompoundCode'] + '.png'
                compoundCIF = xtal + '_' + self.db_dict['CompoundCode'] + '.cif'
                residuePlot = xtal + '_' + ligand + '.png'
                pdb = xtal + '.pdb'
                event = xtal + '_' + ligand + '_event.ccp4'
                thumbNail = xtal + '_' + ligand + '_thumb.png'
                resoHigh = self.db_dict['DataProcessingResolutionHigh']
                spg = self.db_dict['RefinementSpaceGroup']
                unitCell = self.db_dict['DataProcessingUnitCell']
                t = ''
                for ax in unitCell.split():
                    t += str(round(float(ax),1)) + ' '
                unitCell = t[:-1]
                os.chdir(os.path.join(self.projectDir,xtal))
#                FWT = xtal + '-' + ligand + '_2fofc.ccp4'
#                self.cut_and_copy_map(xtal, ligand + '.pdb', '2fofc.map', FWT,'FWT','PHWT')
#                DELFWT = xtal + '-' + ligand + '_fofc.ccp4'
#                self.cut_and_copy_map(xtal, ligand + '.pdb', 'fofc.map', DELFWT,'DELFWT','PHDELWT')
                if os.path.isfile('refine.mtz'):
                    self.Logfile.insert('%s: found refine.mtz' %xtal)
                    FWT, PHWT, DELFWT, PHDELWT = XChemUtils.mtztools_gemmi('refine.mtz').get_map_labels()
                    XChemUtils.maptools().calculate_map('refine.mtz',FWT,PHWT)
                    XChemUtils.maptools().cut_map_around_ligand('refine.ccp4',ligand+'.pdb','7')
                    os.system('/bin/mv %s %s_%s_2fofc.ccp4' %('refine_mapmask.ccp4',xtal,ligand))
                    FWTmap = xtal + '-' + ligand + '_2fofc.ccp4'
                    XChemUtils.maptools().calculate_map('refine.mtz',DELFWT,PHDELWT)
                    XChemUtils.maptools().cut_map_around_ligand('refine.ccp4',ligand+'.pdb','7')
                    os.system('/bin/mv %s %s_%s_fofc.ccp4' %('refine_mapmask.ccp4',xtal,ligand))
                    DELFWTmap = xtal + '-' + ligand + '_fofc.ccp4'
                    self.copy_electron_density(xtal,ligand,eventMap)
                ligConfidence = self.db.get_ligand_confidence_for_ligand(xtal, ligChain, ligNumber, ligName)
                if ligConfidence.startswith('0'):
                    self.Logfile.warning('%s: ligand confidence of %s-%s-%s is %s; ignoring...' %(xtal,ligChain,ligNumber,ligName,ligConfidence))
                    self.Logfile.warning('%s: this seems unlikely because this structure is apparently ready for deposition' %xtal)
                    self.Logfile.warning('%s: will set it to "not assigned" for now, but please update in soakDB' %xtal)
                    ligConfidence = 'not assigned'
                modelStatus = self.db_dict['RefinementOutcome']
                if firstFile:
                    html += XChemMain.html_ngl(pdb,eventMap.replace(self.projectDir,''),FWTmap,DELFWTmap,ligand)
                    html += XChemMain.html_download(self.protein_name)
                    html += XChemMain.html_guide()
                    html += XChemMain.html_table_header()
                    firstFile = False

                html += XChemMain.html_table_row(xtal,pdbID,ligand,compoundImage,residuePlot,pdb,event,
                                                 thumbNail,resoHigh,spg,unitCell,FWTmap,DELFWTmap,ligConfidence,modelStatus)
                self.make_thumbnail(xtal,x,y,z,ligand,eventMap)
                self.prepare_for_download(xtal, pdb, event, compoundCIF, ligand)
        self.prepare_zip_archives()
#        html = XChemMain.html_download_all_section(html,self.protein_name)
        self.write_html_file(html)
        self.Logfile.insert('======== finished preparing HTML summary ========')


    def prepare_zip_archives(self):
        os.chdir(os.path.join(self.htmlDir,'files'))
        self.Logfile.insert('%s: preparing ZIP archive of all PDB files' %self.protein_name)
        os.system('zip %s_allPDBs.zip *.pdb' %self.protein_name)
        self.Logfile.insert('%s: preparing ZIP archive of all PanDDA event maps' %self.protein_name)
        os.system('zip %s_allEVENTmaps.zip *event.ccp4' %self.protein_name)
        self.Logfile.insert('%s: preparing ZIP archive of all CIF files' %self.protein_name)
        os.system('zip %s_allCIFs.zip *.cif' %self.protein_name)
        self.Logfile.insert('%s: preparing ZIP archive of all MTZ files' %self.protein_name)
        os.system('zip %s_allMTZs.zip *.mtz' %self.protein_name)


    def prepare_for_download(self,xtal,pdb,event,compoundCIF,ligID):
        os.chdir(os.path.join(self.htmlDir,'tmp'))
        self.Logfile.insert('%s: preparing files for download' %xtal)
        zip_in = ''

        if os.path.isfile('../files/%s' %pdb):
            os.system('/bin/cp ../files/%s .' %pdb)
            zip_in += pdb + ' '
        else:
            self.Logfile.error('%s: cannot find %s' %(xtal,pdb))

        if os.path.isfile('../files/%s' %event):
            os.system('/bin/cp ../files/%s .' %event)
            zip_in += event + ' '
        else:
            self.Logfile.error('%s: cannot find %s' %(xtal,event))

        if os.path.isfile('../files/%s' %compoundCIF):
            os.system('/bin/cp ../files/%s .' %compoundCIF)
            zip_in += compoundCIF + ' '
        else:
            self.Logfile.error('%s: cannot find %s' %(xtal,compoundCIF))

        if zip_in != '':
            self.Logfile.insert('%s: preparing zip file -> zip %s_%s.zip %s' %(xtal,xtal,ligID,zip_in))
            os.system('zip %s_%s.zip %s' %(xtal,ligID,zip_in))
            os.system('/bin/mv %s_%s.zip ../download' %(xtal,ligID))
            os.system('/bin/rm -f *')
        else:
            self.Logfile.error('%s: cannot find any input files for creating of zip archive of %s_%s' %(xtal,xtal,ligID))



    def copy_jscss(self):
        os.chdir(self.htmlDir)
        self.Logfile.insert('copying css and js files to ' + self.htmlDir)
        os.system('/bin/cp -r %s/web/jscss/css .' %os.getenv('XChemExplorer_DIR'))
        os.system('/bin/cp -r %s/web/jscss/js .' %os.getenv('XChemExplorer_DIR'))


    def make_thumbnail(self,xtal,x,y,z,ligID,eventMap):
        self.Logfile.insert('%s: making thumbnail of for %s and %s' %(xtal,ligID,eventMap))
        sampleDir = os.path.join(self.projectDir,xtal)
        os.chdir(sampleDir)
        if not os.path.isfile('%s_%s_thumb.png' %(xtal,ligID)):
            self.Logfile.insert('%s: preparing thumbnail image of %s' %(xtal,ligID))
            XChemMain.coot_prepare_input(x, y, z, ligID, sampleDir, eventMap)
            XChemMain.coot_write_raster_file(ligID,sampleDir)
            XChemMain.render_scene(xtal,ligID,sampleDir)
            XChemMain.make_thumbnail(xtal, ligID, sampleDir)
        if os.path.isfile('%s_%s_thumb.png' %(xtal,ligID)):
            self.Logfile.insert('%s: managed to prepare %s_%s_thumb.png' %(xtal,xtal,ligID))
            self.copy_thumbnail(xtal,sampleDir,ligID)
        else:
            self.Logfile.error('%s: could not generate %s_%s_thumb.png' %(xtal,xtal,ligID))

    def copy_thumbnail(self,xtal,sampleDir,ligID):
        os.chdir(os.path.join(self.htmlDir, 'png'))
        self.Logfile.insert('%s: copying %s_%s_thumb.png to html png' %(xtal,xtal,ligID))
        os.system('/bin/cp %s/%s_%s_thumb.png .' %(sampleDir,xtal,ligID))


    def makeFolders(self):
        self.Logfile.insert('preparing folders in html directory')
        os.chdir(self.htmlDir)
#        if not os.path.isdir('js'):
#            os.mkdir('js')
        if not os.path.isdir('tmp'):
            os.mkdir('tmp')
        if not os.path.isdir('png'):
            os.mkdir('png')
        if not os.path.isdir('files'):
            os.mkdir('files')
        if not os.path.isdir('download'):
            os.mkdir('download')

    def copy_pdb(self,xtal):
        os.chdir(os.path.join(self.htmlDir, 'files'))
        self.pdb = None
        if os.path.isfile(os.path.join(self.projectDir,xtal,'refine.split.bound-state.pdb')):
            self.pdb = pdbtools(os.path.join(self.projectDir,xtal,'refine.split.bound-state.pdb'))
            self.Logfile.insert('%s: copying refine.split.bound-state.pdb to html directory' %xtal)
            os.system('/bin/cp %s/refine.split.bound-state.pdb %s.pdb' %(os.path.join(self.projectDir,xtal),xtal))
        else:
            self.Logfile.error('%s: cannot find refine.split.bound-state.pdb' %xtal)

    def copy_mtz(self,xtal):
        os.chdir(os.path.join(self.htmlDir, 'files'))
        if os.path.isfile(os.path.join(self.projectDir,xtal,xtal+'.free.mtz')):
            self.Logfile.insert('%s: copying %s.free.mtz to html directory' %(xtal,xtal))
            os.system('/bin/cp %s/%s.free.mtz .' %(os.path.join(self.projectDir,xtal),xtal))
        else:
            self.Logfile.error('%s: cannot find %s.free.mtz' %(xtal,xtal))

    def copy_electron_density(self,xtal,ligand,eventMap):
        os.chdir(os.path.join(self.htmlDir, 'files'))
        if os.path.isfile(os.path.join(self.projectDir,xtal,xtal+'-'+ligand+'_2fofc.ccp4')):
            emap = xtal+'-'+ligand+'_2fofc.ccp4'
            self.Logfile.insert('%s: copying %s to html directory' %(xtal,emap))
            os.system('/bin/cp %s/%s %s' %(os.path.join(self.projectDir,xtal),emap,emap))
        else:
            self.Logfile.error('%s: cannot find %s' %(xtal,emap))

        if os.path.isfile(os.path.join(self.projectDir,xtal,xtal+'-'+ligand+'_fofc.ccp4')):
            emap = xtal+'-'+ligand+'_fofc.ccp4'
            self.Logfile.insert('%s: copying %s to html directory' %(xtal,emap))
            os.system('/bin/cp %s/%s %s' %(os.path.join(self.projectDir,xtal),emap,emap))
        else:
            self.Logfile.error('%s: cannot find %s' %(xtal,emap))

        if os.path.isfile(os.path.join(self.projectDir,xtal,eventMap)):
            self.Logfile.insert('%s: copying %s to html directory' %(xtal,eventMap))
            os.system('/bin/cp %s/%s %s' %(os.path.join(self.projectDir,xtal),eventMap,eventMap))
        else:
            self.Logfile.error('%s: cannot find %s' %(xtal,eventMap))

    def copy_ligand_files(self,xtal):
        os.chdir(os.path.join(self.htmlDir,'files'))

        if os.path.isfile(os.path.join(self.projectDir,xtal,self.db_dict['CompoundCode']+'.cif')):
            self.Logfile.insert('%s: copying compound cif file' %xtal)
            os.system('/bin/cp %s %s_%s' %(os.path.join(self.projectDir,xtal,self.db_dict['CompoundCode']+'.cif'),xtal,self.db_dict['CompoundCode']+'.cif'))
        else:
            self.Logfile.error('%s: cannot find compound cif file -> %s' %(xtal,self.db_dict['CompoundCode']+'.cif'))

        os.chdir(os.path.join(self.htmlDir,'png'))

        if os.path.isfile(os.path.join(self.projectDir,xtal,self.db_dict['CompoundCode']+'.png')):
            self.Logfile.insert('%s: copying compound png file' %xtal)
            os.system('/bin/cp %s %s_%s' %(os.path.join(self.projectDir,xtal,self.db_dict['CompoundCode']+'.png'),xtal,self.db_dict['CompoundCode']+'.png'))
        else:
            self.Logfile.error('%s: cannot find compound png file -> %s' %(xtal,self.db_dict['CompoundCode']+'.png'))

    def copy_spider_plot(self,xtal,ligID):
        os.chdir(os.path.join(self.htmlDir, 'png'))
        refPDB = os.path.join(self.projectDir,xtal,'refine.pdb')
        self.Logfile.insert('%s: looking for spider plots...' %xtal)
        if os.path.isfile(refPDB):
            refPDBreal = os.path.realpath(refPDB)[:os.path.realpath(refPDB).rfind('/')]
            plot = os.path.join(refPDBreal,'residue_plots',ligID.replace('LIG-','')+'.png')
            self.Logfile.insert(xtal + ': looking for ' + plot)
            if os.path.isfile(plot):
                self.Logfile.insert('%s: found %s' % (xtal, plot))
                self.Logfile.insert('%s: copying spider plot for %s' % (xtal, ligID.replace('LIG-', '') + '.png'))
                os.system('/bin/cp %s %s_%s.png' % (plot, xtal, ligID))
            else:
                self.Logfile.error('%s: cannot find spider plot for %s' %(xtal,ligID.replace('LIG-','')+'.png'))
                self.Logfile.warning('%s: using %s instead' %(xtal,'NO_SPIDER_PLOT_AVAILABLE.png'))
                os.system('/bin/cp %s %s_%s.png' % (os.path.join(os.getenv('XChemExplorer_DIR'), 'image', 'NO_SPIDER_PLOT_AVAILABLE.png'), xtal, ligID))
                #
        else:
            self.Logfile.error('%s: cannot find refine.pdb, i.e. cannot start looking for spider plots...' %xtal)


#    def ligands_in_pdbFile(self,xtal):
#        os.chdir(os.path.join(self.projectDir,xtal))
#        ligPDB = []
#        ligList = []
#        self.Logfile.insert('%s: reading ligands to type LIG in refine.split.bound-state.pdb' %xtal)
#        if os.path.isfile('refine.split.bound-state.pdb'):
#            ligPDB = self.pdb.save_residues_with_resname(os.path.join(self.projectDir, xtal), 'LIG')
#        else:
#            self.Logfile.error('%s: cannot find refine.split.bound-state.pdb' %xtal)
#        if not ligPDB:
#            self.Logfile.error('%s; cannot find any ligands of type LIG in refine.split.bound-state.pdb' %xtal)
#        else:
#            for lig in sorted(ligPDB):
#                ligList.append(lig.replace('.pdb', ''))
#        return ligList


#    def find_matching_event_map_from_database(self,xtal,ligID):
#        ligName = ligID.split('-')[0]
#        ligChain = ligID.split('-')[1]
#        ligNumber = ligID.split('-')[2]
#        eventMAP = self.db.get_event_map_for_ligand(xtal, ligChain, ligNumber, ligName)
#        self.Logfile.insert('%s: the database thinks the following event map belongs to %s: %s' %(xtal,ligID,eventMAP))
##        print 'event map', eventMAP
#        if eventMAP == '' or 'none' in str(eventMAP).lower():
#            self.Logfile.warning('%s: the respective field in the DB is apparently emtpy' %xtal)
#            self.Logfile.warning('%s: will try to determine ligand - event map relationship by checking CC...' %xtal)
#            eventMAP = self.find_matching_event_map(xtal,ligID)
#        elif not os.path.isfile(eventMAP):
#            self.Logfile.warning('%s: event map file does not exist!' %xtal)
#            self.Logfile.warning('%s: will try to determine ligand - event map relationship by checking CC...' %xtal)
#            eventMAP = self.find_matching_event_map(xtal,ligID)
#        else:
#            self.Logfile.insert('%s: found matching event map!' %xtal)
#        return eventMAP

#    def find_matching_event_map(self,xtal,ligID):
#        os.chdir(os.path.join(self.projectDir, xtal))
#        eventMAP = []
#        self.Logfile.insert('%s: trying to find fitting event maps for modelled ligands' %xtal)
#        if os.path.isfile('no_pandda_analysis_performed'):
#            self.Logfile.warning('%s: no pandda analysis performed; skipping this step...' %xtal)
#            return
#        ligCC = []
#        for mtz in sorted(glob.glob('*event*.native*P1.mtz')):
#            self.get_lig_cc(xtal, mtz, ligID+'.pdb')
#            cc = self.check_lig_cc(mtz.replace('.mtz', '_CC'+ligID+'.log'))
#            self.Logfile.insert('%s: %s -> CC = %s for %s' %(xtal,ligID,cc,mtz))
#            try:
#                ligCC.append([mtz,float(cc)])
#            except ValueError:
#                ligCC.append([mtz, 0.00])
#        try:
#            highestCC = max(ligCC, key=lambda x: x[0])[1]
#        except ValueError:
#            self.Logfile.error('%s: event maps are not yet converted to mtz files...' %xtal)
#            return
#        if highestCC == 0.00 or ligCC is []:
#            self.Logfile.error('%s: best CC of ligand %s for any event map is 0!' %(xtal,ligID))
#        else:
#            self.Logfile.insert('%s: selected event map -> CC(%s) = %s for %s' %(xtal,ligID,highestCC,mtz[mtz.rfind('/')+1:]))
#            eventMAP = mtz[mtz.rfind('/')+1:].replace('.P1.mtz','.ccp4')
#            if not os.path.isfile(eventMAP):
#                eventMAP = []
#            else:
#                self.cut_eventMAP(xtal,ligID,eventMAP)
#        return eventMAP

#    def cut_and_copy_map(self,xtal,pdbCentre,mapin,mapout,F,PHI):
#        os.chdir(os.path.join(self.projectDir, xtal))
#        self.Logfile.insert('%s: cutting density of %s around %s' %(xtal,mapin.replace('.ccp4','.P1.ccp4'),pdbCentre))
#        if os.path.isfile(mapout):
#            self.Logfile.warning('%s: removing map -> %s' %(xtal,mapout))
#            os.system('/bin/rm '+mapout)
##        else:
#
#        if mapin.endswith('.map') or mapin.endswith('.ccp4'):
#            cmd = (
#                'mapmask mapin %s mapout %s xyzin %s << eof\n'  %(mapin.replace('.ccp4','.P1.ccp4'),mapout,pdbCentre) +
#                ' border 12\n'
#                ' end\n'
#                'eof'
#            )

#            cmd = (
#                'cmapcut -mapin %s -pdbin %s -mapout %s' %(mtzin,pdbCentre,mapout)
#            )

#            cmd = (
#                "phenix.cut_out_density %s %s map_coeff_labels='%s,%s' cutout_model_radius=6 cutout_map_file_name=%s cutout_as_map=True" %(pdbCentre,mtzin,F,PHI,mapout)
#            )
#            self.Logfile.insert(xtal+': running command:\n'+cmd)
#            os.system(cmd)
#            if os.path.isfile(mapout):
#                self.Logfile.insert(xtal+': reduced event map successfully created')
#                self.Logfile.insert('%s: copying %s to %s/files' % (xtal, mapout, self.htmlDir))
#                os.system('/bin/cp %s %s/files' % (mapout, self.htmlDir))
#            else:
#                self.Logfile.error(xtal+': creation of event map failed')



#    def cut_eventMAP(self,xtal,ligID,eventMAP):
#        os.chdir(os.path.join(self.projectDir, xtal))
#        self.Logfile.insert('%s: cutting event map around ligand %s' %(xtal,ligID))
#        ligMAP = xtal + '_' + ligID + '.ccp4'
#        cmd = (
#            'mapmask mapin %s mapout %s xyzin %s << eof\n'  %(eventMAP,ligMAP,ligID+'.pdb') +
#            ' border 10\n'
#            ' end\n'
#            'eof'
#        )
#        os.system(cmd)
#        self.copy_eventMap(xtal, ligID, eventMAP)

#    def copy_eventMap(self,xtal,ligID,eventMAP):
#        os.chdir(os.path.join(self.htmlDir,'files'))
#        self.Logfile.insert('%s: copying event map for %s' %(xtal,ligID))
#        os.system('/bin/mv %s/%s_%s.ccp4 .' %(os.path.join(self.projectDir,xtal),xtal,ligID))

#    def get_lig_cc(self, xtal, mtz, lig):
#        ligID = lig.replace('.pdb','')
#        ccLog = mtz.replace('.mtz', '_CC'+ligID+'.log')
#        self.Logfile.insert('%s: calculating CC for %s in %s' %(xtal,lig,mtz))
##        if os.path.isfile(mtz.replace('.mtz', '_CC'+ligID+'.log')):
#        if os.path.isfile(ccLog) and os.path.getsize(ccLog) != 0:
#            self.Logfile.warning('logfile of CC analysis exists; skipping...')
#            return
##        cmd = ( 'module load phenix\n'
##                'phenix.get_cc_mtz_pdb %s %s > %s' % (mtz, lig, mtz.replace('.mtz', '_CC'+ligID+'.log')) )
#        os.system('/bin/rm %s' %mtz.replace('.mtz', '_CC'+ligID+'.log') )
#        cmd = ( 'phenix.get_cc_mtz_pdb %s %s > %s' % (mtz, lig, mtz.replace('.mtz', '_CC'+ligID+'.log')) )
#        os.system(cmd)

#    def check_lig_cc(self,log):
#        cc = 'n/a'
#        if os.path.isfile(log):
#            for line in open(log):
#                if line.startswith('local'):
#                    cc = line.split()[len(line.split()) - 1]
#        else:
#            self.Logfile.error('logfile does not exist: %s' %log)
#        return cc

    def write_html_file(self,html):
        os.chdir(self.htmlDir)
        self.Logfile.insert('writing index.html')
        html += XChemMain.html_footer()
        if os.path.isfile('index.html'):
            os.system('/bin/rm -f index.html')
        f = open('index.html','w')
        f.write(html)
        f.close()

