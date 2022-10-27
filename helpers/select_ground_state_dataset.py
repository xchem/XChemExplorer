# last edited: 03/08/2017, 15:00

from XChemUtils import mtztools
from XChemUtils import parse
import os
import sys
import glob
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'lib'))


# - select datasets with highest resolution
# - select only those without an event map
# - take the one with the lowest Rfree


def find_highest_resolution_datasets(panddaDir):
    found = False
    datasetList = []
    for logFile in glob.glob(os.path.join(panddaDir, 'logs', '*.log')):
        for n, line in enumerate(open(logFile)):
            if line.startswith('Statistical Electron Density Characterisation') and len(line.split()) == 6:
                resolution = line.split()[5]
                found = True
                foundLine = n
            if found and n >= foundLine+3:
                if line.startswith('---'):
                    break
                else:
                    tmpLine = line.replace(' ', '').replace(
                        '\t', '').replace('\n', '').replace('\r', '')
                    for item in tmpLine.split(','):
                        if item != '':
                            datasetList.append(item)
    print datasetList
    return datasetList


def get_datasets_without_event_map(panddaDir, datasetList):
    datasetListwithoutEvent = []
    for dataset in datasetList:
        noEvent = True
        for files in glob.glob(os.path.join(panddaDir, 'processed_datasets', dataset, '*')):
            if 'event' in files:
                noEvent = False
                break
        if noEvent:
            datasetListwithoutEvent.append(dataset)
    print datasetListwithoutEvent
    return datasetListwithoutEvent


def select_dataset_with_lowest_Rfree(panddaDir, datasetListwithoutEvent):
    datasetList = []
    lowestRfree = ''
    for dataset in datasetListwithoutEvent:
        if os.path.isfile(os.path.join(panddaDir, 'processed_datasets', dataset, dataset+'-pandda-input.pdb')):
            stats = parse().PDBheader(os.path.join(
                panddaDir, 'processed_datasets', dataset, dataset+'-pandda-input.pdb'))
            Rfree = stats['Rfree']
            try:
                print dataset, Rfree, stats['ResolutionHigh']
                datasetList.append([dataset, float(Rfree)])
            except ValueError:
                pass
    if datasetList:
        lowestRfree = min(datasetList, key=lambda x: x[1])[0]
    return lowestRfree


def link_pdb_mtz_files(panddaDir, lowestRfree):
    targetDir = '/'.join(panddaDir.split('/')[:len(panddaDir.split('/'))-1])
    panddaFolder = panddaDir.split('/')[len(panddaDir.split('/'))-1]
    print targetDir
    print panddaFolder
    os.chdir(targetDir)
    if os.path.isfile(os.path.join(panddaDir, 'processed_datasets', lowestRfree, lowestRfree+'-pandda-input.pdb')):
        os.system('/bin/rm %s-ground-state.pdb 2> /dev/null' % lowestRfree)
        os.symlink(os.path.join(panddaFolder, 'processed_datasets', lowestRfree,
                                lowestRfree+'-pandda-input.pdb'), lowestRfree+'-ground-state.pdb')
    if os.path.isfile(os.path.join(panddaDir, 'processed_datasets', lowestRfree, lowestRfree+'-pandda-input.mtz')):
        os.system('/bin/rm %s-ground-state.free.mtz 2> /dev/null' %
                  lowestRfree)
        os.symlink(os.path.join(panddaFolder, 'processed_datasets', lowestRfree,
                                lowestRfree+'-pandda-input.mtz'), lowestRfree+'-ground-state.free.mtz')

    if os.path.isfile(os.path.join(panddaDir, 'processed_datasets', lowestRfree, lowestRfree+'-ground-state-mean-map.native.ccp4')):
        os.system(
            '/bin/rm %s-ground-state-mean-map.native.ccp4 2> /dev/null' % lowestRfree)
        os.symlink(os.path.join(panddaFolder, 'processed_datasets', lowestRfree, lowestRfree +
                                '-ground-state-mean-map.native.ccp4'), lowestRfree+'-ground-state-mean-map.native.ccp4')
    elif os.path.isfile(os.path.join(panddaDir, 'processed_datasets', lowestRfree, lowestRfree+'-ground-state-average-map.native.ccp4')):
        os.system(
            '/bin/rm %s-ground-state-mean-map.native.ccp4 2> /dev/null' % lowestRfree)
        os.symlink(os.path.join(panddaFolder, 'processed_datasets', lowestRfree, lowestRfree +
                                '-ground-state-average-map.native.ccp4'), lowestRfree+'-ground-state-mean-map.native.ccp4')

    convert_mean_map_to_mtz(
        lowestRfree+'-ground-state-mean-map.native.ccp4', lowestRfree+'-ground-state.free.mtz')


def convert_mean_map_to_mtz(emap, mtz):
    print 'converting ground-state-mean-map to MTZ'
    cmd = ('mapmask MAPIN %s MAPOUT %s << eof\n' % (emap, emap.replace('.ccp4', '.P1.ccp4')) +
           ' XYZLIM CELL\n'
           ' PAD 0.0\n'
           ' SYMMETRY 1\n'
           'eof\n')
    print cmd
    os.system(cmd)
    print '--->', mtz
    reso = mtztools(mtz).get_dmin()
    print '-> resolution:', reso
    cmd = ('module load phenix\n'
           'phenix.map_to_structure_factors %s d_min=%s\n' % (emap.replace('.ccp4', '.P1.ccp4'), reso) +
           '/bin/mv map_to_structure_factors.mtz %s' % emap.replace('.ccp4', '.mtz'))
    print cmd
    os.system(cmd)


if __name__ == '__main__':
    panddaDir = sys.argv[1]
    datasetList = find_highest_resolution_datasets(panddaDir)
    datasetListwithoutEvent = get_datasets_without_event_map(
        panddaDir, datasetList)
    lowestRfree = select_dataset_with_lowest_Rfree(
        panddaDir, datasetListwithoutEvent)
    link_pdb_mtz_files(panddaDir, lowestRfree)
