import sys

if __name__ =="__main__":
    sys.path.insert(0, "/dls_sw/i04-1/software/XChemExplorer/dist/xce-1.5.0-py2.7.egg")
    from xce.XChemExplorer import XChemExplorer
    XChemExplorer(sys.argv[1:])

