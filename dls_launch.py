import os
import sys

if __name__ == "__main__":
    dir_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "dist", "xce-1.5.0-py2.7.egg"
    )
    sys.path.insert(0, dir_path)
    from xce.XChemExplorer import XChemExplorer

    XChemExplorer(sys.argv[1:])
