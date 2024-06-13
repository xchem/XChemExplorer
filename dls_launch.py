import os
import sys

if __name__ == "__main__":
    sys.path.insert(
        0, os.path.join(os.environ["XChemExplorer_DIR"], "dist", "xce-1.5.0-py2.7.egg")
    )
    from xce.XChemExplorer import XChemExplorer

    XChemExplorer(sys.argv[1:])
