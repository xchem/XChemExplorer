import os
import glob
import sys


def merge_cifs(cpdDir):
    os.chdir(cpdDir)
    out = ""
    for cif in glob.glob("*.cif"):
        for line in open(cif):
            out += line
    f = open("merged.cif", "w")
    f.write(out)
    f.close()


if __name__ == "__main__":
    cpdDir = sys.argv[1]
    merge_cifs(cpdDir)
