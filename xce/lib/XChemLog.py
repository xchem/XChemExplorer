import os
from datetime import datetime


class startLog:
    def __init__(self, logfile):
        self.logfile = logfile

    def create_logfile(self, version):
        pasteVersion = version
        for i in range(0, 20 - len(version)):
            pasteVersion += " "

        message = (
            "\n\n"
            "     ###################################################################\n"
            "     #                                                                 #\n"
            "     # XCHEMEXPLORER - multi dataset analysis                          #\n"
            "     #                                                                 #\n"
            "     # Version: %s                                   #\n"
            "     #                                                                 #\n"
            "     # Date: 01/08/2024                                               #\n"
            "     #                                                                 #\n"
            "     ###################################################################\n"
            "\n" % pasteVersion
        )

        if not os.path.isfile(self.logfile):
            os.system("touch " + self.logfile)
            message += (
                "creating new logfile for the current XChemExplorer ("
                + version
                + ") session:\n"
                + self.logfile
                + "\n"
            )

        else:
            message += (
                "writing into existing logfile for current XChemExplorer ("
                + version
                + ") session:\n"
                + self.logfile
                + "\n"
            )
        updateLog(self.logfile).insert(message)


class updateLog:
    def __init__(self, logfile):
        self.logfile = open(logfile, "a")

    def insert(self, message):
        present_time = datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S.%f")[:-4]
        self.logfile.write(str(present_time) + " ==> XCE: " + message)
        print("==> XCE: " + message)

    def warning(self, message):
        present_time = datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S.%f")[:-4]
        self.logfile.write(str(present_time) + " ==> XCE: WARNING! " + message)
        print("==> XCE: WARNING! " + message)

    def error(self, message):
        present_time = datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S.%f")[:-4]
        self.logfile.write(str(present_time) + " ==> XCE: ERROR!!! " + message)
        print("==> XCE: ERROR!!! " + message)

    def hint(self, message):
        present_time = datetime.strftime(datetime.now(), "%Y-%m-%d %H:%M:%S.%f")[:-4]
        self.logfile.write(str(present_time) + " ==> XCE: HINT -> " + message)
        print("==> XCE: HINT -> " + message)
