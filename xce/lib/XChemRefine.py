import getpass
import glob
import os
from datetime import datetime

import gtk
import pygtk

from xce.lib import XChemLog
from xce.lib import XChemUtils
from xce.lib.cluster import sge, slurm

pygtk.require("2.0")


def GetSerial(ProjectPath, xtalID):
    # check if there were already previous refinements
    # if no: create a folder Refine_1
    # if yes: create a folder Refine_<max+1>
    temp = []
    found = 0
    if os.path.isdir(os.path.join(ProjectPath, xtalID)):
        for item in glob.glob(os.path.join(ProjectPath, xtalID, "*")):
            if item.startswith(os.path.join(ProjectPath, xtalID, "Refine_")):
                if item.endswith("-report"):
                    continue
                try:
                    temp.append(int(item[item.rfind("_") + 1 :]))
                except ValueError:
                    continue
                found = 1
    if found:
        Serial = max(temp) + 1
    else:
        Serial = 1
    return Serial


class RefineParams(object):
    def __init__(self, ProjectPath, xtalID, compoundID, datasource):
        self.ProjectPath = ProjectPath
        self.xtalID = xtalID
        self.compoundID = compoundID
        self.datasource = datasource

    def RefmacRefinementParams(self, RefmacParams):
        self.RefmacParams = RefmacParams
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", gtk.main_quit)
        self.window.set_border_width(10)
        self.window.set_title("Refmac Parameters")
        self.vbox = gtk.VBox()

        self.hbox1 = gtk.HBox()
        self.hbox1.add(gtk.Label("Refine"))
        self.cb = gtk.combo_box_new_text()
        self.cb.connect("changed", self.ChooseBfacRefinement)
        for item in ["isotropic", "anisotropic"]:
            self.cb.append_text(item)
        if "ISOT" in self.RefmacParams["BREF"]:
            self.cb.set_active(0)
        if "ANIS" in self.RefmacParams["BREF"]:
            self.cb.set_active(1)
        self.hbox1.add(self.cb)
        self.hbox1.add(gtk.Label("temperature factors"))
        self.vbox.add(self.hbox1)

        self.hbox2 = gtk.HBox()
        self.hbox2.add(gtk.Label("Number of Cycles: "))
        self.Ncycles = gtk.Entry()
        self.Ncycles.add_events(gtk.gdk.KEY_RELEASE_MASK)
        self.Ncycles.connect("key-release-event", self.on_key_release_Ncycles)
        self.Ncycles.set_text(self.RefmacParams["NCYCLES"])
        self.hbox2.add(self.Ncycles)
        self.vbox.add(self.hbox2)

        self.hbox3 = gtk.HBox()
        self.hbox3.add(gtk.Label("MATRIX WEIGHT: "))
        self.MATRIX_WEIGHT = gtk.Entry()
        self.MATRIX_WEIGHT.add_events(gtk.gdk.KEY_RELEASE_MASK)
        self.MATRIX_WEIGHT.connect(
            "key-release-event", self.on_key_release_MATRIX_WEIGHT
        )
        self.MATRIX_WEIGHT.set_text(self.RefmacParams["MATRIX_WEIGHT"])
        self.hbox3.add(self.MATRIX_WEIGHT)
        self.vbox.add(self.hbox3)

        self.TLS = gtk.CheckButton("TLS (find TLS groups with phenix.find_tls_groups)")
        self.TLS.connect("toggled", self.TLSCallback)
        if self.RefmacParams["TLS"] == "refi tlsc 10\n":
            self.TLS.set_active(True)
        self.vbox.pack_start(self.TLS, False)

        self.NCS = gtk.CheckButton("NCS (if applicable")
        self.NCS.connect("toggled", self.NCSCallback)
        if self.RefmacParams["NCS"] == "NCSR LOCAL\n":
            self.NCS.set_active(True)
        self.vbox.pack_start(self.NCS, False)

        self.TWIN = gtk.CheckButton("Twin?")
        self.TWIN.connect("toggled", self.TWINCallback)
        if self.RefmacParams["TWIN"] == "TWIN\n":
            self.TWIN.set_active(True)
        self.vbox.pack_start(self.TWIN, False)

        self.WATER = gtk.CheckButton("Update waters (BUSTER)")
        self.WATER.connect("toggled", self.WATERCallback)
        if self.RefmacParams["WATER"] == "WATER\n":
            self.WATER.set_active(True)
        self.vbox.pack_start(self.WATER, False)

        self.LIGOCC = gtk.CheckButton("Refine ligand occupancy (BUSTER)")
        self.LIGOCC.connect("toggled", self.LIGOCCCallback)
        if self.RefmacParams["LIGOCC"] == "LIGOCC\n":
            self.LIGOCC.set_active(True)
        self.vbox.pack_start(self.LIGOCC, False)

        self.SANITY = gtk.CheckButton("Ignore sanity checks (BUSTER)")
        self.SANITY.connect("toggled", self.SANITYCallback)
        if self.RefmacParams["SANITY"] == "off\n":
            self.SANITY.set_active(True)
        self.vbox.pack_start(self.SANITY, False)

        self.OKbutton = gtk.Button(label="OK")
        self.OKbutton.connect("clicked", self.OK)
        self.vbox.add(self.OKbutton)

        self.window.add(self.vbox)
        self.window.show_all()
        return self.RefmacParams

    def TLSCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["TLS"] = "refi tlsc 10\n"
            self.RefmacParams["TLSIN"] = "refmac.tls\n"
            self.RefmacParams["TLSOUT"] = "out.tls\n"
            self.RefmacParams["TLSADD"] = "TLSO ADDU\n"
        else:
            self.RefmacParams["TLS"] = ""
            self.RefmacParams["TLSIN"] = ""
            self.RefmacParams["TLSOUT"] = ""
            self.RefmacParams["TLSADD"] = ""
        return self.RefmacParams

    def NCSCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["NCS"] = "NCSR LOCAL\n"
        else:
            self.RefmacParams["NCS"] = ""
        return self.RefmacParams

    def ChooseBfacRefinement(self, widget):
        if widget.get_active_text() == "isotropic":
            self.RefmacParams["BREF"] = "    bref ISOT\n"
        if widget.get_active_text() == "anisotropic":
            self.RefmacParams["BREF"] = "    bref ANIS\n"
        return self.RefmacParams

    def on_key_release_Ncycles(self, widget, event):
        print(widget.get_text())
        self.RefmacParams["NCYCLES"] = widget.get_text()
        return self.RefmacParams

    def on_key_release_MATRIX_WEIGHT(self, widget, event):
        self.RefmacParams["MATRIX_WEIGHT"] = widget.get_text()
        return self.RefmacParams

    def TWINCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["TWIN"] = "TWIN\n"
        else:
            self.RefmacParams["TWIN"] = ""
        return self.RefmacParams

    def WATERCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["WATER"] = "WATER\n"
        else:
            self.RefmacParams["WATER"] = ""
        return self.RefmacParams

    def LIGOCCCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["LIGOCC"] = "LIGOCC\n"
        else:
            self.RefmacParams["LIGOCC"] = ""
        return self.RefmacParams

    def SANITYCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["SANITY"] = "off\n"
        else:
            self.RefmacParams["SANITY"] = ""
        return self.RefmacParams

    def OK(self, widget):
        self.window.destroy()

    def ParamsFromPreviousCycle(self, Serial):
        RefmacParams = {
            "HKLIN": "",
            "HKLOUT": "",
            "XYZIN": "",
            "XYZOUT": "",
            "LIBIN": "",
            "LIBOUT": "",
            "TLSIN": "",
            "TLSOUT": "",
            "TLSADD": "",
            "NCYCLES": "10",
            "MATRIX_WEIGHT": "AUTO",
            "BREF": "    bref ISOT\n",
            "TLS": "",
            "NCS": "",
            "TWIN": "",
            "WATER": "",
            "LIGOCC": "",
            "SANITY": "",
        }

        if os.path.isfile(
            self.ProjectPath
            + "/"
            + self.xtalID
            + "/Refine_"
            + str(Serial)
            + "/refmac.csh"
        ):
            for line in open(
                self.ProjectPath
                + "/"
                + self.xtalID
                + "/Refine_"
                + str(Serial)
                + "/refmac.csh"
            ):
                if line.startswith("refi tlsc"):
                    RefmacParams["TLS"] = line
                if line.startswith("TLSO"):
                    RefmacParams["TLSADD"] = line
                if line.startswith("NCSR LOCAL"):
                    RefmacParams["NCS"] = line
                if line.startswith("    bref "):
                    RefmacParams["BREF"] = line
                if line.startswith("ncyc"):
                    RefmacParams["Ncycles"] = line.split()[1]
                if line.startswith("weight"):
                    RefmacParams["MATRIX_WEIGHT"] = line.split()[len(line.split()) - 1]
                if line.startswith("TWIN"):
                    RefmacParams["TWIN"] = line
        elif os.path.isfile(
            os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                "refmac.csh",
            )
        ):
            for line in open(
                os.path.join(
                    self.ProjectPath,
                    self.xtalID,
                    "cootOut",
                    "Refine_" + str(Serial),
                    "refmac.csh",
                )
            ):
                if line.startswith("refi tlsc"):
                    RefmacParams["TLS"] = line
                if line.startswith("TLSO"):
                    RefmacParams["TLSADD"] = line
                if line.startswith("NCSR LOCAL"):
                    RefmacParams["NCS"] = line
                if line.startswith("    bref "):
                    RefmacParams["BREF"] = line
                if line.startswith("ncyc"):
                    RefmacParams["Ncycles"] = line.split()[1]
                if line.startswith("weight"):
                    RefmacParams["MATRIX_WEIGHT"] = line.split()[len(line.split()) - 1]
                if line.startswith("TWIN"):
                    RefmacParams["TWIN"] = line

        return RefmacParams

    def GetRefinementHistory(self):
        RefinementCycle = []
        RcrystList = []
        RfreeList = []

        found = False
        for item in glob.glob(os.path.join(self.ProjectPath, self.xtalID, "*")):
            if item.startswith(os.path.join(self.ProjectPath, self.xtalID, "Refine_")):
                print(item[item.rfind("_") + 1 :])
                RefinementCycle.append(int(item[item.rfind("_") + 1 :]))
                found = True
        if found:
            for cycle in sorted(RefinementCycle):
                try:
                    found_Rcryst = False
                    found_Rfree = False
                    newestPDB = max(
                        glob.iglob(
                            self.ProjectPath
                            + "/"
                            + self.xtalID
                            + "/Refine_"
                            + str(cycle)
                            + "/refine_"
                            + str(cycle)
                            + ".pdb"
                        ),
                        key=os.path.getctime,
                    )
                    for line in open(newestPDB):
                        if line.startswith(
                            "REMARK   3   R VALUE     (WORKING + TEST SET) :"
                        ):
                            Rcryst = line.split()[9]
                            RcrystList.append(float(Rcryst))
                            found_Rcryst = True
                        if line.startswith(
                            "REMARK   3   FREE R VALUE                     :"
                        ):
                            Rfree = line.split()[6]
                            RfreeList.append(float(Rfree))
                            found_Rfree = True
                    if not found_Rcryst:
                        RcrystList.append(0)
                    if not found_Rfree:
                        RfreeList.append(0)
                except ValueError:
                    RcrystList.append(0)
                    RfreeList.append(0)
        else:
            RefinementCycle = [0]
            RcrystList = [0]
            RfreeList = [0]
        print(RefinementCycle, RcrystList, RfreeList)
        return (sorted(RefinementCycle), RcrystList, RfreeList)


class Refine(object):
    def __init__(self, ProjectPath, xtalID, compoundID, datasource):
        self.ProjectPath = ProjectPath
        self.xtalID = xtalID
        self.compoundID = compoundID
        self.datasource = datasource
        self.error = False
        self.Logfile = None

    def GetSerial(self):
        # check if there were already previous refinements
        # if no: create a folder Refine_1
        # if yes: create a folder Refine_<max+1>
        temp = []
        found = 0
        if os.path.isdir(os.path.join(self.ProjectPath, self.xtalID)):
            for item in glob.glob(os.path.join(self.ProjectPath, self.xtalID, "*")):
                if item.startswith(
                    os.path.join(self.ProjectPath, self.xtalID, "Refine_")
                ):
                    print(int(item[item.rfind("_") + 1 :]))
                    temp.append(int(item[item.rfind("_") + 1 :]))
                    found = 1
        if found:
            Serial = max(temp) + 1
        else:
            Serial = 1
        return Serial

    def make_covalent_link(self, covLinkAtomSpec, Logfile):
        residuesToModify = ["CYS", "SER"]

        atom1 = covLinkAtomSpec[1][5]
        atom2 = covLinkAtomSpec[2][5]
        residue_1 = covLinkAtomSpec[3]
        residue_2 = covLinkAtomSpec[4]

        if residue_2 in residuesToModify:
            bond_text = (
                "LINK: RES-NAME-1 %s FILE-1 %s_acedrg.cif ATOM-NAME-1  %s  RES-NAME-2"
                " %s ATOM-NAME-2  %s"
                % (residue_1, self.compoundID, atom1, residue_2, atom2)
            )
        else:
            bond_text = (
                "LINK: RES-NAME-1 %s FILE-1 %s_acedrg.cif ATOM-NAME-1  %s  RES-NAME-2"
                " %s ATOM-NAME-2  %s"
                % (residue_2, self.compoundID, atom2, residue_1, atom1)
            )

        os.chdir(os.path.join(self.ProjectPath, self.xtalID))
        f = open("covalent_bond.txt", "w")
        f.write(bond_text)
        f.close()

        # seems somewhat unncessary, but acedrg does not like some certain descriptions
        # generated by phenix.elbow
        os.system("acedrg -c %s.cif -o %s_acedrg" % (self.compoundID, self.compoundID))

        os.system("acedrg -L covalent_bond.txt -o %s" % self.compoundID)

        if os.path.isfile("%s_link.cif" % self.compoundID):
            Logfile.insert("succefully created link file %s_link.cif" % self.compoundID)
            os.system("/bin/rm %s.cif" % self.compoundID)
            os.system("ln -s %s_link.cif %s.cif" % (self.compoundID, self.compoundID))
        else:
            Logfile.error(
                "could not create linkfile; please check logs in %s/%s"
                % (self.ProjectPath, self.xtalID)
            )

    def get_hklin_hklout(self, Serial):
        #######################################################
        # HKLIN & HKLOUT
        hklin = None
        hklout = None
        if os.path.isfile(
            os.path.join(self.ProjectPath, self.xtalID, self.xtalID + ".free.mtz")
        ):
            hklin = os.path.join(
                self.ProjectPath, self.xtalID, self.xtalID + ".free.mtz"
            )
        else:
            self.Logfile.error("%s: cannot find HKLIN for refinement" % self.xtalID)
            self.error = True
        hklout = os.path.join(
            self.ProjectPath,
            self.xtalID,
            "Refine_" + Serial,
            "refine_" + Serial + ".mtz",
        )
        return hklin, hklout

    def get_xyzin_xyzout(self, Serial):
        #######################################################
        # XYZIN & XYZOUT
        xyzin = None
        xyzout = None
        if os.path.isfile(
            os.path.join(self.ProjectPath, self.xtalID, "Refine_" + Serial, "in.pdb")
        ):
            xyzin = os.path.join(
                self.ProjectPath, self.xtalID, "Refine_" + Serial, "in.pdb"
            )
            xyzout = os.path.join(
                self.ProjectPath,
                self.xtalID,
                "Refine_" + Serial,
                "refine_" + Serial + ".pdb",
            )
        elif os.path.isfile(
            os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                "in.pdb",
            )
        ):
            xyzin = os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                "in.pdb",
            )
        else:
            self.Logfile.error("%s: cannot find XYZIN for refinement" % self.xtalID)
            self.error = True
        return xyzin, xyzout

    def get_libin_libout(self, Serial):
        #######################################################
        # LIBIN & LIBOUT
        libin = None
        libout = None
        if os.path.isfile(
            os.path.join(self.ProjectPath, self.xtalID, self.compoundID + ".cif")
        ):
            libin = os.path.join(
                self.ProjectPath, self.xtalID, self.compoundID + ".cif"
            )
            libout = os.path.join(
                self.ProjectPath,
                self.xtalID,
                "Refine_" + Serial,
                "refine_" + Serial + ".cif",
            )
        else:
            self.Logfile.error("%s: cannot find CIF file for refinement" % self.xtalID)
            self.error = True
        return libin, libout

    def write_refinement_in_progress(self):
        #######################################################
        # we write 'REFINEMENT_IN_PROGRESS' immediately to avoid unncessary refiment
        os.chdir(os.path.join(self.ProjectPath, self.xtalID))
        os.system("touch REFINEMENT_IN_PROGRESS")

    def clean_up_before_refinement(self):
        #######################################################
        # clean up!
        # and remove all files which will be re-created by current refinement cycle
        os.chdir(os.path.join(self.ProjectPath, self.xtalID))
        os.system(
            "/bin/rm refine.pdb refine.mtz refine.split.bound-state.pdb"
            " validation_summary.txt validate_ligands.txt 2fofc.map fofc.map"
            " refine_molprobity.log"
        )

    def get_shebang(self):
        return "#!" + os.getenv("SHELL") + "\n"

    def load_modules(self, cmd):
        #######################################################
        # PHENIX stuff (if working at DLS)
        if os.getcwd().startswith("/dls"):
            cmd += "module load phenix/1.20\n"
            cmd += "module load buster/20240123\n"
        return cmd

    def get_source_line(self, cmd):
        if "bash" in os.getenv("SHELL"):
            cmd += 'export XChemExplorer_DIR="' + os.getenv("XChemExplorer_DIR") + '"\n'
        return cmd

    def set_refinement_status(self, cmd):
        cmd += (
            "$CCP4/bin/ccp4-python "
            " $XChemExplorer_DIR/xce/helpers/update_status_flag.py "
            + self.datasource
            + " "
            + self.xtalID
            + " RefinementStatus running\n"
        )
        return cmd

    def add_buster_command(
        self,
        cmd,
        xyzin,
        hklin,
        libin,
        Serial,
        anisotropic_Bfactor,
        update_water,
        refine_ligand_occupancy,
        ignore_sanity_check,
    ):
        cmd += (
            "refine "
            " -p %s" % xyzin
            + " -m %s" % hklin
            + " -l %s" % libin
            + " -autoncs"
            + anisotropic_Bfactor
            + update_water
            + refine_ligand_occupancy
            + ignore_sanity_check
            + " -d Refine_%s\n" % str(Serial)
        )
        return cmd

    def run_giant_score_model(self, cmd, cycle):
        if os.path.isfile(
            os.path.join(
                self.ProjectPath, self.xtalID, self.xtalID + "-pandda-model.pdb"
            )
        ):
            cmd += (
                "cd "
                + self.ProjectPath
                + "/"
                + self.xtalID
                + "/Refine_"
                + cycle
                + "\n"
                + "module load ccp4/7.1.018\n"
                + "giant.score_model "
                " pdb1=../%s-pandda-model.pdb " % self.xtalID + " mtz1=../dimple.mtz "
                " pdb2=refine.pdb "
                " mtz2=refine.mtz\n"
            )
        return cmd

    def add_validation(self, cmd, cycle, program):
        buster_report = ""
        if program == "refmac":
            Serial = "_" + cycle
        elif program == "buster":
            Serial = ""
            if os.getcwd().startswith("/dls"):
                buster_report += (
                    "cd "
                    + self.ProjectPath
                    + "/"
                    + self.xtalID
                    + "/Refine_"
                    + cycle
                    + "\n"
                )
                buster_report += (
                    "corr -p refine.pdb -m refine.mtz -F 2FOFCWT -P PH2FOFCWT\n"
                )
                buster_report += "cd " + self.ProjectPath + "/" + self.xtalID + "\n"
                buster_report += (
                    "export BDG_TOOL_MOGUL=/dls_sw/apps/ccdc/CSD_2020/bin/mogul\n"
                )
                buster_report += "buster-report -d Refine_%s\n" % cycle

        cmd += (
            "\n" + buster_report + "cd " + self.ProjectPath + "/" + self.xtalID + "\n"
            "ln -s ./Refine_%s/refine%s.pdb refine.pdb\n" % (cycle, Serial)
            + "ln -s ./Refine_%s/refine%s.mtz refine.mtz\n" % (cycle, Serial)
            + "ln -s refine.pdb refine.split.bound-state.pdb\n"
            "\n"
            "phenix.molprobity refine%s.pdb refine%s.mtz\n" % (Serial, Serial)
            + "/bin/mv molprobity.out refine_molprobity.log\n"
            "mmtbx.validate_ligands refine%s.pdb refine%s.mtz LIG"
            " > validate_ligands.txt\n" % (Serial, Serial) + "\n"
            "mmtbx.validation_summary refine.pdb > validation_summary.txt\n"
            "\n"
        )
        return cmd

    def calculate_maps(self, cmd, program):
        if program == "refmac":
            FWT = "FWT"
            PHWT = "PHWT"
            DELFWT = "DELFWT"
            PHDELWT = "PHDELWT"
        elif program == "buster":
            FWT = "2FOFCWT"
            PHWT = "PH2FOFCWT"
            DELFWT = "FOFCWT"
            PHDELWT = "PHFOFCWT"
        cmd += (
            "cd " + self.ProjectPath + "/" + self.xtalID + "\n"
            "\n"
            "fft hklin refine.mtz mapout 2fofc.map << EOF\n"
            " labin F1=%s PHI=%s\n" % (FWT, PHWT) + " grid samp 4.5\n"
            "EOF\n"
            "\n"
            "fft hklin refine.mtz mapout fofc.map << EOF\n"
            " labin F1=%s PHI=%s\n" % (DELFWT, PHDELWT) + " grid samp 4.5\n"
            "EOF\n"
            "\n"
        )
        return cmd

    def run_edstats(self, cmd, resh, resl):
        cmd += (
            "\n"
            "edstats XYZIN refine.pdb MAPIN1 2fofc.map MAPIN2 fofc.map OUT"
            " refine.edstats << eof\n"
            " resl={0!s}\n".format(resl) + " resh={0!s}\n".format(resh) + "eof\n"
            "\n"
        )
        return cmd

    def update_database(self, cmd, Serial):
        date = datetime.strftime(datetime.now(), "%Y-%m-%d_%H-%M-%S.%f")[:-4]
        user = getpass.getuser()
        ccp4_module = ""  # need to source again, because giant.score_model still needs
        # outdated version which does not have gemmi
        if self.ProjectPath.startswith("/dls"):
            ccp4_module = "module load ccp4/7.1.018"
        cmd += "\n" + ccp4_module + "\n" "$CCP4/bin/ccp4-python " + os.path.join(
            os.getenv("XChemExplorer_DIR"),
            "xce",
            "helpers",
            "update_data_source_after_refinement.py",
        ) + " %s %s %s %s %s %s\n" % (
            self.datasource,
            self.xtalID,
            self.ProjectPath,
            os.path.join(self.ProjectPath, self.xtalID, "Refine_" + Serial),
            user,
            date,
        )
        cmd += "/bin/rm %s/%s/REFINEMENT_IN_PROGRESS\n" % (
            self.ProjectPath,
            self.xtalID,
        )

        return cmd

    def write_refinement_script(self, cmd, program):
        os.chdir(os.path.join(self.ProjectPath, self.xtalID))
        f = open(program + ".sh", "w")
        f.write(cmd)
        f.close()

    def run_script(self, program, external_software, Serial, xce_logfile):
        os.chdir(os.path.join(self.ProjectPath, self.xtalID))
        if external_software["slurm"]:
            slurm.submit_cluster_job(
                "xce_{!s}".format(program), "{!s}.sh".format(program), xce_logfile
            )
        elif external_software["qsub"]:
            sge.submit_cluster_job(
                "xce_{!s}".format(program), "{!s}.sh".format(program), xce_logfile
            )
        else:
            self.Logfile.insert("starting refinement with command: ./%s.sh &" % program)
            os.system("chmod +x %s.sh" % program)
            os.system("./%s.sh &" % program)

    def RunBuster(
        self, Serial, RefmacParams, external_software, xce_logfile, covLinkAtomSpec
    ):
        if RefmacParams is None:
            anisotropic_Bfactor = " -M ADP "
            update_water = ""
        else:
            if "ANIS" in RefmacParams["BREF"]:
                anisotropic_Bfactor = " -M ADP "
            else:
                anisotropic_Bfactor = " -M TLSbasic "

            if "WATER" in RefmacParams["WATER"]:
                update_water = " -WAT "
            else:
                update_water = ""

            if "off" in RefmacParams["SANITY"]:
                ignore_sanity_check = " StopOnGellySanityCheckError=no "
            else:
                ignore_sanity_check = ""

        self.error = False

        if os.path.isfile(xce_logfile):
            self.Logfile = XChemLog.updateLog(xce_logfile)
        Serial = str(Serial)

        if covLinkAtomSpec is not None:
            self.make_covalent_link(covLinkAtomSpec, self.Logfile)

        # first check if refinement is ongoing and exit if yes
        if os.path.isfile(
            os.path.join(self.ProjectPath, self.xtalID, "REFINEMENT_IN_PROGRESS")
        ):
            self.Logfile.insert(
                "cannot start new refinement for %s: *** REFINEMENT IN PROGRESS ***"
                % self.xtalID
            )
            return None

        hklin, hklout = self.get_hklin_hklout(Serial)

        resh, resl = XChemUtils.mtztools_gemmi(hklin).get_high_low_resolution_limits()

        xyzin, xyzout = self.get_xyzin_xyzout(Serial)

        refine_ligand_occupancy = ""
        if RefmacParams is not None:
            if "LIGOCC" in RefmacParams["LIGOCC"]:
                ligand_info = XChemUtils.pdbtools(xyzin).get_residues_with_resname(
                    "LIG"
                )
                if self.prepare_gelly_dat(ligand_info):
                    refine_ligand_occupancy = " -B user -Gelly gelly.dat "

        libin, libout = self.get_libin_libout(Serial)

        self.clean_up_before_refinement()

        self.write_refinement_in_progress()

        cmd = self.get_shebang()

        cmd = self.load_modules(cmd)

        cmd = self.get_source_line(cmd)

        cmd = self.set_refinement_status(cmd)

        cmd = self.add_buster_command(
            cmd,
            xyzin,
            hklin,
            libin,
            Serial,
            anisotropic_Bfactor,
            update_water,
            refine_ligand_occupancy,
            ignore_sanity_check,
        )

        cmd = self.add_validation(cmd, Serial, "buster")

        cmd = self.run_giant_score_model(cmd, Serial)

        cmd = self.calculate_maps(cmd, "buster")

        cmd = self.run_edstats(cmd, resh, resl)

        cmd = self.update_database(cmd, Serial)

        if self.error:
            self.Logfile.error(
                "%s: cannot start refinement; check error message above" % self.xtalID
            )
        else:
            self.Logfile.insert(
                "%s: writing buster.sh file in %s"
                % (self.xtalID, os.path.join(self.ProjectPath, self.xtalID))
            )
            self.write_refinement_script(cmd, "buster")
            self.Logfile.insert("%s: starting refinement..." % self.xtalID)
            self.run_script("buster", external_software, Serial, xce_logfile)

    def prepare_gelly_dat(self, ligand_info):
        found_ligand = False
        gelly = "NOTE BUSTER_SET Ligand = "
        for n, lig in enumerate(ligand_info):
            resseq = lig[1]
            chain = lig[2]
            if n == 0:
                gelly += "{ %s|%s }" % (chain.replace(" ", ""), resseq.replace(" ", ""))
            else:
                gelly += " + { %s|%s }" % (
                    chain.replace(" ", ""),
                    resseq.replace(" ", ""),
                )
            found_ligand = True
        if "{" in gelly:
            gelly += "\n"
            gelly += (
                "NOTE BUSTER_SET NotLigand = \\ Ligand\n"
                "NOTE BUSTER_CONSTANT OCC NotLigand \n"
                "NOTE BUSTER_COMBINE OCC Ligand \n"
                "NOTE BUSTER_COMBINE B Ligand      \n"
            )
            os.chdir(os.path.join(self.ProjectPath, self.xtalID))
            f = open("gelly.dat", "w")
            f.write(gelly)
            f.close()
        return found_ligand

    def RunRefmac(
        self, Serial, RefmacParams, external_software, xce_logfile, covLinkAtomSpec
    ):
        if os.path.isfile(xce_logfile):
            Logfile = XChemLog.updateLog(xce_logfile)
        Serial = str(Serial)

        if covLinkAtomSpec is not None:
            self.make_covalent_link(covLinkAtomSpec, Logfile)

        # first check if refinement is ongoing and exit if yes
        if os.path.isfile(
            os.path.join(self.ProjectPath, self.xtalID, "REFINEMENT_IN_PROGRESS")
        ):
            #            coot.info_dialog('*** REFINEMENT IN PROGRESS ***')
            Logfile.insert(
                "cannot start new refinement for %s: *** REFINEMENT IN PROGRESS ***"
                % self.xtalID
            )
            return None

        #######################################################
        # HKLIN & HKLOUT
        if os.path.isfile(
            os.path.join(self.ProjectPath, self.xtalID, self.xtalID + ".free.mtz")
        ):
            RefmacParams["HKLIN"] = "HKLIN " + os.path.join(
                self.ProjectPath, self.xtalID, self.xtalID + ".free.mtz \\\n"
            )
        else:
            Logfile.error(
                "%s: cannot find HKLIN for refinement; aborting..." % self.xtalID
            )
            return None
        RefmacParams["HKLOUT"] = "HKLOUT " + os.path.join(
            self.ProjectPath,
            self.xtalID,
            "Refine_" + Serial,
            "refine_" + Serial + ".mtz \\\n",
        )

        #######################################################
        # XYZIN & XYZOUT
        RefmacParams["XYZIN"] = "XYZIN " + os.path.join(
            self.ProjectPath, self.xtalID, "Refine_" + Serial, "in.pdb \\\n"
        )
        RefmacParams["XYZOUT"] = "XYZOUT " + os.path.join(
            self.ProjectPath,
            self.xtalID,
            "Refine_" + Serial,
            "refine_" + Serial + ".pdb \\\n",
        )

        #######################################################
        # LIBIN & LIBOUT
        if os.path.isfile(
            os.path.join(self.ProjectPath, self.xtalID, self.compoundID + ".cif")
        ):
            RefmacParams["LIBIN"] = (
                "LIBIN "
                + self.ProjectPath
                + "/"
                + self.xtalID
                + "/"
                + self.compoundID
                + ".cif \\\n"
            )
            RefmacParams["LIBOUT"] = (
                "LIBOUT "
                + self.ProjectPath
                + "/"
                + self.xtalID
                + "/Refine_"
                + Serial
                + "/refine_"
                + Serial
                + ".cif \\\n"
            )

        #######################################################
        # TLSIN & TLSOUT
        findTLS = "\n"
        if RefmacParams["TLS"].startswith("refi"):
            if external_software["phenix.find_tls_groups"]:
                findTLS = (
                    os.path.join(
                        os.getenv("XChemExplorer_DIR"),
                        "xce",
                        "helpers",
                        "phenix_find_TLS_groups.py",
                    )
                    + " in.pdb\n"
                )
                RefmacParams["TLSIN"] = (
                    "TLSIN "
                    + self.ProjectPath
                    + "/"
                    + self.xtalID
                    + "/Refine_"
                    + Serial
                    + "/refmac.tls \\\n"
                )
                RefmacParams["TLSOUT"] = (
                    "TLSOUT "
                    + self.ProjectPath
                    + "/"
                    + self.xtalID
                    + "/Refine_"
                    + Serial
                    + "/refine.tls \\\n"
                )
            else:
                RefmacParams["TLS"] = "\n"

        print("==> XCE: assembling refmac.csh")

        #######################################################
        # we write 'REFINEMENT_IN_PROGRESS' immediately to avoid unncessary refiment
        os.chdir(os.path.join(self.ProjectPath, self.xtalID))
        os.system("touch REFINEMENT_IN_PROGRESS")

        #######################################################
        # Database updates:
        # no DB will be specified when a reference model is built and refined
        refinementStatus = ""
        updateDB = ""
        date = datetime.strftime(datetime.now(), "%Y-%m-%d_%H-%M-%S.%f")[:-4]
        user = getpass.getuser()
        if os.path.isfile(self.datasource):
            refinementStatus = (
                "$CCP4/bin/ccp4-python"
                " $XChemExplorer_DIR/xce/helpers/update_status_flag.py %s %s %s %s\n"
                % (self.datasource, self.xtalID, "RefinementStatus", "running")
            )
            updateDB = (
                "$CCP4/bin/ccp4-python "
                + os.path.join(
                    os.getenv("XChemExplorer_DIR"),
                    "xce",
                    "helpers",
                    "update_data_source_after_refinement.py",
                )
                + " %s %s %s %s %s %s\n"
                % (
                    self.datasource,
                    self.xtalID,
                    self.ProjectPath,
                    os.path.join(self.ProjectPath, self.xtalID, "Refine_" + Serial),
                    user,
                    date,
                )
            )

        #######################################################
        # clean up!
        # and remove all files which will be re-created by current refinement cycle
        os.system(
            "/bin/rm refine.pdb refine.mtz refine.split.bound-state.pdb"
            " validation_summary.txt validate_ligands.txt 2fofc.map fofc.map"
            " refine_molprobity.log"
        )

        #######################################################
        # weight
        if str(RefmacParams["MATRIX_WEIGHT"]).lower() == "auto":
            weight = "weight AUTO\n"
        else:
            weight = "weight matrix " + str(RefmacParams["MATRIX_WEIGHT"]) + "\n"

        #######################################################
        # PHENIX stuff

        spider_plot = ""
        if os.path.isfile(
            os.path.join(
                self.ProjectPath, self.xtalID, self.xtalID + "-ensemble-model.pdb"
            )
        ):
            if os.path.isfile(
                os.path.join(
                    self.ProjectPath, self.xtalID, self.xtalID + "-pandda-input.mtz"
                )
            ):
                pdb_two = os.path.join(
                    self.ProjectPath, self.xtalID, self.xtalID + "-ensemble-model.pdb"
                )
                mtz_two = os.path.join(
                    self.ProjectPath, self.xtalID, self.xtalID + "-pandda-input.mtz"
                )
                pdb_one = os.path.join(
                    self.ProjectPath,
                    self.xtalID,
                    "Refine_" + str(Serial),
                    "refine_" + str(Serial) + ".pdb",
                )
                mtz_one = os.path.join(
                    self.ProjectPath,
                    self.xtalID,
                    "Refine_" + str(Serial),
                    "refine_" + str(Serial) + ".mtz",
                )
                spider_plot = (
                    "cd "
                    + self.ProjectPath
                    + "/"
                    + self.xtalID
                    + "/Refine_"
                    + Serial
                    + "\n"
                    "\n"
                    "giant.score_model pdb1=%s mtz1=%s pdb2=%s mtz2=%s"
                    " res_names=LIG,UNL,DRG,FRG\n"
                    % (pdb_one, mtz_one, pdb_two, mtz_two)
                )

        refmacCmds = (
            "#!"
            + os.getenv("SHELL")
            + "\n"
            + +'export XChemExplorer_DIR="'
            + os.getenv("XChemExplorer_DIR")
            + '"\n'
            + "module load phenix/1.20\n"
            + "module load buster/20240123\n"
            + "cd "
            + self.ProjectPath
            + "/"
            + self.xtalID
            + "/Refine_"
            + Serial
            + "\n"
            "\n"
            + refinementStatus
            + "\n"
            + findTLS
            + "refmac5 "
            + RefmacParams["HKLIN"]
            + RefmacParams["HKLOUT"]
            + RefmacParams["XYZIN"]
            + RefmacParams["XYZOUT"]
            + RefmacParams["LIBIN"]
            + RefmacParams["LIBOUT"]
            + RefmacParams["TLSIN"]
            + RefmacParams["TLSOUT"]
            + " << EOF > refmac.log\n"
            "make -\n"
            "    hydrogen ALL -\n"
            "    hout NO -\n"
            "    peptide NO -\n"
            "    cispeptide YES -\n"
            "    ssbridge YES -\n"
            "    symmetry YES -\n"
            "    sugar YES -\n"
            "    connectivity NO -\n"
            "    link NO\n" + RefmacParams["NCS"] + "refi -\n"
            "    type REST -\n"
            "    resi MLKF -\n"
            "    meth CGMAT -\n"
            + RefmacParams["BREF"]
            + RefmacParams["TLS"]
            + RefmacParams["TWIN"]
            + "ncyc "
            + RefmacParams["NCYCLES"]
            + "\n"
            "scal -\n"
            "    type SIMP -\n"
            "    LSSC -\n"
            "    ANISO -\n"
            "    EXPE\n" + weight + "solvent YES\n"
            "monitor MEDIUM -\n"
            "    torsion 10.0 -\n"
            "    distance 10.0 -\n"
            "    angle 10.0 -\n"
            "    plane 10.0 -\n"
            "    chiral 10.0 -\n"
            "    bfactor 10.0 -\n"
            "    bsphere 10.0 -\n"
            "    rbond 10.0 -\n"
            "    ncsr 10.0\n"
            "labin  FP=F SIGFP=SIGF FREE=FreeR_flag\n"
            "labout  FC=FC FWT=FWT PHIC=PHIC PHWT=PHWT DELFWT=DELFWT PHDELWT=PHDELWT"
            " FOM=FOM\n" + RefmacParams["TLSADD"] + "\n"
            "DNAME " + self.xtalID + "\n"
            "END\n"
            "EOF\n"
            "\n"
            "phenix.molprobity refine_%s.pdb refine_%s.mtz\n" % (Serial, Serial)
            + "/bin/mv molprobity.out refine_molprobity.log\n"
            "mmtbx.validate_ligands refine_%s.pdb refine_%s.mtz LIG"
            " > validate_ligands.txt\n" % (Serial, Serial)
            + "cd "
            + self.ProjectPath
            + "/"
            + self.xtalID
            + "\n"
            "#ln -s %s/%s/Refine_%s/refine_%s.pdb refine.pdb\n"
            % (self.ProjectPath, self.xtalID, Serial, Serial)
            + "#ln -s %s/%s/Refine_%s/refine_%s.mtz refine.mtz\n"
            % (self.ProjectPath, self.xtalID, Serial, Serial)
            + "ln -s ./Refine_%s/refine_%s.pdb refine.pdb\n" % (Serial, Serial)
            + "ln -s ./Refine_%s/refine_%s.mtz refine.mtz\n" % (Serial, Serial)
            + "ln -s refine.pdb refine.split.bound-state.pdb\n"
            "\n"
            "ln -s Refine_%s/validate_ligands.txt .\n" % Serial
            + "ln -s Refine_%s/refine_molprobity.log .\n" % Serial
            + "mmtbx.validation_summary refine.pdb > validation_summary.txt\n"
            "\n"
            "fft hklin refine.mtz mapout 2fofc.map << EOF\n"
            "labin F1=FWT PHI=PHWT\n"
            "EOF\n"
            "\n"
            "fft hklin refine.mtz mapout fofc.map << EOF\n"
            "labin F1=DELFWT PHI=PHDELWT\n"
            "EOF\n"
            "\n" + updateDB + "\n"
            "/bin/rm %s/%s/REFINEMENT_IN_PROGRESS\n" % (self.ProjectPath, self.xtalID)
            + "\n"
            "cd " + self.ProjectPath + "/" + self.xtalID + "/Refine_" + Serial + "\n"
            "\n" + spider_plot + "\n"
        )

        if os.path.isfile(xce_logfile):
            Logfile.insert(
                "writing refinement shell script to"
                + os.path.join(
                    self.ProjectPath, self.xtalID, "Refine_" + Serial, "refmac.csh"
                )
            )
        cmd = open(
            os.path.join(
                self.ProjectPath, self.xtalID, "Refine_" + Serial, "refmac.csh"
            ),
            "w",
        )
        cmd.write(refmacCmds)
        cmd.close()

        os.chdir(os.path.join(self.ProjectPath, self.xtalID, "Refine_" + Serial))
        Logfile.insert(
            "changing directory to %s"
            % (os.path.join(self.ProjectPath, self.xtalID, "Refine_" + Serial))
        )

        if external_software["slurm"]:
            slurm.submit_cluster_job("xce_refmac", "refmac.csh", xce_logfile)
        elif external_software["qsub"]:
            sge.submit_cluster_job("xce_refmac", "refmac.csh", xce_logfile)
        else:
            os.system("chmod +x refmac.csh")
            if os.path.isfile(xce_logfile):
                Logfile.insert("starting refinement on local machine")
            os.system("./refmac.csh &")

    def RefinementParams(self, RefmacParams):
        self.RefmacParams = RefmacParams
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", gtk.main_quit)
        self.window.set_border_width(10)
        self.window.set_title("Refmac Parameters")
        self.vbox = gtk.VBox()

        self.hbox1 = gtk.HBox()
        self.hbox1.add(gtk.Label("Refine"))
        self.cb = gtk.combo_box_new_text()
        self.cb.connect("changed", self.ChooseBfacRefinement)
        for item in ["isotropic", "anisotropic"]:
            self.cb.append_text(item)
        if "ISOT" in self.RefmacParams["BREF"]:
            self.cb.set_active(0)
        if "ANIS" in self.RefmacParams["BREF"]:
            self.cb.set_active(1)
        self.hbox1.add(self.cb)
        self.hbox1.add(gtk.Label("temperature factors"))
        self.vbox.add(self.hbox1)

        self.hbox2 = gtk.HBox()
        self.hbox2.add(gtk.Label("Number of Cycles: "))
        self.Ncycles = gtk.Entry()
        self.Ncycles.add_events(gtk.gdk.KEY_RELEASE_MASK)
        self.Ncycles.connect("key-release-event", self.on_key_release_Ncycles)
        self.Ncycles.set_text(self.RefmacParams["NCYCLES"])
        self.hbox2.add(self.Ncycles)
        self.vbox.add(self.hbox2)

        self.hbox3 = gtk.HBox()
        self.hbox3.add(gtk.Label("MATRIX WEIGHT: "))
        self.MATRIX_WEIGHT = gtk.Entry()
        self.MATRIX_WEIGHT.add_events(gtk.gdk.KEY_RELEASE_MASK)
        self.MATRIX_WEIGHT.connect(
            "key-release-event", self.on_key_release_MATRIX_WEIGHT
        )
        self.MATRIX_WEIGHT.set_text(self.RefmacParams["MATRIX_WEIGHT"])
        self.hbox3.add(self.MATRIX_WEIGHT)
        self.vbox.add(self.hbox3)

        self.TLS = gtk.CheckButton("TLS (find TLS groups with phenix.find_tls_groups)")
        self.TLS.connect("toggled", self.TLSCallback)
        if self.RefmacParams["TLS"] == "refi tlsc 10\n":
            self.TLS.set_active(True)
        self.vbox.pack_start(self.TLS, False)

        self.NCS = gtk.CheckButton("NCS (if applicable")
        self.NCS.connect("toggled", self.NCSCallback)
        if self.RefmacParams["NCS"] == "NCSR LOCAL\n":
            self.NCS.set_active(True)
        self.vbox.pack_start(self.NCS, False)

        self.TWIN = gtk.CheckButton("Twin?")
        self.TWIN.connect("toggled", self.TWINCallback)
        if self.RefmacParams["TWIN"] == "TWIN\n":
            self.TWIN.set_active(True)
        self.vbox.pack_start(self.TWIN, False)

        self.OKbutton = gtk.Button(label="OK")
        self.OKbutton.connect("clicked", self.OK)
        self.vbox.add(self.OKbutton)

        self.window.add(self.vbox)
        self.window.show_all()
        return self.RefmacParams

    def TLSCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["TLS"] = "refi tlsc 10\n"
            self.RefmacParams["TLSIN"] = "refmac.tls\n"
            self.RefmacParams["TLSOUT"] = "out.tls\n"
            self.RefmacParams["TLSADD"] = "TLSO ADDU\n"
        else:
            self.RefmacParams["TLS"] = ""
            self.RefmacParams["TLSIN"] = ""
            self.RefmacParams["TLSOUT"] = ""
            self.RefmacParams["TLSADD"] = ""
        return self.RefmacParams

    def NCSCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["NCS"] = "NCSR LOCAL\n"
        else:
            self.RefmacParams["NCS"] = ""
        return self.RefmacParams

    def ChooseBfacRefinement(self, widget):
        if widget.get_active_text() == "isotropic":
            self.RefmacParams["BREF"] = "    bref ISOT\n"
        if widget.get_active_text() == "anisotropic":
            self.RefmacParams["BREF"] = "    bref ANIS\n"
        return self.RefmacParams

    def on_key_release_Ncycles(self, widget, event):
        print(widget.get_text())
        self.RefmacParams["NCYCLES"] = widget.get_text()
        return self.RefmacParams

    def on_key_release_MATRIX_WEIGHT(self, widget, event):
        self.RefmacParams["MATRIX_WEIGHT"] = widget.get_text()
        return self.RefmacParams

    def TWINCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["TWIN"] = "TWIN\n"
        else:
            self.RefmacParams["TWIN"] = ""
        return self.RefmacParams

    def OK(self, widget):
        self.window.destroy()

    def ParamsFromPreviousCycle(self, Serial):
        RefmacParams = {
            "HKLIN": "",
            "HKLOUT": "",
            "XYZIN": "",
            "XYZOUT": "",
            "LIBIN": "",
            "LIBOUT": "",
            "TLSIN": "",
            "TLSOUT": "",
            "TLSADD": "",
            "NCYCLES": "10",
            "MATRIX_WEIGHT": "AUTO",
            "BREF": "    bref ISOT\n",
            "TLS": "",
            "NCS": "",
            "TWIN": "",
            "WATER": "",
            "LIGOCC": "",
            "SANITY": "",
        }

        if os.path.isfile(
            self.ProjectPath
            + "/"
            + self.xtalID
            + "/Refine_"
            + str(Serial)
            + "/refmac.csh"
        ):
            for line in open(
                self.ProjectPath
                + "/"
                + self.xtalID
                + "/Refine_"
                + str(Serial)
                + "/refmac.csh"
            ):
                if line.startswith("refi tlsc"):
                    RefmacParams["TLS"] = line
                if line.startswith("TLSO"):
                    RefmacParams["TLSADD"] = line
                if line.startswith("NCSR LOCAL"):
                    RefmacParams["NCS"] = line
                if line.startswith("    bref "):
                    RefmacParams["BREF"] = line
                if line.startswith("ncyc"):
                    RefmacParams["Ncycles"] = line.split()[1]
                if line.startswith("weight"):
                    RefmacParams["MATRIX_WEIGHT"] = line.split()[len(line.split()) - 1]
                if line.startswith("TWIN"):
                    RefmacParams["TWIN"] = line
        elif os.path.isfile(
            os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                "multi-state-restraints.refmac.params",
            )
        ):
            for line in open(
                os.path.join(
                    self.ProjectPath,
                    self.xtalID,
                    "cootOut",
                    "Refine_" + str(Serial),
                    "multi-state-restraints.refmac.params",
                )
            ):
                if "refi tlsc" in line:
                    RefmacParams["TLS"] = line
                if "TLSO" in line:
                    RefmacParams["TLSADD"] = line
                if "NCSR LOCAL" in line:
                    RefmacParams["NCS"] = line
                if "bref" in line:
                    RefmacParams["BREF"] = line
                if "ncyc" in line:
                    RefmacParams["Ncycles"] = line.split()[1]
                if "weight" in line:
                    RefmacParams["MATRIX_WEIGHT"] = line.split()[len(line.split()) - 1]
                if "TWIN" in line:
                    RefmacParams["TWIN"] = line

        return RefmacParams

    def GetRefinementHistory(self):
        RefinementCycle = []
        RcrystList = []
        RfreeList = []

        found = False
        for item in glob.glob(os.path.join(self.ProjectPath, self.xtalID, "*")):
            if item.startswith(os.path.join(self.ProjectPath, self.xtalID, "Refine_")):
                print(item[item.rfind("_") + 1 :])
                RefinementCycle.append(int(item[item.rfind("_") + 1 :]))
                found = True
        if found:
            for cycle in sorted(RefinementCycle):
                try:
                    found_Rcryst = False
                    found_Rfree = False
                    newestPDB = max(
                        glob.iglob(
                            self.ProjectPath
                            + "/"
                            + self.xtalID
                            + "/Refine_"
                            + str(cycle)
                            + "/refine_"
                            + str(cycle)
                            + ".pdb"
                        ),
                        key=os.path.getctime,
                    )
                    for line in open(newestPDB):
                        if line.startswith(
                            "REMARK   3   R VALUE     (WORKING + TEST SET) :"
                        ):
                            Rcryst = line.split()[9]
                            RcrystList.append(float(Rcryst))
                            found_Rcryst = True
                        if line.startswith(
                            "REMARK   3   FREE R VALUE                     :"
                        ):
                            Rfree = line.split()[6]
                            RfreeList.append(float(Rfree))
                            found_Rfree = True
                    if not found_Rcryst:
                        RcrystList.append(0)
                    if not found_Rfree:
                        RfreeList.append(0)
                except ValueError:
                    RcrystList.append(0)
                    RfreeList.append(0)
        else:
            RefinementCycle = [0]
            RcrystList = [0]
            RfreeList = [0]
        print(RefinementCycle, RcrystList, RfreeList)
        return (sorted(RefinementCycle), RcrystList, RfreeList)


class panddaRefine(object):
    def __init__(self, ProjectPath, xtalID, compoundID, datasource):
        self.ProjectPath = ProjectPath
        self.xtalID = xtalID
        self.compoundID = compoundID
        self.datasource = datasource

    def make_covalent_link(self, covLinkAtomSpec, Logfile):
        residuesToModify = ["CYS", "SER"]

        atom1 = covLinkAtomSpec[1][5]
        atom2 = covLinkAtomSpec[2][5]
        residue_1 = covLinkAtomSpec[3]
        residue_2 = covLinkAtomSpec[4]

        if residue_2 in residuesToModify:
            bond_text = (
                "LINK: RES-NAME-1 %s FILE-1 %s_acedrg.cif ATOM-NAME-1  %s  RES-NAME-2"
                " %s ATOM-NAME-2  %s"
                % (residue_1, self.compoundID, atom1, residue_2, atom2)
            )
        else:
            bond_text = (
                "LINK: RES-NAME-1 %s FILE-1 %s_acedrg.cif ATOM-NAME-1  %s  RES-NAME-2"
                " %s ATOM-NAME-2  %s"
                % (residue_2, self.compoundID, atom2, residue_1, atom1)
            )

        os.chdir(os.path.join(self.ProjectPath, self.xtalID))
        f = open("covalent_bond.txt", "w")
        f.write(bond_text)
        f.close()

        # seems somewhat unncessary, but acedrg does not like some certain descriptions
        # generated by phenix.elbow
        os.system("acedrg -c %s.cif -o %s_acedrg" % (self.compoundID, self.compoundID))

        os.system("acedrg -L covalent_bond.txt -o %s" % self.compoundID)

        if os.path.isfile("%s_link.cif" % self.compoundID):
            Logfile.insert("succefully created link file %s_link.cif" % self.compoundID)
            os.system("/bin/rm %s.cif" % self.compoundID)
            os.system("ln -s %s_link.cif %s.cif" % (self.compoundID, self.compoundID))
        else:
            Logfile.error(
                "could not create linkfile; please check logs in %s/%s"
                % (self.ProjectPath, self.xtalID)
            )

    def RunQuickRefine(
        self,
        Serial,
        RefmacParams,
        external_software,
        xce_logfile,
        refinementProtocol,
        covLinkAtomSpec,
    ):
        Logfile = XChemLog.updateLog(xce_logfile)
        Logfile.insert("preparing files for giant.quick_refine")
        # panddaSerial because giant.quick_refine writes Refine_0001 instead of Refine_1
        panddaSerial = (4 - len(str(Serial))) * "0" + str(Serial)

        make_all_links = True
        add_links_line = ""

        if covLinkAtomSpec is not None:
            self.make_covalent_link(covLinkAtomSpec, Logfile)

        # first check if refinement is ongoing and exit if yes
        if os.path.isfile(
            os.path.join(self.ProjectPath, self.xtalID, "REFINEMENT_IN_PROGRESS")
        ):
            Logfile.insert(
                "cannot start new refinement for %s: *** REFINEMENT IN PROGRESS ***"
                % self.xtalID
            )
            return None

        #######################################################
        # HKLIN & HKLOUT
        if os.path.isfile(
            os.path.join(self.ProjectPath, self.xtalID, self.xtalID + ".free.mtz")
        ):
            RefmacParams["HKLIN"] = os.path.join(
                self.ProjectPath, self.xtalID, self.xtalID + ".free.mtz"
            )
        elif os.path.isfile(
            os.path.join(
                self.ProjectPath, self.xtalID, self.xtalID + "-pandda-input.mtz"
            )
        ):
            RefmacParams["HKLIN"] = os.path.join(
                self.ProjectPath, self.xtalID, self.xtalID + "-pandda-input.mtz"
            )
        else:
            Logfile.error(
                "%s: cannot find HKLIN for refinement; aborting..." % self.xtalID
            )
            return None

        #######################################################
        # LIBIN & LIBOUT
        found_cif = False
        if os.path.isfile(
            os.path.join(self.ProjectPath, self.xtalID, self.compoundID + ".cif")
        ):
            found_cif = True
            RefmacParams["LIBIN"] = os.path.join(
                self.ProjectPath, self.xtalID, self.compoundID + ".cif"
            )
        if not found_cif:
            # this should actually not be necessary, but the following scenario can
            # happen if a new data source is created from a file system, but smiles and
            # compoundID where not updated; so the ligand may still be in the structure,
            # but since the compoundID is unknown to the datasource its restraints
            # won't be read in and refmac will fail
            for file in glob.glob(os.path.join(self.ProjectPath, self.xtalID, "*")):
                if file.endswith(".cif"):
                    RefmacParams["LIBIN"] = file
                    break

        #######################################################
        # giant_merge_conformations
        Logfile.insert("trying to merge modified bound state with ground state")
        if os.path.isfile(
            os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                self.xtalID + "-ensemble-model.pdb",
            )
        ):
            Logfile.insert(
                "seems to be an initial refinement after pandda.export,"
                " no need to merge the conformations"
            )
            os.chdir(
                os.path.join(
                    self.ProjectPath, self.xtalID, "cootOut", "Refine_" + str(Serial)
                )
            )
            Logfile.insert(
                "running giant.make_restraints %s:" % self.xtalID
                + "-ensemble-model.pdb"
            )
            cmd = (
                'export XChemExplorer_DIR="%s"\n' % os.getenv("XChemExplorer_DIR")
                + "module load buster/20240123\n"
                + "module load ccp4/7.1.018\n"
                + "giant.make_restraints %s-ensemble-model.pdb" % self.xtalID
            )
            Logfile.insert(cmd + "\n")
            os.system(cmd)
        elif os.path.isfile(
            os.path.join(self.ProjectPath, self.xtalID, "refine.split.ground-state.pdb")
        ):
            Logfile.insert(
                "found model of ground state: "
                + os.path.join(
                    self.ProjectPath, self.xtalID, "refine.split.ground-state.pdb"
                )
            )
            if os.path.isfile(
                os.path.join(
                    self.ProjectPath,
                    self.xtalID,
                    "cootOut",
                    "Refine_" + str(Serial),
                    "refine.modified.pdb",
                )
            ):
                Logfile.insert("found model of modified bound state")
                os.chdir(
                    os.path.join(
                        self.ProjectPath,
                        self.xtalID,
                        "cootOut",
                        "Refine_" + str(Serial),
                    )
                )
                ground_state = os.path.join(
                    self.ProjectPath, self.xtalID, "refine.split.ground-state.pdb"
                )
                bound_state = "refine.modified.pdb"
                Logfile.insert(
                    "running giant.merge_conformations major=%s minor=%s"
                    % (ground_state, bound_state)
                )
                if os.getcwd().startswith("/dls"):
                    cmd = (
                        'export XChemExplorer_DIR="%s"\n'
                        % os.getenv("XChemExplorer_DIR")
                        + "module load ccp4/7.1.018\n"
                        "giant.merge_conformations major=%s minor=%s"
                        " reset_all_occupancies=False options.major_occupancy=1.0"
                        " options.minor_occupancy=1.0" % (ground_state, bound_state)
                    )
                else:
                    cmd = (
                        "giant.merge_conformations major=%s minor=%s"
                        " reset_all_occupancies=False options.major_occupancy=1.0"
                        " options.minor_occupancy=1.0" % (ground_state, bound_state)
                    )
                Logfile.insert(cmd + "\n")
                os.system(cmd)
                link_line = ""
                for line in open(bound_state):
                    if line.startswith("LINK"):
                        Logfile.insert("found LINK: " + line[:-1])
                        link_line += line
                if link_line != "":
                    if os.path.isfile("multi-state-model.pdb"):
                        for line in open("multi-state-model.pdb"):
                            if line.startswith("CRYST"):
                                Logfile.insert(
                                    "adding LINK lines to multi-state-model.pdb"
                                )
                                out = link_line + line
                            else:
                                out += line
                        Logfile.insert("writing updated multi-state-model.pdb")
                        f = open("multi-state-model.pdb", "w")
                        f.write(out)
                        f.close()
                    else:
                        Logfile.error("cannot find multi-state-model.pdb")
            else:
                Logfile.error(
                    "cannot find modified version of bound state in %s"
                    % os.path.join(
                        self.ProjectPath,
                        self.xtalID,
                        "cootOut",
                        "Refine_" + str(Serial),
                    )
                )
                return None
        else:
            if os.path.isfile(
                os.path.join(
                    self.ProjectPath,
                    self.xtalID,
                    "cootOut",
                    "Refine_" + str(Serial),
                    "refine.modified.pdb",
                )
            ):
                os.chdir(
                    os.path.join(
                        self.ProjectPath,
                        self.xtalID,
                        "cootOut",
                        "Refine_" + str(Serial),
                    )
                )
                Logfile.warning(
                    "%s: ground-state model does not exist,"
                    " but refine.modified.pdb does exist" % self.xtalID
                )
                Logfile.warning(
                    "This may be a case where there are no differences"
                    " between bound and ground state"
                )
                Logfile.warning(
                    "creating symbolic link:"
                    " ln -s refine.modified.pdb %s-ensemble-model.pdb" % self.xtalID
                )
                os.system(
                    "ln -s refine.modified.pdb %s-ensemble-model.pdb" % self.xtalID
                )
                # note: after first cycle of refinement, REFMAC will remove alternate
                # conformer from the ligand
                # i.e. need to link refine.pdb to refine.split.bound-state.pdb
                make_all_links = False
                add_links_line = "ln -s refine.pdb refine.split.bound-state.pdb"
                Logfile.insert("trying to continue with refinement")
            else:
                Logfile.error(
                    "cannot find any suitable PDB file for refinement, aborting..."
                )
                return None

        #######################################################
        # checking if input PDB files are present
        Logfile.insert("checking if input PDB files for REFMAC are present")
        if os.path.isfile(
            os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                self.xtalID + "-ensemble-model.pdb",
            )
        ):
            RefmacParams["XYZIN"] = os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                self.xtalID + "-ensemble-model.pdb",
            )
        elif os.path.isfile(
            os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                "multi-state-model.pdb",
            )
        ):
            RefmacParams["XYZIN"] = os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                "multi-state-model.pdb",
            )
        else:
            Logfile.error(
                "cannot find multi-state-model.pdb in %s; aborting..."
                % os.path.join(
                    self.ProjectPath, self.xtalID, "cootOut", "Refine_" + str(Serial)
                )
            )
            return None

        #######################################################
        # checking if multi-state-restraints.refmac.params file is present
        if os.path.isfile(
            os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                "multi-state-restraints.refmac.params",
            )
        ):
            # add REFMAC keywords to multi-state-restraints.refmac.params
            with open("multi-state-restraints.refmac.params", "a") as refmacParams:
                refmacParams.write(RefmacParams["BREF"] + "\n")
                refmacParams.write(RefmacParams["TLS"] + "\n")
                refmacParams.write(RefmacParams["TWIN"] + "\n")
                refmacParams.write("ncyc " + RefmacParams["NCYCLES"] + "\n")
                if str(RefmacParams["MATRIX_WEIGHT"]).lower() == "auto":
                    refmacParams.write("weight AUTO" + "\n")
                else:
                    refmacParams.write(
                        "weight matrix " + str(RefmacParams["MATRIX_WEIGHT"]) + "\n"
                    )
                refmacParams.write(RefmacParams["TLSADD"] + "\n")
        else:
            Logfile.warning(
                "cannot find multi-state-restraints.refmac.params in %s!"
                % os.path.join(
                    self.ProjectPath, self.xtalID, "cootOut", "Refine_" + str(Serial)
                )
            )
            try:
                os.chdir(
                    os.path.join(
                        self.ProjectPath,
                        self.xtalID,
                        "cootOut",
                        "Refine_" + str(Serial),
                    )
                )
                Logfile.insert(
                    "checking if %s-ensemble-model.pdb contains residue of type"
                    " LIG, DRG, FRG, UNK or UNL" % self.xtalID
                )
                knowLigandIDs = ["LIG", "DRG", "FRG", "UNK", "UNL"]
                ligandsInFile = XChemUtils.pdbtools(
                    self.xtalID + "-ensemble-model.pdb"
                ).find_ligands()
                found_lig = False
                for lig in ligandsInFile:
                    if lig[0] in knowLigandIDs:
                        Logfile.insert("found ligand of type: " + lig[0])
                        found_lig = True
                if found_lig:
                    Logfile.warning(
                        "giant.make_restraints was not able to create"
                        " multi-state-restraints.refmac.params."
                        " Something may have gone wrong, but it could be that ligand"
                        " binding did not lead to displacement of water molecules or"
                        " rearrangement of protein side-chains."
                        " Hence, there is no difference between the bound-state and"
                        " the ground-state."
                        " We will create an empty multi-state-restraints.refmac.params"
                        " which may contain additional REFMAC keywords and otherwise"
                        " try to continue with refinement"
                    )
                    os.system("touch multi-state-restraints.refmac.params")
                else:
                    Logfile.error(
                        "%s-ensemble-model.pdb does not contain any modelled ligand;"
                        " aborting refinement" % self.xtalID
                    )
                    return None
            except OSError:
                Logfile.error(
                    "directory does not exisit: %s"
                    % os.path.join(
                        self.ProjectPath,
                        self.xtalID,
                        "cootOut",
                        "Refine_" + str(Serial),
                    )
                )
                Logfile.error("aborting refinement...")
                return None

        #######################################################
        # we write 'REFINEMENT_IN_PROGRESS' immediately to avoid unncessary refinement
        os.chdir(os.path.join(self.ProjectPath, self.xtalID))
        os.system("touch REFINEMENT_IN_PROGRESS")

        #######################################################
        # clean up!
        # and remove all files which will be re-created by current refinement cycle
        os.chdir(os.path.join(self.ProjectPath, self.xtalID))
        files_to_remove = (
            "refine.pdb "
            "refine.mtz "
            "refine.split.bound-state.pdb "
            "refine.split.ground-state.pdb "
            "validation_summary.txt "
            "validate_ligands.txt "
            "2fofc.map "
            "fofc.map "
            "refine_molprobity.log"
        )
        os.system("/bin/rm %s" % files_to_remove)

        #######################################################
        # PANDDA validation @ spider plot
        spider_plot = ""
        if os.path.isfile(
            os.path.join(
                self.ProjectPath, self.xtalID, self.xtalID + "-ensemble-model.pdb"
            )
        ):
            if os.path.isfile(
                os.path.join(
                    self.ProjectPath, self.xtalID, self.xtalID + "-pandda-input.mtz"
                )
            ):
                pdb_two = os.path.join(
                    self.ProjectPath, self.xtalID, self.xtalID + "-ensemble-model.pdb"
                )
                mtz_two = os.path.join(
                    self.ProjectPath, self.xtalID, self.xtalID + "-pandda-input.mtz"
                )
                pdb_one = os.path.join(
                    self.ProjectPath,
                    self.xtalID,
                    "Refine_" + str(panddaSerial),
                    "refine_" + str(Serial) + ".pdb",
                )
                mtz_one = os.path.join(
                    self.ProjectPath,
                    self.xtalID,
                    "Refine_" + str(panddaSerial),
                    "refine_" + str(Serial) + ".mtz",
                )
                spider_plot += (
                    "giant.score_model pdb1=%s mtz1=%s pdb2=%s mtz2=%s"
                    " res_names=LIG,UNL,DRG,FRG\n"
                    % (pdb_one, mtz_one, pdb_two, mtz_two)
                )

        #######################################################
        # CCP4 & PHENIX stuff (if working at DLS)
        module_load = ""
        if os.getcwd().startswith("/dls"):
            module_load = "module load ccp4/7.1.018\n"
            module_load += "module load phenix/1.20\n"

        # 2017-07-20: for the time being this will explicitly source pandda since
        # version 0.2 really only works at DLS
        # 2020-11-02: will start using default ccp4 installation at DLS otherwise gemmi
        # does not work
        source = ""
        if "bash" in os.getenv("SHELL"):
            source = (
                'export XChemExplorer_DIR="' + os.getenv("XChemExplorer_DIR") + '"\n'
                "\n"
            )

        if refinementProtocol == "pandda_refmac":
            refinementProgram = "refmac"
            refinementParams = os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                "multi-state-restraints.refmac.params",
            )
            mapCalculation = (
                "fft hklin refine.mtz mapout 2fofc.map << EOF\n"
                "labin F1=FWT PHI=PHWT\n"
                "EOF\n"
                "\n"
                "fft hklin refine.mtz mapout fofc.map << EOF\n"
                "labin F1=DELFWT PHI=PHDELWT\n"
                "EOF\n"
            )
        elif refinementProtocol == "pandda_phenix":
            if os.getcwd().startswith("/dls"):
                module_load = "module load phenix/1.20\n"

            refinementProgram = "phenix"
            refinementParams = os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                "multi-state-restraints.phenix.params",
            )
            mapCalculation = (
                "fft hklin refine.mtz mapout 2fofc.map << EOF\n"
                "labin F1=2FOFCWT PHI=PH2FOFCWT\n"
                "EOF\n"
                "\n"
                "fft hklin refine.mtz mapout fofc.map << EOF\n"
                "labin F1=FOFCWT PHI=PHFOFCWT\n"
                "EOF\n"
            )

        if make_all_links:
            add_links_line = (
                "ln -s Refine_%s/refine_%s.split.bound-state.pdb"
                " ./refine.split.bound-state.pdb\n"
                % (
                    panddaSerial,
                    Serial,
                )
                + "ln -s Refine_%s/refine_%s.split.ground-state.pdb"
                " ./refine.split.ground-state.pdb\n"
                % (
                    panddaSerial,
                    Serial,
                )
            )

        date = datetime.strftime(datetime.now(), "%Y-%m-%d_%H-%M-%S.%f")[:-4]
        user = getpass.getuser()

        refmacCmds = (
            "#!"
            + os.getenv("SHELL")
            + "\n"
            + module_load
            + "\n"
            + source
            + "cd "
            + self.ProjectPath
            + "/"
            + self.xtalID
            + "\n"
            "\n"
            "$CCP4/bin/ccp4-python $XChemExplorer_DIR/xce/helpers/update_status_flag.py"
            " %s %s %s %s\n"
            % (self.datasource, self.xtalID, "RefinementStatus", "running")
            + "\n"
            "giant.quick_refine"
            " input.pdb=%s" % RefmacParams["XYZIN"]
            + " mtz=%s" % RefmacParams["HKLIN"]
            + " cif=%s" % RefmacParams["LIBIN"]
            + " program=%s" % refinementProgram
            + " params=%s" % refinementParams
            + " dir_prefix='Refine_'"
            " out_prefix='refine_%s'" % str(Serial) + " split_conformations='False'"
            "\n"
            "cd "
            + self.ProjectPath
            + "/"
            + self.xtalID
            + "/Refine_"
            + str(panddaSerial)
            + "\n"
            "ln -s "
            + self.ProjectPath
            + "/"
            + self.xtalID
            + "/Refine_"
            + str(panddaSerial)
            + "/refine_"
            + str(Serial)
            + "_001.pdb "
            + self.ProjectPath
            + "/"
            + self.xtalID
            + "/Refine_"
            + str(panddaSerial)
            + "/refine_"
            + str(Serial)
            + ".pdb"
            + "\n"
            "ln -s "
            + self.ProjectPath
            + "/"
            + self.xtalID
            + "/Refine_"
            + str(panddaSerial)
            + "/refine_"
            + str(Serial)
            + "_001.mtz "
            + self.ProjectPath
            + "/"
            + self.xtalID
            + "/Refine_"
            + str(panddaSerial)
            + "/refine_"
            + str(Serial)
            + ".mtz"
            + "\n"
            "giant.split_conformations"
            " input.pdb='refine_%s.pdb'" % str(Serial) + " reset_occupancies=False"
            " suffix_prefix=split"
            "\n"
            "giant.split_conformations"
            " input.pdb='refine_%s.pdb'" % str(Serial) + " reset_occupancies=True"
            " suffix_prefix=output "
            "\n" + spider_plot + "\n"
            "phenix.molprobity refine_%s.pdb refine_%s.mtz\n" % (Serial, Serial)
            + "/bin/mv molprobity.out refine_molprobity.log\n"
            "module load phenix/1.20\n"
            "mmtbx.validate_ligands refine_%s.pdb refine_%s.mtz LIG"
            " > validate_ligands.txt\n" % (Serial, Serial)
            + "cd "
            + self.ProjectPath
            + "/"
            + self.xtalID
            + "\n"
            "\n"
            "ln -s Refine_%s/validate_ligands.txt .\n" % panddaSerial
            + "ln -s Refine_%s/refine_molprobity.log .\n" % panddaSerial
            + "\n"
            + add_links_line
            + "\n"
            "mmtbx.validation_summary refine.pdb > validation_summary.txt\n"
            "\n" + mapCalculation + "\n"
            "$CCP4/bin/ccp4-python "
            + os.path.join(
                os.getenv("XChemExplorer_DIR"),
                "xce",
                "helpers",
                "update_data_source_after_refinement.py",
            )
            + " %s %s %s %s %s %s\n"
            % (
                self.datasource,
                self.xtalID,
                self.ProjectPath,
                os.path.join(
                    self.ProjectPath, self.xtalID, "Refine_" + str(panddaSerial)
                ),
                user,
                date,
            )
            + "\n"
            "/bin/rm %s/%s/REFINEMENT_IN_PROGRESS\n" % (self.ProjectPath, self.xtalID)
            + "\n"
        )

        Logfile.insert(
            "writing refinement shell script to"
            + os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                "refmac.csh",
            )
        )
        cmd = open(
            os.path.join(
                self.ProjectPath,
                self.xtalID,
                "cootOut",
                "Refine_" + str(Serial),
                "refmac.csh",
            ),
            "w",
        )
        cmd.write(refmacCmds)
        cmd.close()

        os.chdir(
            os.path.join(
                self.ProjectPath, self.xtalID, "cootOut", "Refine_" + str(Serial)
            )
        )
        Logfile.insert(
            "changing directory to %s"
            % (
                os.path.join(
                    self.ProjectPath, self.xtalID, "cootOut", "Refine_" + str(Serial)
                )
            )
        )

        if external_software["slurm"]:
            slurm.submit_cluster_job("refmac", "refmac.csh", xce_logfile)
        elif external_software["qsub"]:
            sge.submit_cluster_job("refmac", "refmac.csh", xce_logfile)
        else:
            Logfile.insert("changing permission of refmac.csh: chmod +x refmac.csh")
            os.system("chmod +x refmac.csh")
            Logfile.insert("starting refinement on local machine")
            os.system("./refmac.csh &")

    def RefinementParams(self, RefmacParams):
        self.RefmacParams = RefmacParams
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", gtk.main_quit)
        self.window.set_border_width(10)
        self.window.set_title("Refmac Parameters")
        self.vbox = gtk.VBox()

        self.hbox1 = gtk.HBox()
        self.hbox1.add(gtk.Label("Refine"))
        self.cb = gtk.combo_box_new_text()
        self.cb.connect("changed", self.ChooseBfacRefinement)
        for item in ["isotropic", "anisotropic"]:
            self.cb.append_text(item)
        if "ISOT" in self.RefmacParams["BREF"]:
            self.cb.set_active(0)
        if "ANIS" in self.RefmacParams["BREF"]:
            self.cb.set_active(1)
        self.hbox1.add(self.cb)
        self.hbox1.add(gtk.Label("temperature factors"))
        self.vbox.add(self.hbox1)

        self.hbox2 = gtk.HBox()
        self.hbox2.add(gtk.Label("Number of Cycles: "))
        self.Ncycles = gtk.Entry()
        self.Ncycles.add_events(gtk.gdk.KEY_RELEASE_MASK)
        self.Ncycles.connect("key-release-event", self.on_key_release_Ncycles)
        self.Ncycles.set_text(self.RefmacParams["NCYCLES"])
        self.hbox2.add(self.Ncycles)
        self.vbox.add(self.hbox2)

        self.hbox3 = gtk.HBox()
        self.hbox3.add(gtk.Label("MATRIX WEIGHT: "))
        self.MATRIX_WEIGHT = gtk.Entry()
        self.MATRIX_WEIGHT.add_events(gtk.gdk.KEY_RELEASE_MASK)
        self.MATRIX_WEIGHT.connect(
            "key-release-event", self.on_key_release_MATRIX_WEIGHT
        )
        self.MATRIX_WEIGHT.set_text(self.RefmacParams["MATRIX_WEIGHT"])
        self.hbox3.add(self.MATRIX_WEIGHT)
        self.vbox.add(self.hbox3)

        self.TLS = gtk.CheckButton("TLS (find TLS groups with phenix.find_tls_groups)")
        self.TLS.connect("toggled", self.TLSCallback)
        if self.RefmacParams["TLS"] == "refi tlsc 10\n":
            self.TLS.set_active(True)
        self.vbox.pack_start(self.TLS, False)

        self.NCS = gtk.CheckButton("NCS (if applicable")
        self.NCS.connect("toggled", self.NCSCallback)
        if self.RefmacParams["NCS"] == "NCSR LOCAL\n":
            self.NCS.set_active(True)
        self.vbox.pack_start(self.NCS, False)

        self.TWIN = gtk.CheckButton("Twin?")
        self.TWIN.connect("toggled", self.TWINCallback)
        if self.RefmacParams["TWIN"] == "TWIN\n":
            self.TWIN.set_active(True)
        self.vbox.pack_start(self.TWIN, False)

        self.OKbutton = gtk.Button(label="OK")
        self.OKbutton.connect("clicked", self.OK)
        self.vbox.add(self.OKbutton)

        self.window.add(self.vbox)
        self.window.show_all()
        return self.RefmacParams

    def TLSCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["TLS"] = "refi tlsc 10\n"
            self.RefmacParams["TLSIN"] = "refmac.tls\n"
            self.RefmacParams["TLSOUT"] = "out.tls\n"
            self.RefmacParams["TLSADD"] = "TLSO ADDU\n"
        else:
            self.RefmacParams["TLS"] = ""
            self.RefmacParams["TLSIN"] = ""
            self.RefmacParams["TLSOUT"] = ""
            self.RefmacParams["TLSADD"] = ""
        return self.RefmacParams

    def NCSCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["NCS"] = "NCSR LOCAL\n"
        else:
            self.RefmacParams["NCS"] = ""
        return self.RefmacParams

    def ChooseBfacRefinement(self, widget):
        if widget.get_active_text() == "isotropic":
            self.RefmacParams["BREF"] = "    bref ISOT\n"
        if widget.get_active_text() == "anisotropic":
            self.RefmacParams["BREF"] = "    bref ANIS\n"
        return self.RefmacParams

    def on_key_release_Ncycles(self, widget, event):
        print(widget.get_text())
        self.RefmacParams["NCYCLES"] = widget.get_text()
        return self.RefmacParams

    def on_key_release_MATRIX_WEIGHT(self, widget, event):
        self.RefmacParams["MATRIX_WEIGHT"] = widget.get_text()
        return self.RefmacParams

    def TWINCallback(self, widget):
        if widget.get_active():
            self.RefmacParams["TWIN"] = "TWIN\n"
        else:
            self.RefmacParams["TWIN"] = ""
        return self.RefmacParams

    def OK(self, widget):
        self.window.destroy()

    def ParamsFromPreviousCycle(self, Serial):
        RefmacParams = {
            "HKLIN": "",
            "HKLOUT": "",
            "XYZIN": "",
            "XYZOUT": "",
            "LIBIN": "",
            "LIBOUT": "",
            "TLSIN": "",
            "TLSOUT": "",
            "TLSADD": "",
            "NCYCLES": "10",
            "MATRIX_WEIGHT": "AUTO",
            "BREF": "    bref ISOT\n",
            "TLS": "",
            "NCS": "",
            "TWIN": "",
        }

        if os.path.isfile(
            self.ProjectPath
            + "/"
            + self.xtalID
            + "/Refine_"
            + str(Serial)
            + "/refmac.csh"
        ):
            for line in open(
                self.ProjectPath
                + "/"
                + self.xtalID
                + "/Refine_"
                + str(Serial)
                + "/refmac.csh"
            ):
                if line.startswith("refi tlsc"):
                    RefmacParams["TLS"] = line
                if line.startswith("TLSO"):
                    RefmacParams["TLSADD"] = line
                if line.startswith("NCSR LOCAL"):
                    RefmacParams["NCS"] = line
                if line.startswith("    bref "):
                    RefmacParams["BREF"] = line
                if line.startswith("ncyc"):
                    RefmacParams["Ncycles"] = line.split()[1]
                if line.startswith("weight"):
                    RefmacParams["MATRIX_WEIGHT"] = line.split()[len(line.split()) - 1]
                if line.startswith("TWIN"):
                    RefmacParams["TWIN"] = line

        return RefmacParams

    def GetRefinementHistory(self):
        RefinementCycle = []
        RcrystList = []
        RfreeList = []

        found = False
        for item in glob.glob(os.path.join(self.ProjectPath, self.xtalID, "*")):
            if item.startswith(os.path.join(self.ProjectPath, self.xtalID, "Refine_")):
                print(item[item.rfind("_") + 1 :])
                RefinementCycle.append(int(item[item.rfind("_") + 1 :]))
                found = True
        if found:
            for cycle in sorted(RefinementCycle):
                try:
                    found_Rcryst = False
                    found_Rfree = False
                    newestPDB = max(
                        glob.iglob(
                            self.ProjectPath
                            + "/"
                            + self.xtalID
                            + "/Refine_"
                            + str(cycle)
                            + "/refine_"
                            + str(cycle)
                            + ".pdb"
                        ),
                        key=os.path.getctime,
                    )
                    for line in open(newestPDB):
                        if line.startswith(
                            "REMARK   3   R VALUE     (WORKING + TEST SET) :"
                        ):
                            Rcryst = line.split()[9]
                            RcrystList.append(float(Rcryst))
                            found_Rcryst = True
                        if line.startswith(
                            "REMARK   3   FREE R VALUE                     :"
                        ):
                            Rfree = line.split()[6]
                            RfreeList.append(float(Rfree))
                            found_Rfree = True
                    if not found_Rcryst:
                        RcrystList.append(0)
                    if not found_Rfree:
                        RfreeList.append(0)
                except ValueError:
                    RcrystList.append(0)
                    RfreeList.append(0)
        else:
            RefinementCycle = [0]
            RcrystList = [0]
            RfreeList = [0]
        print(RefinementCycle, RcrystList, RfreeList)
        return (sorted(RefinementCycle), RcrystList, RfreeList)
