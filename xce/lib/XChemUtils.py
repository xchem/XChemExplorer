import bz2
import gzip
import json
import math
import os
import subprocess
import time
from distutils.spawn import find_executable

from rdkit import Chem

import gemmi
import iotbx.pdb
from xce.lib import XChemDB
from xce.lib import XChemLog
from iotbx import mtz
from iotbx.reflection_file_reader import any_reflection_file


def open_decompress_file(filename, mode="rb"):
    """Open file for reading, decompressing silently if necessary."""

    if filename.endswith(".bz2"):
        return bz2.BZ2File(filename, mode=mode)

    elif filename.endswith(".gz"):
        return gzip.GzipFile(filename, mode=mode)

    return open(filename, mode=mode)


class helpers:
    def make_png(
        self,
        initial_model_directory,
        sample,
        compoundID,
        smiles,
        external_software,
        database_directory,
        data_source_file,
        ccp4_scratch_directory,
        counter,
        xce_logfile,
        restraints_program,
    ):
        Logfile = XChemLog.updateLog(xce_logfile)

        os.system("touch RESTRAINTS_IN_PROGRESS")

        header = "#!" + os.getenv("SHELL") + "\n"
        if external_software["qsub"]:
            if not external_software["qsub_array"]:
                header = "#PBS -joe -N xce_acedrg\n"

        # check for combisoaks/ cocktails
        # note: compoundIDs and smiles are semi-colon separated
        if len(compoundID.split(";")) > 1:
            Logfile.insert(
                "looks like you are working with cocktails;"
                " found the following IDs and smiles:"
            )
            if len(compoundID.split(";")) != len(smiles.split(";")):
                Logfile.error("Number of compoundIDs and SMILES does not match:")
                Logfile.error(
                    "N(compoundID): {0!s} -> {1!s}".format(
                        len(compoundID.split(";")), compoundID.split(";")
                    )
                )
                Logfile.error(
                    "N(SMILES):     {0!s} -> {1!s}".format(
                        len(smiles.split(";")), smiles.split(";")
                    )
                )
                Logfile.error("aborting...")
                return
            else:
                software = ""
                if restraints_program == "acedrg":
                    if os.getcwd().startswith("/dls"):
                        software += "module load ccp4/7.1.018\n"
                elif restraints_program == "phenix.elbow":
                    if os.getcwd().startswith("/dls"):
                        software += "module load ccp4/7.1.018\n"
                        software += "module load phenix/1.20\n"
                elif restraints_program == "grade":
                    if os.getcwd().startswith("/dls"):
                        software += "module load ccp4/7.1.018\n"
                        software += "module load buster/20211020\n"
                        software += (
                            "export BDG_TOOL_MOGUL="
                            "/dls_sw/apps/ccdc/CSD_2020/bin/mogul\n"
                        )
                    software += "export BDG_TOOL_OBABEL='none'\n"

                for i in range(len(compoundID.split(";"))):
                    Logfile.insert(
                        "{0!s} - {1!s}".format(
                            compoundID.split(";")[i], smiles.split(";")[i]
                        )
                    )
                    cID = "L" + (2 - len(str(i))) * "0" + str(i)
                    if restraints_program == "acedrg":
                        software += 'acedrg --res {0!s} -i "{1!s}" -o {2!s}\n'.format(
                            cID,
                            smiles.split(";")[i],
                            compoundID.split(";")[i].replace(" ", ""),
                        )
                    elif restraints_program == "phenix.elbow":
                        software += (
                            "phenix.elbow"
                            + '--smiles="{0!s}" --id {1!s} --output {2!s}\n'.format(
                                smiles.split(";")[i],
                                cID,
                                compoundID.split(";")[i].replace(" ", ""),
                            )
                        )
                    elif restraints_program == "grade":
                        if external_software["mogul"]:
                            mogul = ""
                        else:
                            mogul = "-nomogul"
                        software += (
                            "grade"
                            + ' -resname {0!s} {1!s} "{2!s}"'.format(
                                cID, mogul, smiles.split(";")[i]
                            )
                            + " -ocif {0!s}.cif -opdb {1!s}.pdb\n".format(
                                compoundID.split(";")[i].replace(" ", ""),
                                compoundID.split(";")[i].replace(" ", ""),
                            )
                        )

        else:
            # check if CompoundSMILEScovalent field is not Null
            # CompoundSMILESproduct can be used to create only a CIF file
            # for the product to make fitting easier
            # however, the complete smiles string will be used to make the png file
            productSmiles = None
            db = XChemDB.data_source(os.path.join(database_directory, data_source_file))
            sql = (
                "select CompoundSMILESproduct from mainTable where CrystalName = '%s'"
                % sample
            )
            query = db.execute_statement(sql)
            productSmiles = query[0][0]
            if str(productSmiles).replace(" ", "") == "":
                productSmiles = smiles
            elif "none" in str(productSmiles).lower():
                productSmiles = smiles
            elif "null" in str(productSmiles).lower():
                productSmiles = smiles
            else:
                productSmiles = str(productSmiles)

            software = ""
            if restraints_program == "acedrg":
                if os.getcwd().startswith("/dls"):
                    software += "module load ccp4/7.1.018\n"
                if os.path.isfile(
                    os.path.join(initial_model_directory, sample, "old.cif")
                ):
                    software += "acedrg --res LIG -c ../old.cif -o {0!s}\n".format(
                        (compoundID.replace(" ", ""))
                    )
                else:
                    software += 'acedrg --res LIG -i "{0!s}" -o {1!s}\n'.format(
                        productSmiles, compoundID.replace(" ", "")
                    )
            elif restraints_program == "phenix.elbow":
                if os.getcwd().startswith("/dls"):
                    software += "module load ccp4/7.1.018\n"
                    software += "module load phenix/1.20\n"
                if os.path.isfile(
                    os.path.join(initial_model_directory, sample, "old.cif")
                ):
                    software += (
                        "phenix.elbow --file=../old.cif --id LIG"
                        + " --output {0!s}\n".format((compoundID.replace(" ", "")))
                    )
                else:
                    software += (
                        "phenix.elbow"
                        + ' --smiles="{0!s}" --id LIG --output {1!s}\n'.format(
                            productSmiles, compoundID.replace(" ", "")
                        )
                    )
            elif restraints_program == "grade":
                if os.getcwd().startswith("/dls"):
                    software += "module load ccp4/7.1.018\n"
                    software += "module load buster/20211020\n"
                    software += (
                        "export BDG_TOOL_MOGUL=/dls_sw/apps/ccdc/CSD_2020/bin/mogul\n"
                    )
                software += "export BDG_TOOL_OBABEL='none'\n"
                if external_software["mogul"]:
                    mogul = ""
                else:
                    mogul = "-nomogul"
                if os.path.isfile(
                    os.path.join(initial_model_directory, sample, "old.cif")
                ):
                    software += "grade -resname LIG {0!s}".format(
                        mogul
                    ) + " -in ../old.cif -ocif {0!s}.cif -opdb {1!s}.pdb\n".format(
                        compoundID.replace(" ", ""), compoundID.replace(" ", "")
                    )
                else:
                    software += 'grade -resname LIG {0!s} "{1!s}"'.format(
                        mogul, productSmiles
                    ) + " -ocif {0!s}.cif -opdb {1!s}.pdb\n".format(
                        compoundID.replace(" ", ""),
                        compoundID.replace(" ", ""),
                    )

            # merge all compound CIFs into 1 file called merged.cif
            software += (
                "$CCP4/bin/ccp4-python"
                " $XChemExplorer_DIR/xce/helpers/"
                "merge_ligand_cif_files.py {0!s}\n".format(
                    os.path.join(initial_model_directory, sample, "compound")
                )
            )

        # Removal of the hydrogen atoms in PDB files is required for REFMAC 5 run.
        # With hydrogens some ligands fail to pass the external restraints in
        # pandda.giant.make_restraints.
        # Copy the file with hydrogens to retain in case needed
        Cmds = (
            header + "\n"
            'export XChemExplorer_DIR="' + os.getenv("XChemExplorer_DIR") + '"' + "\n"
            "$CCP4/bin/ccp4-python $XChemExplorer_DIR/xce/helpers/update_status_flag.py"
            " {0!s} {1!s} {2!s} {3!s}".format(
                os.path.join(database_directory, data_source_file),
                sample,
                "RefinementCIFStatus",
                "running",
            )
            + "\n"
            '$CCP4/bin/ccp4-python {0!s} "{1!s}" {2!s} {3!s} {4!s}\n'.format(
                os.path.join(
                    os.getenv("XChemExplorer_DIR"),
                    "xce",
                    "helpers",
                    "create_png_of_compound.py",
                ),
                smiles,
                compoundID.replace(" ", ""),
                sample,
                initial_model_directory,
            )
            + "\n"
            "cd "
            + os.path.join(initial_model_directory, sample, "compound")
            + "\n"
            + software
            + "\n"
            "cd " + os.path.join(initial_model_directory, sample) + "\n"
            "ln -s compound/%s.cif .\n" % compoundID.replace(" ", "")
            + "ln -s compound/{0!s}.pdb .\n".format(compoundID.replace(" ", ""))
            + "ln -s compound/{0!s}.png .\n".format(compoundID.replace(" ", ""))
            + "\n"
            "$CCP4/bin/ccp4-python "
            + os.path.join(
                os.getenv("XChemExplorer_DIR"),
                "xce",
                "helpers",
                "update_data_source_for_new_cif_files.py",
            )
            + " {0!s} {1!s} {2!s} {3!s}\n".format(
                os.path.join(database_directory, data_source_file),
                sample,
                initial_model_directory,
                compoundID.replace(" ", ""),
            )
            + "\n"
            "/bin/rm -f compound/RESTRAINTS_IN_PROGRESS\n"
        )

        os.chdir(ccp4_scratch_directory)
        Logfile.insert(
            "creating ACEDRG shell script for {0!s},{1!s} in {2!s}".format(
                sample, compoundID, ccp4_scratch_directory
            )
        )
        f = open("xce_{0!s}_{1!s}.sh".format(restraints_program, str(counter)), "w")
        f.write(Cmds)
        f.close()
        os.system(
            "chmod +x xce_{0!s}_{1!s}.sh".format(restraints_program, str(counter))
        )


class parse:
    def __init__(self):
        self.space_group_dict = {
            "triclinic": ["P1"],
            "monoclinic_P": ["P2", "P21", "P1211", "P121"],
            "monoclinic_C": ["C2", "C121"],
            "orthorhombic": [
                "P222",
                "P2122",
                "P2212",
                "P2221",
                "P21212",
                "P21221",
                "P22121",
                "P212121",
                "C222",
                "C2221",
                "F222",
                "I222",
                "I212121",
            ],
            "tetragonal": [
                "P4",
                "P41",
                "P42",
                "P43",
                "I4",
                "I41",
                "P422",
                "P4212",
                "P4122",
                "P41212",
                "P4222",
                "P42212",
                "P4322",
                "P43212",
                "I422",
                "I4122",
            ],
            "hexagonal": [
                "P3",
                "P31",
                "P32",
                "P312",
                "P321",
                "P3112",
                "P3121",
                "P3212",
                "P3221",
                "P6",
                "P61",
                "P65",
                "P62",
                "P64",
                "P63",
                "P622",
                "P6122",
                "P6522",
                "P6222",
                "P6422",
                "P6322",
            ],
            "rhombohedral": ["H3", "H32"],
            "cubic": [
                "P23",
                "F23",
                "I23",
                "P213",
                "I213",
                "P432",
                "P4232",
                "F432",
                "F4132",
                "I432",
                "P4332",
                "P4132",
                "I4132",
            ],
        }

        self.point_group_dict = {
            "1": ["P1"],
            "2": ["P2", "P21", "C121", "P1211", "P121", "C2"],
            "222": [
                "P222",
                "P2122",
                "P2212",
                "P2221",
                "P21212",
                "P21221",
                "P22121",
                "P212121",
                "C222",
                "C2221",
                "F222",
                "I222",
                "I212121",
            ],
            "4": ["P4", "P41", "P42", "P43", "I4", "I41"],
            "422": [
                "P422",
                "P4212",
                "P4122",
                "P41212",
                "P4222",
                "P42212",
                "P4322",
                "P43212",
                "I422",
                "I4122",
            ],
            "3": ["P3", "P31", "P32", "H3"],
            "32": ["P312", "P321", "P3112", "P3121", "P3212", "P3221", "H32"],
            "6": ["P6", "P61", "P65", "P62", "P64", "P63"],
            "622": ["P622", "P6122", "P6522", "P6222", "P6422", "P6322"],
            "23": ["P23", "F23", "I23", "P213", "I213"],
            "432": [
                "P432",
                "P4232",
                "F432",
                "F4132",
                "I432",
                "P4332",
                "P4132",
                "I4132",
            ],
        }

        self.nr_asu_in_unitcell_for_point_group = {
            "1": 1,
            "2": 2,
            "222": 4,
            "4": 4,
            "422": 8,
            "3": 3,
            "32": 6,
            "6": 6,
            "622": 12,
            "23": 12,
            "432": 24,
        }

        self.aimless = {
            "DataProcessingProgram": "n/a",
            "DataProcessingSpaceGroup": "n/a",
            "DataProcessingUnitCell": "n/a",
            "DataProcessingA": "n/a",
            "DataProcessingB": "n/a",
            "DataProcessingC": "n/a",
            "DataProcessingAlpha": "n/a",
            "DataProcessingBeta": "n/a",
            "DataProcessingGamma": "n/a",
            "DataProcessingResolutionLow": "n/a",
            "DataProcessingResolutionLowInnerShell": "n/a",
            "DataProcessingResolutionHigh": "n/a",
            "DataProcessingResolutionHighOuterShell": "n/a",
            "DataProcessingResolutionOverall": "n/a",
            "DataProcessingRmergeOverall": "n/a",
            "DataProcessingRmergeLow": "n/a",
            "DataProcessingRmergeHigh": "n/a",
            "DataProcessingIsigOverall": "n/a",
            "DataProcessingIsigLow": "n/a",
            "DataProcessingIsigHigh": "n/a",
            "DataProcessingCompletenessOverall": "n/a",
            "DataProcessingCompletenessLow": "n/a",
            "DataProcessingCompletenessHigh": "n/a",
            "DataProcessingMultiplicityOverall": "n/a",
            "DataProcessingMultiplicityLow": "n/a",
            "DataProcessingMultiplicityHigh": "n/a",
            "DataProcessingCChalfOverall": "n/a",
            "DataProcessingCChalfLow": "n/a",
            "DataProcessingCChalfHigh": "n/a",
            "DataProcessingResolutionHigh15sigma": "n/a",
            "DataProcessingUniqueReflectionsLow": "n/a",
            "DataProcessingUniqueReflectionsHigh": "n/a",
            "DataProcessingUniqueReflectionsOverall": "n/a",
            "DataProcessingLattice": "n/a",
            "DataProcessingPointGroup": "n/a",
            "DataProcessingUnitCellVolume": 0,
            "DataProcessingAlert": "#FF0000",
            "DataCollectionWavelength": "n/a",
            "DataProcessingScore": "n/a",
        }

    def read_aimless_logfile(self, logfile):
        # essentially same as above, but compatible with datasource
        # will hopefully supersede function above

        a = "n/a"
        b = "n/a"
        c = "n/a"
        alpha = "n/a"
        beta = "n/a"
        gamma = "n/a"

        if "fast_dp" in logfile:
            self.aimless["DataProcessingProgram"] = "fast_dp"
        elif "3d-run" in logfile:
            self.aimless["DataProcessingProgram"] = "xia2-3d"
        elif "3dii-run" in logfile:
            self.aimless["DataProcessingProgram"] = "xia2-3dii"
        elif "xia2-3dii" in logfile:
            self.aimless["DataProcessingProgram"] = "xia2-3dii"
        elif "dials-run" in logfile:
            self.aimless["DataProcessingProgram"] = "dials"
        elif "dials" in logfile:
            self.aimless["DataProcessingProgram"] = "dials"
        elif "autoPROC" in logfile:
            self.aimless["DataProcessingProgram"] = "autoPROC"
        elif "staraniso" in logfile:
            self.aimless["DataProcessingProgram"] = "aP_staraniso"

        # get run number from logfile
        # only works if file is in original directory, but not once it moved to
        # 'inital_model' folder

        if logfile.endswith(".log") or logfile.endswith(".table1"):
            self.aimless_logile(logfile)
        elif logfile.endswith(".json"):
            self.json_logfile(logfile)

        if (
            self.aimless["DataProcessingA"] != "n/a"
            and self.aimless["DataProcessingB"] != "n/a"
            and self.aimless["DataProcessingC"] != "n/a"
            and self.aimless["DataProcessingAlpha"] != "n/a"
            and self.aimless["DataProcessingBeta"] != "n/a"
            and self.aimless["DataProcessingGamma"] != "n/a"
            and self.aimless["DataProcessingLattice"] != "n/a"
        ):
            a = self.aimless["DataProcessingA"]
            b = self.aimless["DataProcessingB"]
            c = self.aimless["DataProcessingC"]
            alpha = self.aimless["DataProcessingAlpha"]
            beta = self.aimless["DataProcessingBeta"]
            gamma = self.aimless["DataProcessingGamma"]
            self.aimless["DataProcessingUnitCellVolume"] = str(
                self.calc_unitcell_volume_from_logfile(
                    float(a),
                    float(b),
                    float(c),
                    math.radians(float(alpha)),
                    math.radians(float(beta)),
                    math.radians(float(gamma)),
                    self.aimless["DataProcessingLattice"],
                )
            )
            try:
                high_symmetry_boost = self.nr_asu_in_unitcell_for_point_group[
                    self.aimless["DataProcessingPointGroup"]
                ]
                self.aimless["DataProcessingScore"] = (
                    float(self.aimless["DataProcessingUniqueReflectionsOverall"])
                    * float(self.aimless["DataProcessingCompletenessOverall"])
                    * high_symmetry_boost
                    * float(self.aimless["DataProcessingIsigOverall"])
                ) / float(self.aimless["DataProcessingUnitCellVolume"])
            # When P-6 was accidentally used self.aimless['DataProcessingPointGroup']
            # through a KeyError, so handling this
            except (ValueError, KeyError):
                self.aimless["DataProcessingScore"] = 0.0
        self.aimless["DataProcessingUnitCell"] = (
            str(a)
            + " "
            + str(b)
            + " "
            + str(c)
            + " "
            + str(alpha)
            + " "
            + str(beta)
            + " "
            + str(gamma)
        )
        self.aimless["DataProcessingResolutionOverall"] = (
            str(self.aimless["DataProcessingResolutionLow"])
            + " - "
            + str(self.aimless["DataProcessingResolutionHigh"])
        )

        if (
            self.aimless["DataProcessingResolutionHigh"] == "n/a"
            or self.aimless["DataProcessingRmergeLow"] == "n/a"
        ):
            self.aimless["DataProcessingAlert"] = "#FF0000"
        else:
            if (
                float(self.aimless["DataProcessingResolutionHigh"]) > 3.5
                or float(self.aimless["DataProcessingRmergeLow"]) > 0.1
            ):
                self.aimless["DataProcessingAlert"] = "#FF0000"
            if (3.5 >= float(self.aimless["DataProcessingResolutionHigh"]) > 2.5) or (
                0.1 >= float(self.aimless["DataProcessingRmergeLow"]) > 0.05
            ):
                self.aimless["DataProcessingAlert"] = "#FF9900"
            if (
                float(self.aimless["DataProcessingResolutionHigh"]) <= 2.5
                and float(self.aimless["DataProcessingRmergeLow"]) <= 0.05
            ):
                self.aimless["DataProcessingAlert"] = "#00FF00"

        return self.aimless

    def aimless_logile(self, logfile):
        resolution_at_15_sigma_line_overall_found = False
        resolution_at_20_sigma_line_overall_found = False
        for _, line in enumerate(open(logfile)):
            if "Wavelength" in line and len(line.split()) >= 2:
                self.aimless["DataCollectionWavelength"] = line.split()[1]
            if "Low resolution limit" in line and len(line.split()) == 6:
                self.aimless["DataProcessingResolutionLow"] = line.split()[3]
                self.aimless["DataProcessingResolutionHighOuterShell"] = line.split()[5]
            if "High resolution limit" in line and len(line.split()) == 6:
                self.aimless["DataProcessingResolutionHigh"] = line.split()[3]
                self.aimless["DataProcessingResolutionLowInnerShell"] = line.split()[4]
            if "Rmerge  (all I+ and I-)" in line and len(line.split()) == 8:
                self.aimless["DataProcessingRmergeOverall"] = line.split()[5]
                self.aimless["DataProcessingRmergeLow"] = line.split()[6]
                self.aimless["DataProcessingRmergeHigh"] = line.split()[7]
            if "Rmerge  (all I+ & I-)" in line and len(line.split()) == 8:
                self.aimless["DataProcessingRmergeOverall"] = line.split()[5]
                self.aimless["DataProcessingRmergeLow"] = line.split()[6]
                self.aimless["DataProcessingRmergeHigh"] = line.split()[7]
            if "Mean((I)/sd(I))" in line and len(line.split()) == 4:
                self.aimless["DataProcessingIsigOverall"] = line.split()[1]
                self.aimless["DataProcessingIsigHigh"] = line.split()[3]
                self.aimless["DataProcessingIsigLow"] = line.split()[2]
            if "Mean(I)/sd(I)" in line and len(line.split()) == 4:
                self.aimless["DataProcessingIsigOverall"] = line.split()[1]
                self.aimless["DataProcessingIsigHigh"] = line.split()[3]
                self.aimless["DataProcessingIsigLow"] = line.split()[2]
            if line.startswith("Completeness") and len(line.split()) == 4:
                self.aimless["DataProcessingCompletenessOverall"] = line.split()[1]
                self.aimless["DataProcessingCompletenessHigh"] = line.split()[3]
                self.aimless["DataProcessingCompletenessLow"] = line.split()[2]
            if "Completeness (ellipsoidal)" in line and len(line.split()) == 5:
                self.aimless["DataProcessingCompletenessOverall"] = line.split()[2]
                self.aimless["DataProcessingCompletenessHigh"] = line.split()[4]
                self.aimless["DataProcessingCompletenessLow"] = line.split()[3]
            if "Multiplicity" in line and len(line.split()) == 4:
                self.aimless["DataProcessingMultiplicityOverall"] = line.split()[1]
                self.aimless["DataProcessingMultiplicityHigh"] = line.split()[3]
                self.aimless["DataProcessingMultiplicityLow"] = line.split()[3]
            if (
                line.startswith("Mn(I) half-set correlation CC(1/2)")
                and len(line.split()) == 7
            ):
                self.aimless["DataProcessingCChalfOverall"] = line.split()[4]
                self.aimless["DataProcessingCChalfLow"] = line.split()[5]
                self.aimless["DataProcessingCChalfHigh"] = line.split()[6]
            if line.startswith("     CC(1/2)") and len(line.split()) == 4:
                self.aimless["DataProcessingCChalfOverall"] = line.split()[1]
                self.aimless["DataProcessingCChalfLow"] = line.split()[2]
                self.aimless["DataProcessingCChalfHigh"] = line.split()[3]
            if line.startswith("Estimates of resolution limits: overall"):
                resolution_at_15_sigma_line_overall_found = True
                resolution_at_20_sigma_line_overall_found = True
            if resolution_at_15_sigma_line_overall_found:
                if "from Mn(I/sd)" in line and len(line.split()) >= 7:
                    if "1.5" in line.split()[3]:
                        self.aimless[
                            "DataProcessingResolutionHigh15sigma"
                        ] = line.split()[6][:-1]
                        resolution_at_15_sigma_line_overall_found = False
            if resolution_at_20_sigma_line_overall_found:
                if "from Mn(I/sd)" in line and len(line.split()) >= 7:
                    if "2.0" in line.split()[3]:
                        self.aimless[
                            "DataProcessingResolutionHigh20sigma"
                        ] = line.split()[6][:-1]
                        resolution_at_20_sigma_line_overall_found = False
            if (
                line.startswith("Average unit cell:")
                or line.startswith("  Unit cell parameters")
            ) and len(line.split()) == 9:
                tmp = [line.split()]
                a = int(float(tmp[0][3]))
                b = int(float(tmp[0][4]))
                c = int(float(tmp[0][5]))
                alpha = int(float(tmp[0][6]))
                beta = int(float(tmp[0][7]))
                gamma = int(float(tmp[0][8]))
                self.aimless["DataProcessingA"] = str(a)
                self.aimless["DataProcessingB"] = str(b)
                self.aimless["DataProcessingC"] = str(c)
                self.aimless["DataProcessingAlpha"] = str(alpha)
                self.aimless["DataProcessingBeta"] = str(beta)
                self.aimless["DataProcessingGamma"] = str(gamma)
            if "Total number unique" in line and len(line.split()) == 6:
                self.aimless["DataProcessingUniqueReflectionsOverall"] = line.split()[3]
            if line.startswith("Space group:") or line.startswith("  Spacegroup name"):
                if "Laue" in line:
                    continue
                if "Spacegroup name" in line:
                    self.aimless["DataProcessingSpaceGroup"] = line.replace(
                        "  Spacegroup name", ""
                    )[:-1].replace(" ", "")
                else:
                    self.aimless["DataProcessingSpaceGroup"] = line.replace(
                        "Space group: ", ""
                    )[:-1]
                self.aimless[
                    "DataProcessingLattice"
                ] = self.get_lattice_from_space_group(
                    self.aimless["DataProcessingSpaceGroup"]
                )
                self.aimless[
                    "DataProcessingPointGroup"
                ] = self.get_pointgroup_from_space_group(
                    self.aimless["DataProcessingSpaceGroup"]
                )

    def json_logfile(self, logfile):
        with open(logfile, "r") as log:
            data = log.read()
        obj = json.loads(data)
        self.aimless["DataProcessingResolutionLow"] = str(
            round(math.sqrt(1 / float(obj["d_star_sq_max"][0])), 2)
        )
        self.aimless["DataProcessingResolutionLowInnerShell"] = str(
            round(math.sqrt(1 / float(obj["d_star_sq_max"][1])), 2)
        )
        self.aimless["DataProcessingResolutionHigh"] = str(
            round(
                math.sqrt(
                    1 / float(obj["d_star_sq_min"][len(obj["d_star_sq_min"]) - 1])
                ),
                2,
            )
        )
        self.aimless["DataProcessingResolutionHighOuterShell"] = str(
            round(
                math.sqrt(
                    1 / float(obj["d_star_sq_min"][len(obj["d_star_sq_min"]) - 2])
                ),
                2,
            )
        )
        self.aimless["DataProcessingResolutionOverall"] = (
            self.aimless["DataProcessingResolutionLow"]
            + "-"
            + self.aimless["DataProcessingResolutionHigh"]
        )
        self.aimless["DataProcessingRmergeOverall"] = str(
            round(obj["overall"]["r_merge"], 3)
        )
        self.aimless["DataProcessingRmergeLow"] = str(round(obj["r_merge"][0], 3))
        self.aimless["DataProcessingRmergeHigh"] = str(
            round(obj["r_merge"][len(obj["r_merge"]) - 1], 3)
        )
        self.aimless["DataProcessingIsigOverall"] = str(
            round(obj["overall"]["i_over_sigma_mean"], 1)
        )
        self.aimless["DataProcessingIsigLow"] = str(
            round(obj["i_over_sigma_mean"][0], 1)
        )
        self.aimless["DataProcessingIsigHigh"] = str(
            round(obj["i_over_sigma_mean"][len(obj["i_over_sigma_mean"]) - 1], 1)
        )
        self.aimless["DataProcessingCompletenessOverall"] = str(
            round(obj["overall"]["completeness"], 1)
        )
        self.aimless["DataProcessingCompletenessLow"] = str(
            round(obj["completeness"][0], 1)
        )
        self.aimless["DataProcessingCompletenessHigh"] = str(
            round(obj["completeness"][len(obj["completeness"]) - 1], 1)
        )
        self.aimless["DataProcessingMultiplicityOverall"] = str(
            round(obj["overall"]["multiplicity"], 1)
        )
        self.aimless["DataProcessingMultiplicityLow"] = str(
            round(obj["multiplicity"][0], 1)
        )
        self.aimless["DataProcessingMultiplicityHigh"] = str(
            round(obj["multiplicity"][len(obj["multiplicity"]) - 1], 1)
        )
        self.aimless["DataProcessingCChalfOverall"] = str(
            round(obj["overall"]["cc_one_half"], 2)
        )
        self.aimless["DataProcessingCChalfLow"] = str(round(obj["cc_one_half"][0], 2))
        self.aimless["DataProcessingCChalfHigh"] = str(
            round(obj["cc_one_half"][len(obj["cc_one_half"]) - 1], 2)
        )

        self.aimless["DataProcessingUniqueReflectionsLow"] = str(obj["n_uniq"][0])
        self.aimless["DataProcessingUniqueReflectionsHigh"] = str(
            obj["n_uniq"][len(obj["n_uniq"]) - 1]
        )

        self.aimless["DataProcessingUniqueReflectionsOverall"] = str(
            obj["overall"]["n_obs"]
        )
        json_name = logfile[logfile.rfind("/") + 1 :]
        mmcif_file = logfile.replace("LogFiles", "DataFiles").replace(
            json_name, "xia2.mmcif"
        )
        if os.path.isfile(mmcif_file):
            self.read_mmcif(mmcif_file)
        elif os.path.isfile("%s.bz2" % mmcif_file):
            self.read_mmcif("%s.bz2" % mmcif_file)

    def make_pseudo_aimless_log_from_json(self, logfile):
        self.json_logfile(logfile)
        template = (
            "==============================================================\n"
            "\n"
            "<!--SUMMARY_BEGIN--> $TEXT:Result: $$ $$\n"
            "Summary data for"
            "        Project: nt11175v63"
            " Crystal: xPHIPAx17245"
            " Dataset: SAD\n"
            "\n"
            "                                           "
            "Overall  InnerShell  OuterShell\n"
            "Low resolution limit"
            "                       {0!s}     {1!s}     {2!s}\n".format(
                self.aimless["DataProcessingResolutionLow"],
                self.aimless["DataProcessingResolutionLow"],
                self.aimless["DataProcessingResolutionHighOuterShell"],
            )
            + "High resolution limit"
            "                      {0!s}     {1!s}     {2!s}\n".format(
                self.aimless["DataProcessingResolutionHigh"],
                self.aimless["DataProcessingResolutionLowInnerShell"],
                self.aimless["DataProcessingResolutionHigh"],
            )
            + "\n"
            "Rmerge  (within I+/I-)"
            "                     {0!s}     {1!s}     {2!s}\n".format(
                self.aimless["DataProcessingRmergeOverall"],
                self.aimless["DataProcessingRmergeLow"],
                self.aimless["DataProcessingRmergeHigh"],
            )
            + "Rmerge  (all I+ and I-)"
            "                    {0!s}     {1!s}     {2!s}\n".format(
                self.aimless["DataProcessingRmergeOverall"],
                self.aimless["DataProcessingRmergeLow"],
                self.aimless["DataProcessingRmergeHigh"],
            )
            + "Rmeas (within I+/I-)                          -         -         - \n"
            "Rmeas (all I+ & I-)                           -         -         - \n"
            "Rpim (within I+/I-)                           -         -         - \n"
            "Rpim (all I+ & I-)                            -         -         - \n"
            "Rmerge in top intensity bin                   -         -         - \n"
            "Total number of observations                  -         -         - \n"
            "Total number unique"
            "                        {0!s}     {1!s}     {2!s}\n".format(
                self.aimless["DataProcessingUniqueReflectionsOverall"],
                self.aimless["DataProcessingUniqueReflectionsLow"],
                self.aimless["DataProcessingUniqueReflectionsHigh"],
            )
            + "Mean((I)/sd(I))"
            "                            {0!s}     {1!s}     {2!s}\n".format(
                self.aimless["DataProcessingIsigOverall"],
                self.aimless["DataProcessingIsigLow"],
                self.aimless["DataProcessingIsigHigh"],
            )
            + "Mn(I) half-set correlation CC(1/2)"
            "         {0!s}     {1!s}     {2!s}\n".format(
                self.aimless["DataProcessingCChalfOverall"],
                self.aimless["DataProcessingCChalfLow"],
                self.aimless["DataProcessingCChalfHigh"],
            )
            + "Completeness"
            "                               {0!s}     {1!s}     {2!s}\n".format(
                self.aimless["DataProcessingCompletenessOverall"],
                self.aimless["DataProcessingCompletenessLow"],
                self.aimless["DataProcessingCompletenessHigh"],
            )
            + "Multiplicity"
            "                               {0!s}     {1!s}     {2!s}\n".format(
                self.aimless["DataProcessingMultiplicityOverall"],
                self.aimless["DataProcessingMultiplicityLow"],
                self.aimless["DataProcessingMultiplicityHigh"],
            )
        )
        f = open("aimless_dials.log", "w")
        f.write(template)
        f.close()

    def read_mmcif(self, mmcif):
        for line in open_decompress_file(mmcif):
            if line.startswith("_cell.angle_alpha "):
                self.aimless["DataProcessingAlpha"] = line.split()[1]
            elif line.startswith("_cell.angle_beta "):
                self.aimless["DataProcessingBeta"] = line.split()[1]
            elif line.startswith("_cell.angle_gamma "):
                self.aimless["DataProcessingGamma"] = line.split()[1]
            elif line.startswith("_cell.length_a "):
                self.aimless["DataProcessingA"] = line.split()[1]
            elif line.startswith("_cell.length_b "):
                self.aimless["DataProcessingB"] = line.split()[1]
            elif line.startswith("_cell.length_c "):
                self.aimless["DataProcessingC"] = line.split()[1]
            elif line.startswith("_diffrn_radiation_wavelength.wavelength"):
                self.aimless["DataCollectionWavelength"] = line.split()[1]
            elif line.startswith("_symmetry.space_group_name_H-M"):
                self.aimless["DataProcessingSpaceGroup"] = line[
                    line.find("'") + 1 : line.rfind("'")
                ].replace(" ", "")
                if "R32:H" in self.aimless["DataProcessingSpaceGroup"]:
                    self.aimless["DataProcessingSpaceGroup"] = "H32"
                if "R3:H" in self.aimless["DataProcessingSpaceGroup"]:
                    self.aimless["DataProcessingSpaceGroup"] = "H3"
                self.aimless[
                    "DataProcessingPointGroup"
                ] = self.get_pointgroup_from_space_group(
                    self.aimless["DataProcessingSpaceGroup"]
                )
            elif line.startswith("_space_group.crystal_system"):
                self.aimless["DataProcessingLattice"] = line.split()[1]
                if self.aimless["DataProcessingLattice"] == "trigonal":
                    self.aimless["DataProcessingLattice"] = "hexagonal"
            self.aimless["DataProcessingUnitCell"] = (
                self.aimless["DataProcessingA"]
                + " "
                + self.aimless["DataProcessingA"]
                + " "
                + self.aimless["DataProcessingA"]
                + " "
                + self.aimless["DataProcessingA"]
                + " "
                + self.aimless["DataProcessingA"]
                + " "
                + self.aimless["DataProcessingA"]
            )

    def get_lattice_from_space_group(self, logfile_spacegroup):
        lattice_type = "n/a"
        for lattice in self.space_group_dict:
            for spacegroup in self.space_group_dict[lattice]:
                if logfile_spacegroup.replace(" ", "") == spacegroup:
                    lattice_type = lattice
                    break
        return lattice_type

    def get_pointgroup_from_space_group(self, logfile_spacegroup):
        pointgroup = "n/a"
        for pg in self.point_group_dict:
            for spacegroup in self.point_group_dict[pg]:
                if logfile_spacegroup.replace(" ", "") == spacegroup:
                    pointgroup = pg
                    break
        return pointgroup

    def calc_unitcell_volume_from_logfile(self, a, b, c, alpha, beta, gamma, lattice):
        unitcell_volume = 0
        if lattice == "triclinic":
            unitcell_volume = (
                a
                * b
                * c
                * math.sqrt(
                    (
                        1
                        - math.cos(alpha) ** 2
                        - math.cos(beta) ** 2
                        - math.cos(gamma) ** 2
                    )
                    + 2 * (math.cos(alpha) * math.cos(beta) * math.cos(gamma))
                )
            )
        if "monoclinic" in lattice:
            unitcell_volume = round(a * b * c * math.sin(beta), 1)
        if lattice == "orthorhombic" or lattice == "tetragonal" or lattice == "cubic":
            unitcell_volume = round(a * b * c, 1)
        if lattice == "hexagonal" or lattice == "rhombohedral":
            unitcell_volume = round(a * b * c * (math.sin(math.radians(60))), 1)
        return unitcell_volume

    def PDBheader(self, pdbfile):
        PDBinfo = {
            "Rcryst": "n/a",
            "RcrystTL": "gray",
            "Rfree": "n/a",
            "RfreeTL": "gray",
            "SpaceGroup": "n/a",
            "PointGroup": "n/a",
            "UnitCell": "n/a",
            "ResolutionHigh": "n/a",
            "ResolutionColor": "gray",
            "Lattice": "n/a",
            "UnitCellVolume": 0,
            "Alert": "#E0E0E0",
            "rmsdBonds": "n/a",
            "rmsdBondsTL": "gray",
            "rmsdAngles": "n/a",
            "rmsdAnglesTL": "gray",
            "TwinFraction": "n/a",
        }

        a = "n/a"
        b = "n/a"
        c = "n/a"
        alpha = "n/a"
        beta = "n/a"
        gamma = "n/a"

        if os.path.isfile(pdbfile):
            for line in open(pdbfile):
                try:
                    if line.startswith(
                        "REMARK   3   R VALUE     (WORKING + TEST SET) :"
                    ):
                        PDBinfo["Rcryst"] = line.split()[9]
                        if float(PDBinfo["Rcryst"]) > 0.4:
                            PDBinfo["Alert"] = "#FF0000"
                            PDBinfo["RcrystTL"] = "red"
                        if 0.4 >= float(PDBinfo["Rcryst"]) >= 0.3:
                            PDBinfo["Alert"] = "#FF9900"
                            PDBinfo["RcrystTL"] = "orange"
                        if float(PDBinfo["Rcryst"]) < 0.3:
                            PDBinfo["Alert"] = "#00FF00"
                            PDBinfo["RcrystTL"] = "green"
                    if line.startswith(
                        "REMARK   3   FREE R VALUE                     :"
                    ):
                        PDBinfo["Rfree"] = line.split()[6]
                        if float(PDBinfo["Rfree"]) > 0.4:
                            PDBinfo["Alert"] = "#FF0000"
                            PDBinfo["RfreeTL"] = "red"
                        if 0.4 >= float(PDBinfo["Rfree"]) >= 0.3:
                            PDBinfo["Alert"] = "#FF9900"
                            PDBinfo["RfreeTL"] = "orange"
                        if float(PDBinfo["Rfree"]) < 0.3:
                            PDBinfo["Alert"] = "#00FF00"
                            PDBinfo["RfreeTL"] = "green"
                except ValueError:
                    pass

                if line.startswith("REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) :"):
                    PDBinfo["ResolutionHigh"] = line.split()[7]
                    try:
                        if float(line.split()[7]) < 2.4:
                            PDBinfo["ResolutionColor"] = "green"
                        if 2.4 <= float(line.split()[7]) < 2.8:
                            PDBinfo["ResolutionColor"] = "orange"
                        if float(line.split()[7]) >= 2.8:
                            PDBinfo["ResolutionColor"] = "red"
                    except ValueError:
                        pass
                if line.startswith(
                    "REMARK   3   BOND LENGTHS REFINED ATOMS        (A):"
                ):
                    PDBinfo["rmsdBonds"] = line.split()[9]
                    try:
                        if float(line.split()[9]) < 0.02:
                            PDBinfo["rmsdBondsTL"] = "green"
                        if 0.02 <= float(line.split()[9]) < 0.03:
                            PDBinfo["rmsdBondsTL"] = "orange"
                        if float(line.split()[9]) >= 0.03:
                            PDBinfo["rmsdBondsTL"] = "red"
                    except ValueError:
                        pass
                if line.startswith(
                    "REMARK   3   BOND ANGLES REFINED ATOMS   (DEGREES):"
                ):
                    PDBinfo["rmsdAngles"] = line.split()[9]
                    try:
                        if float(line.split()[9]) < 2.0:
                            PDBinfo["rmsdAnglesTL"] = "green"
                        if 2.0 <= float(line.split()[9]) < 3.0:
                            PDBinfo["rmsdAnglesTL"] = "orange"
                        if float(line.split()[9]) >= 3.0:
                            PDBinfo["rmsdAnglesTL"] = "red"
                    except ValueError:
                        pass

                if line.startswith("REMARK   3      TWIN FRACTION"):
                    try:
                        PDBinfo["TwinFraction"] = line.split()[5]
                    except IndexError:
                        pass

                if line.startswith("CRYST1"):
                    a = int(float(line.split()[1]))
                    b = int(float(line.split()[2]))
                    c = int(float(line.split()[3]))
                    alpha = int(float(line.split()[4]))
                    beta = int(float(line.split()[5]))
                    gamma = int(float(line.split()[6]))

                    PDBinfo["UnitCell"] = (
                        line.split()[1]
                        + " "
                        + line.split()[2]
                        + " "
                        + line.split()[3]
                        + " "
                        + line.split()[4]
                        + " "
                        + line.split()[5]
                        + " "
                        + line.split()[6]
                    )

                    PDBinfo["SpaceGroup"] = str(line[55:65]).rstrip()

                    PDBinfo["Lattice"] = self.get_lattice_from_space_group(
                        PDBinfo["SpaceGroup"]
                    )
                    PDBinfo["PointGroup"] = self.get_pointgroup_from_space_group(
                        PDBinfo["SpaceGroup"]
                    )
                    if (
                        a != "n/a"
                        and b != "n/a"
                        and c != "n/a"
                        and alpha != "n/a"
                        and beta != "n/a"
                        and gamma != "n/a"
                        and PDBinfo["Lattice"] != "n/a"
                    ):
                        PDBinfo[
                            "UnitCellVolume"
                        ] = self.calc_unitcell_volume_from_logfile(
                            float(a),
                            float(b),
                            float(c),
                            math.radians(float(alpha)),
                            math.radians(float(beta)),
                            math.radians(float(gamma)),
                            PDBinfo["Lattice"],
                        )

        return PDBinfo

    def dict_for_datasource_update(self, pdbfile):
        pdb = self.PDBheader(pdbfile)
        db_dict = {
            "RefinementPDB_latest": os.path.realpath(pdbfile),
            "RefinementRcryst": pdb["Rcryst"],
            "RefinementRcrystTraficLight": pdb["RcrystTL"],
            "RefinementRfree": pdb["Rfree"],
            "RefinementRfreeTraficLight": pdb["RfreeTL"],
            "RefinementRmsdBonds": pdb["rmsdBonds"],
            "RefinementRmsdBondsTL": pdb["rmsdBondsTL"],
            "RefinementRmsdAngles": pdb["rmsdAngles"],
            "RefinementRmsdAnglesTL": pdb["rmsdAnglesTL"],
            "RefinementSpaceGroup": pdb["SpaceGroup"],
            "RefinementResolution": pdb["ResolutionHigh"],
            "RefinementResolutionTL": pdb["ResolutionColor"],
        }
        return db_dict

    def update_datasource_with_PDBheader(self, xtal, datasource, pdbfile):
        pdb = self.PDBheader(pdbfile)
        db_dict = {
            "RefinementPDB_latest": os.path.realpath(pdbfile),
            "RefinementRcryst": pdb["Rcryst"],
            "RefinementRcrystTraficLight": pdb["RcrystTL"],
            "RefinementRfree": pdb["Rfree"],
            "RefinementRfreeTraficLight": pdb["RfreeTL"],
            "RefinementRmsdBonds": pdb["rmsdBonds"],
            "RefinementRmsdBondsTL": pdb["rmsdBondsTL"],
            "RefinementRmsdAngles": pdb["rmsdAngles"],
            "RefinementRmsdAnglesTL": pdb["rmsdAnglesTL"],
            "RefinementSpaceGroup": pdb["SpaceGroup"],
            "RefinementResolution": pdb["ResolutionHigh"],
            "RefinementResolutionTL": pdb["ResolutionColor"],
        }
        print(db_dict)
        db = XChemDB.data_source(datasource)
        db.update_data_source(xtal, db_dict)

    def update_datasource_with_phenix_validation_summary(
        self, xtal, datasource, validation_summary
    ):
        db_dict = {}
        if os.path.isfile(validation_summary):
            for line in open(validation_summary):
                if "molprobity score" in line.lower():
                    if len(line.split()) >= 4:
                        db_dict["RefinementMolProbityScore"] = line.split()[3]
                        if float(line.split()[3]) < 2:
                            db_dict["RefinementMolProbityScoreTL"] = "green"
                        if 2 <= float(line.split()[3]) < 3:
                            db_dict["RefinementMolProbityScoreTL"] = "orange"
                        if float(line.split()[3]) >= 3:
                            db_dict["RefinementMolProbityScoreTL"] = "red"

                if "ramachandran outliers" in line.lower():
                    if len(line.split()) >= 4:
                        db_dict["RefinementRamachandranOutliers"] = line.split()[3]
                        if float(line.split()[3]) < 0.3:
                            db_dict["RefinementRamachandranOutliersTL"] = "green"
                        if 0.3 <= float(line.split()[3]) < 1:
                            db_dict["RefinementRamachandranOutliersTL"] = "orange"
                        if float(line.split()[3]) >= 1:
                            db_dict["RefinementRamachandranOutliersTL"] = "red"

                if "favored" in line.lower():
                    if len(line.split()) >= 3:
                        db_dict["RefinementRamachandranFavored"] = line.split()[2]
                        if float(line.split()[2]) < 90:
                            db_dict["RefinementRamachandranFavoredTL"] = "red"
                        if 90 <= float(line.split()[2]) < 98:
                            db_dict["RefinementRamachandranFavoredTL"] = "orange"
                        if float(line.split()[2]) >= 98:
                            db_dict["RefinementRamachandranFavoredTL"] = "green"
        else:
            db_dict["RefinementMolProbityScore"] = "-"
            db_dict["RefinementMolProbityScoreTL"] = "gray"
            db_dict["RefinementRamachandranOutliers"] = "-"
            db_dict["RefinementRamachandranOutliersTL"] = "gray"
            db_dict["RefinementRamachandranFavored"] = "-"
            db_dict["RefinementRamachandranFavoredTL"] = "gray"
        db = XChemDB.data_source(datasource)
        db.update_data_source(xtal, db_dict)


class mtztools:
    def __init__(self, mtzfile):
        self.mtzfile = mtzfile
        self.hkl = any_reflection_file(file_name=self.mtzfile)
        self.miller_arrays = self.hkl.as_miller_arrays()
        self.mtz = self.miller_arrays[0]
        self.iotbxMTZ = mtz.object(self.mtzfile)

        self.space_group_dict = {
            "triclinic": [1],
            "monoclinic_P": [3, 4],
            "monoclinic_C": [5],
            "orthorhombic": [16, 17, 18, 19, 20, 21, 22, 23, 24],
            "tetragonal": [
                75,
                76,
                77,
                78,
                79,
                80,
                89,
                90,
                91,
                92,
                93,
                94,
                95,
                96,
                97,
                98,
            ],
            "hexagonal": [
                143,
                144,
                145,
                149,
                150,
                151,
                152,
                153,
                154,
                168,
                169,
                170,
                171,
                172,
                173,
                177,
                178,
                179,
                180,
                181,
                182,
            ],
            "rhombohedral": [146, 155],
            "cubic": [195, 196, 197, 198, 199, 207, 208, 209, 210, 211, 212, 213, 214],
        }

        self.point_group_dict = {
            "1": [1],
            "2": [3, 4, 5],
            "222": [16, 17, 18, 19, 20, 21, 22, 23, 24],
            "4": [75, 76, 77, 78, 79, 80],
            "422": [89, 90, 91, 92, 93, 94, 95, 96, 97, 98],
            "3": [143, 144, 145, 146],
            "32": [149, 150, 151, 152, 153, 154, 155],
            "6": [168, 169, 170, 171, 172, 173],
            "622": [177, 178, 179, 180, 181, 182],
            "23": [195, 196, 197, 198, 199],
            "432": [207, 208, 209, 210, 211, 212, 213, 214],
        }

        self.nr_asu_in_unitcell_for_point_group = {
            "1": 1,
            "2": 2,
            "222": 4,
            "4": 4,
            "422": 8,
            "3": 3,
            "32": 6,
            "6": 6,
            "622": 12,
            "23": 12,
            "432": 24,
        }

        self.aimless = {
            "DataProcessingProgram": "n/a",
            "DataProcessingSpaceGroup": "n/a",
            "DataProcessingUnitCell": "n/a",
            "DataProcessingA": "n/a",
            "DataProcessingB": "n/a",
            "DataProcessingC": "n/a",
            "DataProcessingAlpha": "n/a",
            "DataProcessingBeta": "n/a",
            "DataProcessingGamma": "n/a",
            "DataProcessingResolutionLow": "n/a",
            "DataProcessingResolutionLowInnerShell": "n/a",
            "DataProcessingResolutionHigh": "n/a",
            "DataProcessingResolutionHighOuterShell": "n/a",
            "DataProcessingResolutionOverall": "n/a",
            "DataProcessingRmergeOverall": "n/a",
            "DataProcessingRmergeLow": "n/a",
            "DataProcessingRmergeHigh": "n/a",
            "DataProcessingIsigOverall": "n/a",
            "DataProcessingIsigLow": "n/a",
            "DataProcessingIsigHigh": "n/a",
            "DataProcessingCompletenessOverall": "n/a",
            "DataProcessingCompletenessLow": "n/a",
            "DataProcessingCompletenessHigh": "n/a",
            "DataProcessingMultiplicityOverall": "n/a",
            "DataProcessingMultiplicityLow": "n/a",
            "DataProcessingMultiplicityHigh": "n/a",
            "DataProcessingCChalfOverall": "n/a",
            "DataProcessingCChalfLow": "n/a",
            "DataProcessingCChalfHigh": "n/a",
            "DataProcessingResolutionHigh15sigma": "n/a",
            "DataProcessingUniqueReflectionsOverall": "n/a",
            "DataProcessingLattice": "n/a",
            "DataProcessingPointGroup": "n/a",
            "DataProcessingUnitCellVolume": 0,
            "DataProcessingAlert": "#FF0000",
            "DataCollectionWavelength": "n/a",
            "DataProcessingScore": "n/a",
        }

    def get_dmin(self):
        return str(round(float(self.mtz.d_min()), 2))

    def get_wavelength(self):
        wavelength = 0.0
        for crystal in self.iotbxMTZ.crystals():
            for dataset in crystal.datasets():
                if not dataset.wavelength() == 0.0:
                    wavelength = str(round(dataset.wavelength(), 5))
                    break
        return wavelength

    def get_information_for_datasource(self):
        db_dict = {}
        mtz_dict = self.get_all_values_as_dict()
        pg = self.get_pointgroup_from_mtz()
        if mtz_dict != {}:
            db_dict["DataProcessingResolutionHigh"] = mtz_dict["resolution_high"]
            db_dict["DataProcessingUnitCell"] = mtz_dict["unitcell"]
            db_dict["DataProcessingSpaceGroup"] = mtz_dict["spacegroup"]
            db_dict["DataProcessingUnitCellVolume"] = mtz_dict["unitcell_volume"]
            db_dict["DataProcessingLattice"] = mtz_dict["bravais_lattice"]
        if pg != "":
            db_dict["DataProcessingPointGroup"] = pg
        return db_dict

    def get_bravais_lattice_from_spg_number(self, number):
        lattice = ""
        for bravaislattice in self.space_group_dict:
            for spg_number in self.space_group_dict[bravaislattice]:
                if spg_number == number:
                    lattice = bravaislattice
        return lattice

    def get_point_group_from_spg_number(self, number):
        pointgroup = ""
        for pg in self.point_group_dict:
            for spg_number in self.point_group_dict[pg]:
                if spg_number == number:
                    pointgroup = pg
        return pointgroup

    def get_spg_number_from_mtz(self):
        spg_number = 0
        mtzdmp = subprocess.Popen(["mtzdmp", self.mtzfile], stdout=subprocess.PIPE)
        for n, line in enumerate(iter(mtzdmp.stdout.readline, "")):
            if line.startswith(" * Space group ="):
                spg_number = int(line[line.rfind(")") - 3 : line.rfind(")")])
        return spg_number

    def get_pointgroup_from_mtz(self):
        pointgroup = ""
        spg_number = self.get_spg_number_from_mtz()
        pointgroup = self.get_point_group_from_spg_number(spg_number)
        return pointgroup

    def get_number_measured_reflections(self):
        missing_reflections = "0"
        all_reflections = "0"
        meassured_reflections = "0"
        mtzdmp = subprocess.Popen(["mtzdmp", self.mtzfile], stdout=subprocess.PIPE)
        foundTable = False
        for n, line in enumerate(iter(mtzdmp.stdout.readline, "")):
            if line.startswith(
                " Col Sort    Min    Max    Num      %"
                "     Mean     Mean   Resolution   Type Column"
            ):
                foundTable = True
            if foundTable and len(line.split()) == 12:
                if line.split()[11] == "F":
                    missing_reflections = line.split()[4]
                    foundTable = False
            if line.startswith(" No. of reflections used in FILE STATISTICS"):
                all_reflections = line.split()[7]
                break
        try:
            meassured_reflections = int(all_reflections) - int(missing_reflections)
        except ValueError:
            pass
        return meassured_reflections

    def calculate_correlaton_between_intensities_in_mtzfiles(self, mtzin):
        CC = "0.0"
        errorMessage = ""
        cmd = (
            "pointless hklin %s hklref %s << eof\n" % (mtzin, self.mtzfile)
            + "labref I=IMEAN\n"
            "labin I=IMEAN\n"
            "eof\n"
        )

        pointless = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        foundLine = False
        for line in iter(pointless.stdout.readline, ""):
            if foundLine:
                CC = line.split()[3]
                break
            if "Alternative reindexing        Lklhd      CC" in line:
                foundLine = True
            if "**** Incompatible symmetries ****" in line:
                errorMessage = "**** Incompatible symmetries ****"
                break
            if (
                "Merged test dataset (HKLIN)"
                " has different Laue symmetry to reference set" in line
            ):
                errorMessage = "%s has different Laue symmetry to %s" % (
                    mtzin,
                    self.mtzfile,
                )
                break

        return CC, errorMessage

    def get_all_values_as_dict(self):
        mtz = {
            "resolution_high": "n/a",
            "unitcell": "n/a",
            "spacegroup": "n/a",
            "unitcell_volume": "n/a",
            "bravais_lattice": "n/a",
        }
        a = 0.0
        b = 0.0
        c = 0.0
        alpha_rad = 0.0
        beta_rad = 0.0
        gamma_rad = 0.0
        resolution_line = 1000000
        cell_line = 1000000
        mtzdmp = subprocess.Popen(["mtzdmp", self.mtzfile], stdout=subprocess.PIPE)
        for n, line in enumerate(iter(mtzdmp.stdout.readline, "")):
            if line.startswith(" *  Resolution Range :"):
                resolution_line = n + 2
            if n == resolution_line and len(line.split()) == 8:
                mtz["resolution_high"] = round(float(line.split()[5]), 2)
            if line.startswith(" * Cell Dimensions :"):
                cell_line = n + 2
            if n == cell_line and len(line.split()) == 6:
                a = round(float(line.split()[0]), 1)
                b = round(float(line.split()[1]), 1)
                c = round(float(line.split()[2]), 1)
                alpha = round(float(line.split()[3]), 1)
                beta = round(float(line.split()[4]), 1)
                gamma = round(float(line.split()[5]), 1)
                mtz["unitcell"] = (
                    str(a)
                    + " "
                    + str(b)
                    + " "
                    + str(c)
                    + " "
                    + str(alpha)
                    + " "
                    + str(beta)
                    + " "
                    + str(gamma)
                )
                alpha_rad = math.radians(alpha)
                beta_rad = math.radians(beta)
                gamma_rad = math.radians(gamma)
            if line.startswith(" * Space group ="):
                spg_number = int(line[line.rfind(")") - 3 : line.rfind(")")])
                mtz["bravais_lattice"] = self.get_bravais_lattice_from_spg_number(
                    spg_number
                )
                mtz["spacegroup"] = line[line.find("'") + 1 : line.rfind("'")]
        if mtz["bravais_lattice"] == "triclinic":
            mtz["unitcell_volume"] = (
                a
                * b
                * c
                * math.sqrt(
                    (
                        1
                        - math.cos(alpha_rad) ** 2
                        - math.cos(beta_rad) ** 2
                        - math.cos(gamma_rad) ** 2
                    )
                    + 2
                    * (math.cos(alpha_rad) * math.cos(beta_rad) * math.cos(gamma_rad))
                )
            )
        elif "monoclinic" in mtz["bravais_lattice"]:
            mtz["unitcell_volume"] = round(a * b * c * math.sin(beta_rad), 1)
        elif (
            mtz["bravais_lattice"] == "orthorhombic"
            or mtz["bravais_lattice"] == "tetragonal"
            or mtz["bravais_lattice"] == "cubic"
        ):
            mtz["unitcell_volume"] = round(a * b * c, 1)
        elif (
            mtz["bravais_lattice"] == "hexagonal"
            or mtz["bravais_lattice"] == "rhombohedral"
        ):
            mtz["unitcell_volume"] = round(a * b * c * (math.sin(math.radians(60))), 1)

        return mtz

    def get_all_columns_as_dict(self):
        column_dict = {"F": [], "I": [], "SIG": [], "PHS": [], "FOM": [], "RFREE": []}
        startline = 1000000
        mtzdmp = subprocess.Popen(["mtzdmp", self.mtzfile], stdout=subprocess.PIPE)
        for n, line in enumerate(iter(mtzdmp.stdout.readline, "")):
            if line.startswith(" Col Sort    Min    Max    Num"):
                startline = n + 2
            if n >= startline and len(line.split()) > 10:
                if line.split()[10] == "F":
                    column_dict["F"].append(line.split()[11])
                if line.split()[10] == "J":
                    column_dict["I"].append(line.split()[11])
                if line.split()[10] == "Q":
                    column_dict["SIG"].append(line.split()[11])
                if line.split()[10] == "I":
                    column_dict["RFREE"].append(line.split()[11])
                if line.split()[10] == "P":
                    column_dict["PHS"].append(line.split()[11])
                if line.split()[10] == "W":
                    column_dict["FOM"].append(line.split()[11])

        return column_dict

    def get_all_columns_as_list(self):
        column_list = self.iotbxMTZ.column_labels()
        return column_list


class external_software:
    def __init__(self, xce_logfile):
        self.available_programs = {}
        self.Logfile = XChemLog.updateLog(xce_logfile)

    def log_found_status(self, program_name):
        self.Logfile.insert(
            "{0:50} {1:10}".format(
                "checking for {}:".format(program_name),
                "found" if self.available_programs[program_name] else "not found",
            )
        )

    def check(self):
        self.Logfile.insert("Searching for external software...")

        # default is False; user needs to explicitely set this
        self.available_programs["qsub_remote"] = False

        self.available_programs["qsub"] = self.available_programs["qsub_array"] = (
            find_executable("qsub") is not None
        )
        self.log_found_status("qsub")
        self.log_found_status("qsub_array")

        self.available_programs["refmac5"] = find_executable("refmac5") is not None
        self.log_found_status("refmac5")

        self.available_programs["phenix.molprobity"] = (
            find_executable("phenix.molprobity") is not None
        )
        self.log_found_status("phenix.molprobity")

        self.available_programs["phenix.find_tls_groups"] = (
            find_executable("phenix.find_tls_groups") is not None
        )
        self.log_found_status("phenix.find_tls_groups")

        self.available_programs["mmtbx.validate_ligands"] = (
            find_executable("mmtbx.validate_ligands") is not None
        )
        self.log_found_status("mmtbx.validate_ligands")

        self.available_programs["acedrg"] = find_executable("acedrg") is not None
        self.log_found_status("acedrg")

        self.available_programs["phenix.elbow"] = (
            find_executable("phenix.elbow") is not None
        )
        self.log_found_status("phenix.elbow")

        self.available_programs["grade"] = find_executable("grade") is not None
        self.log_found_status("grade")

        self.available_programs["giant.create_occupancy_params"] = (
            find_executable("giant.create_occupancy_params") is not None
        )
        self.log_found_status("giant.create_occupancy_params")

        self.available_programs["mogul"] = (
            "BDG_TOOL_MOGUL" in os.environ
            and find_executable(os.environ["BDG_TOOL_MOGUL"]) is not None
        )
        self.log_found_status("mogul")

        self.available_programs["gemmi"] = find_executable("gemmi") is not None
        self.log_found_status("gemmi")

        return self.available_programs


class pdbtools(object):
    def __init__(self, pdb):
        self.pdb = pdb
        self.pdb_inp = iotbx.pdb.input(file_name=self.pdb)
        self.hierarchy = self.pdb_inp.construct_hierarchy()

        self.AminoAcids = [
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLU",
            "GLN",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
            "CSO",
            "HYP",
        ]
        self.Solvents = ["DMS", "EDO", "GOL", "HOH"]
        self.Ions = ["NA", "MG", "CL", "K", "SO4", "PO4", "CA"]
        self.xce_ligands = ["LIG", "DRG", "FRS"]

        self.space_group_dict = {
            "triclinic": [1],
            "monoclinic_P": [3, 4],
            "monoclinic_C": [5],
            "orthorhombic": [16, 17, 18, 19, 20, 21, 22, 23, 24],
            "tetragonal": [
                75,
                76,
                77,
                78,
                79,
                80,
                89,
                90,
                91,
                92,
                93,
                94,
                95,
                96,
                97,
                98,
            ],
            "hexagonal": [
                143,
                144,
                145,
                149,
                150,
                151,
                152,
                153,
                154,
                168,
                169,
                170,
                171,
                172,
                173,
                177,
                178,
                179,
                180,
                181,
                182,
            ],
            "rhombohedral": [146, 155],
            "cubic": [195, 196, 197, 198, 199, 207, 208, 209, 210, 211, 212, 213, 214],
        }

        self.point_group_dict = {
            "1": [1],
            "2": [3, 4, 5],
            "222": [16, 17, 18, 19, 20, 21, 22, 23, 24],
            "4": [75, 76, 77, 78, 79, 80],
            "422": [89, 90, 91, 92, 93, 94, 95, 96, 97, 98],
            "3": [143, 144, 145, 146],
            "32": [149, 150, 151, 152, 153, 154, 155],
            "6": [168, 169, 170, 171, 172, 173],
            "622": [177, 178, 179, 180, 181, 182],
            "23": [195, 196, 197, 198, 199],
            "432": [207, 208, 209, 210, 211, 212, 213, 214],
        }

        self.nr_asu_in_unitcell_for_point_group = {
            "1": 1,
            "2": 2,
            "222": 4,
            "4": 4,
            "422": 8,
            "3": 3,
            "32": 6,
            "6": 6,
            "622": 12,
            "23": 12,
            "432": 24,
        }

    def amino_acids(self):
        return self.AminoAcids

    def get_refinement_program(self):
        program = "unknown"
        for remark in self.pdb_inp.remark_section():
            if "PROGRAM" in remark:
                if "refmac" in remark.lower():
                    program = "REFMAC"
                elif "phenix" in remark.lower():
                    program = "PHENIX"
                elif "buster" in remark.lower():
                    program = "BUSTER"
        return program

    def get_residues_with_resname(self, resname):
        ligands = []
        for model in self.hierarchy.models():
            for chain in model.chains():
                for conformer in chain.conformers():
                    for residue in conformer.residues():
                        if residue.resname == resname:
                            if [
                                residue.resname,
                                residue.resseq,
                                chain.id,
                            ] not in ligands:
                                ligands.append(
                                    [residue.resname, residue.resseq, chain.id]
                                )
        return ligands

    def save_residues_with_resname(self, outDir, resname):
        ligands = self.get_residues_with_resname(resname)
        ligList = []
        for ligand in ligands:
            sel_cache = self.hierarchy.atom_selection_cache()
            lig_sel = sel_cache.selection(
                "(resname %s and resseq %s and chain %s)"
                % (ligand[0], ligand[1], ligand[2])
            )
            hierarchy_lig = self.hierarchy.select(lig_sel)
            ligName = (ligand[0] + "-" + ligand[2] + "-" + ligand[1] + ".pdb").replace(
                " ", ""
            )
            ligList.append(ligName)

            try:
                f = open(os.path.join(outDir, ligName), "w")
                f.write(
                    hierarchy_lig.as_pdb_string(
                        crystal_symmetry=self.pdb_inp.crystal_symmetry()
                    )
                )
                f.close()
            except IOError:
                print("ERROR: {0!s} exists; skipping...")

        return ligList

    def GetProteinChains(self):
        chain = []
        for line in open(self.pdb):
            if line.startswith("ATOM"):
                if line[17:20] in self.AminoAcids:
                    if line[21:22] not in chain:
                        chain.append(line[21:22])
        return chain

    def get_bravais_lattice_from_spg_number(self, number):
        lattice = ""
        for bravaislattice in self.space_group_dict:
            for spg_number in self.space_group_dict[bravaislattice]:
                if str(spg_number) == str(number):
                    lattice = bravaislattice
        return lattice

    def find_ligands(self):
        Ligands = []
        for line in open(self.pdb):
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[
                17:20
            ].replace(" ", "") not in self.AminoAcids + self.Solvents + self.Ions:
                if [line[17:20], line[21:22], line[23:26]] not in Ligands:
                    Ligands.append([line[17:20], line[21:22], line[23:26]])
        return Ligands

    def save_ligands_to_pdb(self):
        Ligands = self.find_ligands()
        if not Ligands == []:
            for n, item in enumerate(Ligands):
                pdb = ""
                for line in open(self.pdb):
                    if line.startswith("CRYST"):
                        pdb += line
                    if (
                        (line.startswith("ATOM") or line.startswith("HETATM"))
                        and line[17:20] == item[0]
                        and line[21:22] == item[1]
                        and line[23:26] == item[2]
                    ):
                        pdb = pdb + line
                f = open("ligand_{0!s}.pdb".format(n), "w")
                f.write(pdb)
                f.close()
        return Ligands

    def find_xce_ligand_details(self):
        Ligands = []
        for line in open(self.pdb):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resname = str(line[17:20]).replace(" ", "")
                if resname in self.xce_ligands:
                    chainID = str(line[21:23]).replace(" ", "")
                    resseq = str(line[23:26]).replace(" ", "")
                    altLoc = str(line[16:17]).replace(" ", "")
                    if [resname, chainID, resseq, altLoc] not in Ligands:
                        Ligands.append([resname, chainID, resseq, altLoc])
        return Ligands

    def ligand_details_as_list(self):
        Ligands = []
        for line in open(self.pdb):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resname = str(line[17:20]).replace(" ", "")
                if resname in self.xce_ligands:
                    chainID = str(line[21:23]).replace(" ", "")
                    resseq = str(line[23:26]).replace(" ", "")
                    altLoc = str(line[16:17]).replace(" ", "")
                    occupancy = str(line[56:60]).replace(" ", "")
                    if [resname, chainID, resseq, altLoc, occupancy] not in Ligands:
                        Ligands.append([resname, chainID, resseq, altLoc, occupancy])
        return Ligands

    def get_center_of_gravity_of_residue_ish(self, chain, number):
        print("-> chain:", chain)
        print("-> number:", number)
        X = 0.0
        Y = 0.0
        Z = 0.0
        x_list = []
        y_list = []
        z_list = []
        # pdb definition see:
        # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        for line in open(self.pdb):
            if (
                (line.startswith("ATOM") or line.startswith("HETATM"))
                and line[21:22].replace(" ", "") == chain.replace(" ", "")
                and line[22:26].replace(" ", "") == str(number).replace(" ", "")
            ):
                X = float(line[30:38])
                x_list.append(X)
                Y = float(line[38:46])
                y_list.append(Y)
                Z = float(line[46:54])
                z_list.append(Z)
        # 'ish' because it's not really the centre of gravity
        # but the the middle of the min/max of each x,y,z
        X = ((max(x_list) - min(x_list)) / 2) + min(x_list)
        Y = ((max(y_list) - min(y_list)) / 2) + min(y_list)
        Z = ((max(z_list) - min(z_list)) / 2) + min(z_list)
        return X, Y, Z

    def get_center_of_gravity_of_molecule_ish(self):
        X = 0.0
        Y = 0.0
        Z = 0.0
        x_list = []
        y_list = []
        z_list = []
        # pdb definition see:
        # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        for line in open(self.pdb):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                X = float(line[30:38])
                x_list.append(X)
                Y = float(line[38:46])
                y_list.append(Y)
                Z = float(line[46:54])
                z_list.append(Z)
        # 'ish' because it's not really the centre of gravity
        # but the the middle of the min/max of each x,y,z
        X = ((max(x_list) - min(x_list)) / 2) + min(x_list)
        Y = ((max(y_list) - min(y_list)) / 2) + min(y_list)
        Z = ((max(z_list) - min(z_list)) / 2) + min(z_list)
        return X, Y, Z

    def ElementDict(self, resname, chainID, resseq, altLoc):
        ElementDict = {
            "C": 0,
            "N": 0,
            "O": 0,
            "P": 0,
            "S": 0,
            "BR": 0,
            "CL": 0,
            "I": 0,
            "F": 0,
        }

        for line in open(self.pdb):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resname_line = str(line[17:20]).replace(" ", "")
                chainID_line = str(line[21:23]).replace(" ", "")
                resseq_line = str(line[23:26]).replace(" ", "")
                altLoc_line = str(line[16:17]).replace(" ", "")
                element_line = str(line[66:78]).replace(" ", "")
                if (
                    resname_line == resname
                    and chainID_line == chainID
                    and resseq_line == resseq
                    and altLoc_line == altLoc
                ):
                    if element_line.upper() in ElementDict:
                        ElementDict[element_line.upper()] += 1

        return ElementDict


class logtools:
    def __init__(self, logfile):
        self.logfile = logfile

    def phenix_molprobity(self):
        QualityIndicators = {
            "MolprobityScore": "n/a",
            "MolprobityScoreColor": "gray",
            "RamachandranOutliers": "n/a",
            "RamachandranOutliersColor": "gray",
            "RamachandranFavored": "n/a",
            "RamachandranFavoredColor": "gray",
        }

        # Molprobity = validation_summary.txt
        if os.path.isfile(self.logfile):
            for line in open(self.logfile):
                if "molprobity score" in line.lower():
                    if len(line.split()) >= 4:
                        QualityIndicators["MolprobityScore"] = line.split()[3]
                        try:
                            if float(line.split()[3]) < 2:
                                QualityIndicators["MolprobityScoreColor"] = "green"
                            if 2 <= float(line.split()[3]) < 3:
                                QualityIndicators["MolprobityScoreColor"] = "orange"
                            if float(line.split()[3]) >= 3:
                                QualityIndicators["MolprobityScoreColor"] = "red"
                        except ValueError:
                            pass
                if "ramachandran outliers" in line.lower():
                    if len(line.split()) >= 4:
                        QualityIndicators["RamachandranOutliers"] = line.split()[3]
                        try:
                            if float(line.split()[3]) < 0.3:
                                QualityIndicators["RamachandranOutliersColor"] = "green"
                            if 0.3 <= float(line.split()[3]) < 1:
                                QualityIndicators[
                                    "RamachandranOutliersColor"
                                ] = "orange"
                            if float(line.split()[3]) >= 1:
                                QualityIndicators["RamachandranOutliersColor"] = "red"
                        except ValueError:
                            pass
                if "favored" in line.lower():
                    if len(line.split()) == 4:
                        QualityIndicators["RamachandranFavored"] = line.split()[2]
                        try:
                            if float(line.split()[2]) < 90:
                                QualityIndicators["RamachandranFavoredColor"] = "red"
                            if 90 <= float(line.split()[2]) < 98:
                                QualityIndicators["RamachandranFavoredColor"] = "orange"
                            if float(line.split()[2]) >= 98:
                                QualityIndicators["RamachandranFavoredColor"] = "green"
                        except ValueError:
                            pass

        return QualityIndicators

    def refmac_log(self):
        QualityIndicators = {"MatrixWeight": "n/a"}

        # Matrix Weight
        if os.path.isfile(self.logfile):
            for line in open(self.logfile):
                if line.startswith(" Weight matrix") and len(line.split()) == 3:
                    QualityIndicators["MatrixWeight"] = line.split()[2]

        return QualityIndicators


def calculate_distance_between_coordinates(x1, y1, z1, x2, y2, z2):
    distance = 0.0
    distance = math.sqrt(
        math.pow(float(x1) - float(x2), 2)
        + math.pow(float(y1) - float(y2), 2)
        + math.pow(float(z1) - float(z2), 2)
    )
    return distance


class smilestools(object):
    def __init__(self, smiles):
        self.smiles = smiles

    def ElementDict(self):
        ElementDict = {
            "C": 0,
            "N": 0,
            "O": 0,
            "P": 0,
            "S": 0,
            "BR": 0,
            "CL": 0,
            "I": 0,
            "F": 0,
        }

        m = Chem.MolFromSmiles(self.smiles)
        for atom in m.GetAtoms():
            if str(atom.GetSymbol()).upper() in ElementDict:
                ElementDict[str(atom.GetSymbol()).upper()] += 1

        return ElementDict


class maptools(object):
    def calculate_map(self, mtz, F, PH):
        cmd = (
            "fft hklin %s mapout %s << EOF\n" % (mtz, mtz.replace(".mtz", ".ccp4"))
            + "labin F1=%s PHI=%s\n" % (F, PH)
            + "EOF\n"
        )
        os.system(cmd)

    def cut_map_around_ligand(self, map, ligPDB, border):
        if map.endswith(".map"):
            map_extension = ".map"
        elif map.endswith(".ccp4"):
            map_extension = ".ccp4"
        else:
            map_extension = ""

        cmd = (
            "mapmask mapin %s mapout %s xyzin %s << eof\n"
            % (map, map.replace(map_extension, "_mapmask" + map_extension), ligPDB)
            + " border %s\n" % border
            + " end\n"
            "eof"
        )
        os.system(cmd)


class mtztools_gemmi:
    def __init__(self, mtz):
        self.mtz = gemmi.read_mtz_file(mtz)

    def get_map_labels(self):
        labelList = []
        for column in self.mtz.columns:
            labelList.append(column.label)
        FWT = None
        PHWT = None
        DELFWT = None
        PHDELWT = None

        if "FWT" in labelList and "PHWT" in labelList:
            FWT = "FWT"
            PHWT = "PHWT"

        if "DELFWT" in labelList and "PHDELWT" in labelList:
            DELFWT = "DELFWT"
            PHDELWT = "PHDELWT"

        if "2FOFCWT" in labelList and "PH2FOFCWT" in labelList:
            FWT = "2FOFCWT"
            PHWT = "PH2FOFCWT"

        if "FOFCWT" in labelList and "PHFOFCWT" in labelList:
            DELFWT = "FOFCWT"
            PHDELWT = "PHFOFCWT"

        return FWT, PHWT, DELFWT, PHDELWT

    def get_high_low_resolution_limits(self):
        resl = self.mtz.resolution_low()
        resh = self.mtz.resolution_high()
        return resh, resl


class maptools_gemmi:
    def __init__(self, emap):
        self.emap = emap
        self.emtz = emap.replace(".ccp4", ".mtz").replace(".map", ".mtz")

    def map_to_sf(self, resolution):
        if os.path.isfile(self.emtz):
            print("mtz file of event map exists; skipping...")
            return
        cmd = "gemmi map2sf %s %s FWT PHWT --dmin=%s" % (
            self.emap,
            self.emtz,
            resolution,
        )
        print(("converting map with command:\n" + cmd))
        os.system(cmd)
        if os.path.isfile(self.emtz):
            print("event map to SF conversion successful")
            mtz = gemmi.read_mtz_file(self.emtz)
            mtz.history += ["date created: " + time.ctime(os.path.getmtime(self.emap))]
            mtz.history += ["folder: " + os.getcwd()]
            mtz.history += ["file name: " + self.emap]
            if "BDC" in self.emap:
                mtz.history += [
                    "BDC value: "
                    + self.emap[
                        self.emap.find("BDC")
                        + 4 : self.emap.find("BDC")
                        + 4
                        + self.emap[self.emap.find("BDC") + 4 :].find("_")
                    ]
                ]
        else:
            print("failed to convert event map to SF")


class pdbtools_gemmi:
    def __init__(self, pdb):
        self.pdb = gemmi.read_structure(pdb)

    def get_ligand_models_as_dict(self, ligandID):
        ligandDict = {}
        for model in self.pdb:
            for chain in model:
                for residue in chain:
                    if residue.name == ligandID:
                        if (
                            str(
                                residue.name
                                + "-"
                                + chain.name
                                + "-"
                                + str(residue.seqid.num)
                            )
                            not in ligandDict
                        ):
                            ligandDict[
                                str(
                                    residue.name
                                    + "-"
                                    + chain.name
                                    + "-"
                                    + str(residue.seqid.num)
                                )
                            ] = None
                            m = gemmi.Model("1")
                            m.add_chain(gemmi.Chain("X"))
                            c = m["X"].get_polymer()
                            c.add_residue(gemmi.Residue(), 0)
                            c[0].name = residue.name
                            for n, atom in enumerate(residue):
                                c[0].add_atom(atom, n)
                            ligandDict[
                                str(
                                    residue.name
                                    + "-"
                                    + chain.name
                                    + "-"
                                    + str(residue.seqid.num)
                                )
                            ] = m
        return ligandDict

    def center_of_mass_ligand_dict(self, ligandID):
        ligandDict = self.get_ligand_models_as_dict(ligandID)
        ligandPositionDict = {}
        for ligand in ligandDict:
            pos = ligandDict[ligand].calculate_center_of_mass()
            ligandPositionDict[ligand] = [pos.x, pos.y, pos.z]
        return ligandPositionDict

    def save_ligands_to_pdb(self, ligandID):
        ligandDict = self.get_ligand_models_as_dict(ligandID)
        for ligand in ligandDict:
            s = gemmi.Structure()
            s.add_model(ligandDict[ligand])
            s.write_pdb(ligand + ".pdb")
