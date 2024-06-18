import os
import pickle

import coot
import gobject
import gtk
from matplotlib.figure import Figure

from xce.lib import coot_utils_XChem
from xce.lib import XChemDB
from xce.lib import XChemRefine
from xce.lib import XChemUtils

# had to adapt the original coot_utils.py file
# otherwise unable to import the original file without complaints about missing modules
# etc.
# modified file is now in $XChemExplorer_DIR/lib


class GUI(object):

    """
    main class which opens the actual GUI
    """

    def __init__(self):
        ################################################################################
        # read in settings file from XChemExplorer to set the relevant paths
        print("current dir", os.getcwd())
        self.settings = pickle.load(open(".xce_settings.pkl", "rb"))
        print("setting", self.settings)
        self.database_directory = self.settings["database_directory"]
        self.xce_logfile = self.settings["xce_logfile"]
        self.data_source = self.settings["data_source"]
        self.db = XChemDB.data_source(self.data_source)

        # checking for external software packages
        self.external_software = XChemUtils.external_software(self.xce_logfile).check()

        self.selection_criteria = [
            "0 - All Datasets",
            "1 - Analysis Pending",
            "2 - PANDDA model",
            "3 - In Refinement",
            "4 - CompChem ready",
            "5 - Deposition ready",
            "6 - Deposited",
        ]

        self.experiment_stage = [
            ["Review PANDDA export", "2 - PANDDA model", 65000, 0, 0],
            ["In Refinement", "3 - In Refinement", 65000, 0, 0],
            ["Comp Chem Ready!", "4 - CompChem ready", 65000, 0, 0],
            ["Ready for Deposition!", "5 - Deposition ready", 65000, 0, 0],
            ["In PDB", "6 - Deposited", 65000, 0, 0],
        ]

        self.ligand_confidence_category = [
            "0 - no ligand present",
            "1 - Low Confidence",
            "2 - Correct ligand, weak density",
            "3 - Clear density, unexpected ligand",
            "4 - High Confidence",
        ]

        # this decides which samples will be looked at
        self.selection_mode = ""
        self.pandda_index = -1  # refers to the number of sites
        self.site_index = "0"
        self.event_index = "0"

        # the Folder is kind of a legacy thing because my inital idea was to have
        # separate folders for Data Processing and Refinement
        self.project_directory = self.settings["initial_model_directory"]
        self.Serial = 0
        self.Refine = None
        self.index = -1
        self.Todo = []
        self.siteDict = {}

        self.xtalID = ""
        self.compoundID = ""
        self.spider_plot = ""
        self.ligand_confidence = ""

        self.pdb_style = "refine.pdb"
        self.mtz_style = "refine.mtz"

        # stores imol of currently loaded molecules and maps
        self.mol_dict = {
            "protein": -1,
            "ligand": -1,
            "2fofc": -1,
            "fofc": -1,
            "event": -1,
        }

        # two dictionaries which are flushed when a new crystal is loaded
        # and which contain information to update the data source if necessary
        self.db_dict_mainTable = {}
        self.db_dict_panddaTable = {}

        ################################################################################
        # some COOT settings
        coot.set_map_radius(15)
        coot.set_colour_map_rotation_for_map(0)
        coot.set_colour_map_rotation_on_read_pdb_flag(0)

        self.QualityIndicators = {
            "RefinementRcryst": "-",
            "RefinementRfree": "-",
            "RefinementRfreeTraficLight": "gray",
            "RefinementResolution": "-",
            "RefinementResolutionTL": "gray",
            "RefinementMolProbityScore": "-",
            "RefinementMolProbityScoreTL": "gray",
            "RefinementRamachandranOutliers": "-",
            "RefinementRamachandranOutliersTL": "gray",
            "RefinementRamachandranFavored": "-",
            "RefinementRamachandranFavoredTL": "gray",
            "RefinementRmsdBonds": "-",
            "RefinementRmsdBondsTL": "gray",
            "RefinementRmsdAngles": "-",
            "RefinementRmsdAnglesTL": "gray",
            "RefinementMatrixWeight": "-",
        }

        self.spider_plot_data = {
            "PANDDA_site_ligand_id": "-",
            "PANDDA_site_occupancy": "-",
            "PANDDA_site_B_average": "-",
            "PANDDA_site_B_ratio_residue_surroundings": "-",
            "PANDDA_site_RSCC": "-",
            "PANDDA_site_rmsd": "-",
            "PANDDA_site_RSR": "-",
            "PANDDA_site_RSZD": "-",
        }

        # default refmac parameters
        self.RefmacParams = {
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

    def StartGUI(self):
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", gtk.main_quit)
        self.window.set_border_width(10)
        self.window.set_default_size(400, 800)
        self.window.set_title("XChemExplorer")
        self.vbox = gtk.VBox()  # this is the main container

        ################################################################################
        # --- Sample Selection ---

        frame = gtk.Frame(label="Select Samples")
        self.hbox_select_samples = gtk.HBox()

        self.cb_select_samples = gtk.combo_box_new_text()
        self.cb_select_samples.connect("changed", self.set_selection_mode)
        for citeria in self.selection_criteria:
            self.cb_select_samples.append_text(citeria)
        self.hbox_select_samples.add(self.cb_select_samples)

        self.select_samples_button = gtk.Button(label="GO")
        self.select_samples_button.connect("clicked", self.get_samples_to_look_at)
        self.hbox_select_samples.add(self.select_samples_button)
        frame.add(self.hbox_select_samples)
        self.vbox.pack_start(frame)

        ################################################################################
        # --- status window ---
        frame = gtk.Frame()
        self.status_label = gtk.Label()
        frame.add(self.status_label)
        self.vbox.pack_start(frame)

        # SPACER
        self.vbox.add(gtk.Label(" "))

        ################################################################################
        # --- Refinement Statistics ---
        # next comes a section which displays some global quality indicators
        # a combination of labels and textview widgets, arranged in a table

        RRfreeLabel_frame = gtk.Frame()
        self.RRfreeLabel = gtk.Label("R/Rfree")
        RRfreeLabel_frame.add(self.RRfreeLabel)
        self.RRfreeValue = gtk.Label(
            self.QualityIndicators["RefinementRcryst"]
            + "/"
            + self.QualityIndicators["RefinementRfree"]
        )
        RRfreeBox_frame = gtk.Frame()
        self.RRfreeBox = gtk.EventBox()
        self.RRfreeBox.add(self.RRfreeValue)
        RRfreeBox_frame.add(self.RRfreeBox)

        ResolutionLabel_frame = gtk.Frame()
        self.ResolutionLabel = gtk.Label("Resolution")
        ResolutionLabel_frame.add(self.ResolutionLabel)
        self.ResolutionValue = gtk.Label(self.QualityIndicators["RefinementResolution"])
        ResolutionBox_frame = gtk.Frame()
        self.ResolutionBox = gtk.EventBox()
        self.ResolutionBox.add(self.ResolutionValue)
        ResolutionBox_frame.add(self.ResolutionBox)

        MolprobityScoreLabel_frame = gtk.Frame()
        self.MolprobityScoreLabel = gtk.Label("MolprobityScore")
        MolprobityScoreLabel_frame.add(self.MolprobityScoreLabel)
        self.MolprobityScoreValue = gtk.Label(
            self.QualityIndicators["RefinementMolProbityScore"]
        )
        MolprobityScoreBox_frame = gtk.Frame()
        self.MolprobityScoreBox = gtk.EventBox()
        self.MolprobityScoreBox.add(self.MolprobityScoreValue)
        MolprobityScoreBox_frame.add(self.MolprobityScoreBox)

        RamachandranOutliersLabel_frame = gtk.Frame()
        self.RamachandranOutliersLabel = gtk.Label("Rama Outliers")
        RamachandranOutliersLabel_frame.add(self.RamachandranOutliersLabel)
        self.RamachandranOutliersValue = gtk.Label(
            self.QualityIndicators["RefinementRamachandranOutliers"]
        )
        RamachandranOutliersBox_frame = gtk.Frame()
        self.RamachandranOutliersBox = gtk.EventBox()
        self.RamachandranOutliersBox.add(self.RamachandranOutliersValue)
        RamachandranOutliersBox_frame.add(self.RamachandranOutliersBox)

        RamachandranFavoredLabel_frame = gtk.Frame()
        self.RamachandranFavoredLabel = gtk.Label("Rama Favored")
        RamachandranFavoredLabel_frame.add(self.RamachandranFavoredLabel)
        self.RamachandranFavoredValue = gtk.Label(
            self.QualityIndicators["RefinementRamachandranFavored"]
        )
        RamachandranFavoredBox_frame = gtk.Frame()
        self.RamachandranFavoredBox = gtk.EventBox()
        self.RamachandranFavoredBox.add(self.RamachandranFavoredValue)
        RamachandranFavoredBox_frame.add(self.RamachandranFavoredBox)

        rmsdBondsLabel_frame = gtk.Frame()
        self.rmsdBondsLabel = gtk.Label("rmsd(Bonds)")
        rmsdBondsLabel_frame.add(self.rmsdBondsLabel)
        self.rmsdBondsValue = gtk.Label(self.QualityIndicators["RefinementRmsdBonds"])
        rmsdBondsBox_frame = gtk.Frame()
        self.rmsdBondsBox = gtk.EventBox()
        self.rmsdBondsBox.add(self.rmsdBondsValue)
        rmsdBondsBox_frame.add(self.rmsdBondsBox)

        rmsdAnglesLabel_frame = gtk.Frame()
        self.rmsdAnglesLabel = gtk.Label("rmsd(Angles)")
        rmsdAnglesLabel_frame.add(self.rmsdAnglesLabel)
        self.rmsdAnglesValue = gtk.Label(self.QualityIndicators["RefinementRmsdAngles"])
        rmsdAnglesBox_frame = gtk.Frame()
        self.rmsdAnglesBox = gtk.EventBox()
        self.rmsdAnglesBox.add(self.rmsdAnglesValue)
        rmsdAnglesBox_frame.add(self.rmsdAnglesBox)

        MatrixWeightLabel_frame = gtk.Frame()
        self.MatrixWeightLabel = gtk.Label("Matrix Weight")
        MatrixWeightLabel_frame.add(self.MatrixWeightLabel)
        self.MatrixWeightValue = gtk.Label(
            self.QualityIndicators["RefinementMatrixWeight"]
        )
        MatrixWeightBox_frame = gtk.Frame()
        self.MatrixWeightBox = gtk.EventBox()
        self.MatrixWeightBox.add(self.MatrixWeightValue)
        MatrixWeightBox_frame.add(self.MatrixWeightBox)

        ligandIDLabel_frame = gtk.Frame()
        self.ligandIDLabel = gtk.Label("Ligand ID")
        ligandIDLabel_frame.add(self.ligandIDLabel)
        self.ligandIDValue = gtk.Label(self.spider_plot_data["PANDDA_site_ligand_id"])
        ligandIDBox_frame = gtk.Frame()
        self.ligandIDBox = gtk.EventBox()
        self.ligandIDBox.add(self.ligandIDValue)
        ligandIDBox_frame.add(self.ligandIDBox)

        ligand_occupancyLabel_frame = gtk.Frame()
        self.ligand_occupancyLabel = gtk.Label("occupancy")
        ligand_occupancyLabel_frame.add(self.ligand_occupancyLabel)
        self.ligand_occupancyValue = gtk.Label(
            self.spider_plot_data["PANDDA_site_occupancy"]
        )
        ligand_occupancyBox_frame = gtk.Frame()
        self.ligand_occupancyBox = gtk.EventBox()
        self.ligand_occupancyBox.add(self.ligand_occupancyValue)
        ligand_occupancyBox_frame.add(self.ligand_occupancyBox)

        ligand_BaverageLabel_frame = gtk.Frame()
        self.ligand_BaverageLabel = gtk.Label("B average")
        ligand_BaverageLabel_frame.add(self.ligand_BaverageLabel)
        self.ligand_BaverageValue = gtk.Label(
            self.spider_plot_data["PANDDA_site_B_average"]
        )
        ligand_BaverageBox_frame = gtk.Frame()
        self.ligand_BaverageBox = gtk.EventBox()
        self.ligand_BaverageBox.add(self.ligand_BaverageValue)
        ligand_BaverageBox_frame.add(self.ligand_BaverageBox)

        ligand_BratioSurroundingsLabel_frame = gtk.Frame()
        self.ligand_BratioSurroundingsLabel = gtk.Label("B ratio")
        ligand_BratioSurroundingsLabel_frame.add(self.ligand_BratioSurroundingsLabel)
        self.ligand_BratioSurroundingsValue = gtk.Label(
            self.spider_plot_data["PANDDA_site_B_ratio_residue_surroundings"]
        )
        ligand_BratioSurroundingsBox_frame = gtk.Frame()
        self.ligand_BratioSurroundingsBox = gtk.EventBox()
        self.ligand_BratioSurroundingsBox.add(self.ligand_BratioSurroundingsValue)
        ligand_BratioSurroundingsBox_frame.add(self.ligand_BratioSurroundingsBox)

        ligand_RSCCLabel_frame = gtk.Frame()
        self.ligand_RSCCLabel = gtk.Label("RSCC")
        ligand_RSCCLabel_frame.add(self.ligand_RSCCLabel)
        self.ligand_RSCCValue = gtk.Label(self.spider_plot_data["PANDDA_site_RSCC"])
        ligand_RSCCBox_frame = gtk.Frame()
        self.ligand_RSCCBox = gtk.EventBox()
        self.ligand_RSCCBox.add(self.ligand_RSCCValue)
        ligand_RSCCBox_frame.add(self.ligand_RSCCBox)

        ligand_rmsdLabel_frame = gtk.Frame()
        self.ligand_rmsdLabel = gtk.Label("rmsd")
        ligand_rmsdLabel_frame.add(self.ligand_rmsdLabel)
        self.ligand_rmsdValue = gtk.Label(self.spider_plot_data["PANDDA_site_rmsd"])
        ligand_rmsdBox_frame = gtk.Frame()
        self.ligand_rmsdBox = gtk.EventBox()
        self.ligand_rmsdBox.add(self.ligand_rmsdValue)
        ligand_rmsdBox_frame.add(self.ligand_rmsdBox)

        ligand_RSRLabel_frame = gtk.Frame()
        self.ligand_RSRLabel = gtk.Label("RSR")
        ligand_RSRLabel_frame.add(self.ligand_RSRLabel)
        self.ligand_RSRValue = gtk.Label(self.spider_plot_data["PANDDA_site_RSR"])
        ligand_RSRBox_frame = gtk.Frame()
        self.ligand_RSRBox = gtk.EventBox()
        self.ligand_RSRBox.add(self.ligand_RSRValue)
        ligand_RSRBox_frame.add(self.ligand_RSRBox)

        ligand_RSZDLabel_frame = gtk.Frame()
        self.ligand_RSZDLabel = gtk.Label("RSZD")
        ligand_RSZDLabel_frame.add(self.ligand_RSZDLabel)
        self.ligand_RSZDValue = gtk.Label(self.spider_plot_data["PANDDA_site_RSZD"])
        ligand_RSZDBox_frame = gtk.Frame()
        self.ligand_RSZDBox = gtk.EventBox()
        self.ligand_RSZDBox.add(self.ligand_RSZDValue)
        ligand_RSZDBox_frame.add(self.ligand_RSZDBox)

        outer_frame = gtk.Frame()
        hbox = gtk.HBox()

        frame = gtk.Frame()
        self.table_left = gtk.Table(8, 2, False)
        self.table_left.attach(RRfreeLabel_frame, 0, 1, 0, 1)
        self.table_left.attach(ResolutionLabel_frame, 0, 1, 1, 2)
        self.table_left.attach(MolprobityScoreLabel_frame, 0, 1, 2, 3)
        self.table_left.attach(RamachandranOutliersLabel_frame, 0, 1, 3, 4)
        self.table_left.attach(RamachandranFavoredLabel_frame, 0, 1, 4, 5)
        self.table_left.attach(rmsdBondsLabel_frame, 0, 1, 5, 6)
        self.table_left.attach(rmsdAnglesLabel_frame, 0, 1, 6, 7)
        self.table_left.attach(MatrixWeightLabel_frame, 0, 1, 7, 8)
        self.table_left.attach(RRfreeBox_frame, 1, 2, 0, 1)
        self.table_left.attach(ResolutionBox_frame, 1, 2, 1, 2)
        self.table_left.attach(MolprobityScoreBox_frame, 1, 2, 2, 3)
        self.table_left.attach(RamachandranOutliersBox_frame, 1, 2, 3, 4)
        self.table_left.attach(RamachandranFavoredBox_frame, 1, 2, 4, 5)
        self.table_left.attach(rmsdBondsBox_frame, 1, 2, 5, 6)
        self.table_left.attach(rmsdAnglesBox_frame, 1, 2, 6, 7)
        self.table_left.attach(MatrixWeightBox_frame, 1, 2, 7, 8)
        frame.add(self.table_left)
        hbox.add(frame)

        frame = gtk.Frame()
        self.table_right = gtk.Table(8, 2, False)
        self.table_right.attach(ligandIDLabel_frame, 0, 1, 0, 1)
        self.table_right.attach(ligand_occupancyLabel_frame, 0, 1, 1, 2)
        self.table_right.attach(ligand_BaverageLabel_frame, 0, 1, 2, 3)
        self.table_right.attach(ligand_BratioSurroundingsLabel_frame, 0, 1, 3, 4)
        self.table_right.attach(ligand_RSCCLabel_frame, 0, 1, 4, 5)
        self.table_right.attach(ligand_rmsdLabel_frame, 0, 1, 5, 6)
        self.table_right.attach(ligand_RSRLabel_frame, 0, 1, 6, 7)
        self.table_right.attach(ligand_RSZDLabel_frame, 0, 1, 7, 8)
        self.table_right.attach(ligandIDBox_frame, 1, 2, 0, 1)
        self.table_right.attach(ligand_occupancyBox_frame, 1, 2, 1, 2)
        self.table_right.attach(ligand_BaverageBox_frame, 1, 2, 2, 3)
        self.table_right.attach(ligand_BratioSurroundingsBox_frame, 1, 2, 3, 4)
        self.table_right.attach(ligand_RSCCBox_frame, 1, 2, 4, 5)
        self.table_right.attach(ligand_rmsdBox_frame, 1, 2, 5, 6)
        self.table_right.attach(ligand_RSRBox_frame, 1, 2, 6, 7)
        self.table_right.attach(ligand_RSZDBox_frame, 1, 2, 7, 8)
        frame.add(self.table_right)
        hbox.add(frame)

        outer_frame.add(hbox)
        self.vbox.add(outer_frame)

        button = gtk.Button(label="Show MolProbity to-do list")
        button.connect("clicked", self.show_molprobity_to_do)
        self.vbox.add(button)
        self.vbox.pack_start(frame)

        # SPACER
        self.vbox.add(gtk.Label(" "))

        ################################################################################
        # --- hbox for compound picture & spider_plot (formerly: refinement history) ---
        frame = gtk.Frame()
        self.hbox_for_info_graphics = gtk.HBox()

        # --- compound picture ---
        compound_frame = gtk.Frame()
        pic = gtk.gdk.pixbuf_new_from_file(
            os.path.join(
                os.getenv("XChemExplorer_DIR"),
                "xce",
                "image",
                "NO_COMPOUND_IMAGE_AVAILABLE.png",
            )
        )
        self.pic = pic.scale_simple(190, 190, gtk.gdk.INTERP_BILINEAR)
        self.image = gtk.Image()
        self.image.set_from_pixbuf(self.pic)
        compound_frame.add(self.image)
        self.hbox_for_info_graphics.add(compound_frame)

        # --- Spider Plot ---
        spider_plot_frame = gtk.Frame()
        spider_plot_pic = gtk.gdk.pixbuf_new_from_file(
            os.path.join(
                os.getenv("XChemExplorer_DIR"),
                "xce",
                "image",
                "NO_SPIDER_PLOT_AVAILABLE.png",
            )
        )
        self.spider_plot_pic = spider_plot_pic.scale_simple(
            190, 190, gtk.gdk.INTERP_BILINEAR
        )
        self.spider_plot_image = gtk.Image()
        self.spider_plot_image.set_from_pixbuf(self.spider_plot_pic)
        spider_plot_frame.add(self.spider_plot_image)
        self.hbox_for_info_graphics.add(spider_plot_frame)

        frame.add(self.hbox_for_info_graphics)
        self.vbox.add(frame)

        ################################################################################
        # --- pandda.inspect user comments ---
        outer_frame = gtk.Frame(label="pandda.inspect comments")
        vbox = gtk.VBox()
        ligand_site_name_label_frame = gtk.Frame()
        ligand_site_name_label = gtk.Label("Site Name")
        ligand_site_name_label_frame.add(ligand_site_name_label)
        ligand_inspect_confidence_label_frame = gtk.Frame()
        ligand_inspect_confidence_label = gtk.Label("Confidence")
        ligand_inspect_confidence_label_frame.add(ligand_inspect_confidence_label)
        ligand_inspect_interesting_label_frame = gtk.Frame()
        ligand_inspect_interesting_label = gtk.Label("Interesting")
        ligand_inspect_interesting_label_frame.add(ligand_inspect_interesting_label)
        ligand_inspect_comment_label_frame = gtk.Frame()
        ligand_inspect_comment_label = gtk.Label("Comment")
        ligand_inspect_comment_label_frame.add(ligand_inspect_comment_label)
        ligand_site_name_value_frame = gtk.Frame()
        self.ligand_site_name_value = gtk.Label("-")
        ligand_site_name_value_frame.add(self.ligand_site_name_value)
        ligand_inspect_confidence_value_frame = gtk.Frame()
        self.ligand_inspect_confidence_value = gtk.Label("-")
        ligand_inspect_confidence_value_frame.add(self.ligand_inspect_confidence_value)
        ligand_inspect_interesting_value_frame = gtk.Frame()
        self.ligand_inspect_interesting_value = gtk.Label("-")
        ligand_inspect_interesting_value_frame.add(
            self.ligand_inspect_interesting_value
        )
        ligand_inspect_comment_value_frame = gtk.Frame()
        self.ligand_inspect_comment_value = gtk.Label("-")
        ligand_inspect_comment_value_frame.add(self.ligand_inspect_comment_value)

        frame_pandda_inspect_comments_table = gtk.Frame()
        pandda_inspect_comments_table = gtk.Table(2, 6, False)
        pandda_inspect_comments_table.attach(ligand_site_name_label_frame, 0, 1, 0, 1)
        pandda_inspect_comments_table.attach(ligand_site_name_value_frame, 1, 2, 0, 1)
        pandda_inspect_comments_table.attach(
            ligand_inspect_confidence_label_frame, 2, 3, 0, 1
        )
        pandda_inspect_comments_table.attach(
            ligand_inspect_confidence_value_frame, 3, 4, 0, 1
        )
        pandda_inspect_comments_table.attach(
            ligand_inspect_interesting_label_frame, 4, 5, 0, 1
        )
        pandda_inspect_comments_table.attach(
            ligand_inspect_interesting_value_frame, 5, 6, 0, 1
        )
        pandda_inspect_comments_table.attach(
            ligand_inspect_comment_label_frame, 0, 1, 1, 2
        )
        pandda_inspect_comments_table.attach(
            ligand_inspect_comment_value_frame, 1, 6, 1, 2
        )

        frame_pandda_inspect_comments_table.add(pandda_inspect_comments_table)
        vbox.add(frame_pandda_inspect_comments_table)
        outer_frame.add(vbox)
        self.vbox.pack_start(outer_frame)

        #        # SPACER
        self.vbox.add(gtk.Label(" "))

        ################################################################################
        outer_frame = gtk.Frame(label="Sample Navigator")
        hboxSample = gtk.HBox()

        # --- crystal navigator combobox ---
        frame = gtk.Frame()
        self.vbox_sample_navigator = gtk.VBox()
        self.cb = gtk.combo_box_new_text()
        self.cb.connect("changed", self.ChooseXtal)
        self.vbox_sample_navigator.add(self.cb)
        # --- crystal navigator backward/forward button ---
        self.PREVbutton = gtk.Button(label="<<<")
        self.NEXTbutton = gtk.Button(label=">>>")
        self.PREVbutton.connect("clicked", self.ChangeXtal, -1)
        self.NEXTbutton.connect("clicked", self.ChangeXtal, +1)
        hbox = gtk.HBox()
        hbox.pack_start(self.PREVbutton)
        hbox.pack_start(self.NEXTbutton)
        self.vbox_sample_navigator.add(hbox)
        frame.add(self.vbox_sample_navigator)
        hboxSample.add(frame)

        # --- site navigator combobox ---
        frame = gtk.Frame()
        self.vbox_site_navigator = gtk.VBox()
        self.cb_site = gtk.combo_box_new_text()
        self.cb_site.connect("changed", self.ChooseSite)
        self.vbox_site_navigator.add(self.cb_site)
        # --- site navigator backward/forward button ---
        self.PREVbuttonSite = gtk.Button(label="<<<")
        self.NEXTbuttonSite = gtk.Button(label=">>>")
        self.PREVbuttonSite.connect("clicked", self.ChangeSite, -1)
        self.NEXTbuttonSite.connect("clicked", self.ChangeSite, +1)
        hbox = gtk.HBox()
        hbox.pack_start(self.PREVbuttonSite)
        hbox.pack_start(self.NEXTbuttonSite)
        self.vbox_site_navigator.add(hbox)
        frame.add(self.vbox_site_navigator)
        hboxSample.add(frame)

        outer_frame.add(hboxSample)
        self.vbox.add(outer_frame)

        # SPACER
        self.vbox.add(gtk.Label(" "))

        ################################################################################
        # --- current refinement stage ---
        outer_frame = gtk.Frame()
        hbox = gtk.HBox()

        frame = gtk.Frame(label="Analysis Status")
        vbox = gtk.VBox()
        self.experiment_stage_button_list = []
        for n, button in enumerate(self.experiment_stage):
            if n == 0:
                new_button = gtk.RadioButton(None, button[0])
            else:
                new_button = gtk.RadioButton(new_button, button[0])
            new_button.connect(
                "toggled", self.experiment_stage_button_clicked, button[1]
            )
            vbox.add(new_button)
            self.experiment_stage_button_list.append(new_button)
        frame.add(vbox)
        hbox.pack_start(frame)

        # --- ligand confidence ---
        frame = gtk.Frame(label="Ligand Confidence")
        vbox = gtk.VBox()
        self.ligand_confidence_button_list = []
        for n, criteria in enumerate(self.ligand_confidence_category):
            if n == 0:
                new_button = gtk.RadioButton(None, criteria)
            else:
                new_button = gtk.RadioButton(new_button, criteria)
            new_button.connect(
                "toggled", self.ligand_confidence_button_clicked, criteria
            )
            vbox.add(new_button)
            self.ligand_confidence_button_list.append(new_button)
        frame.add(vbox)
        hbox.pack_start(frame)

        outer_frame.add(hbox)
        self.vbox.pack_start(outer_frame)

        # SPACER
        self.vbox.add(gtk.Label(" "))

        # --- ligand modeling ---
        frame = gtk.Frame(label="Ligand Modeling")
        self.hbox_for_modeling = gtk.HBox()
        self.merge_ligand_button = gtk.Button(label="Merge Ligand")
        self.place_ligand_here_button = gtk.Button(label="Place Ligand here")
        self.hbox_for_modeling.add(self.place_ligand_here_button)
        self.place_ligand_here_button.connect("clicked", self.place_ligand_here)
        self.hbox_for_modeling.add(self.merge_ligand_button)
        self.merge_ligand_button.connect("clicked", self.merge_ligand_into_protein)
        frame.add(self.hbox_for_modeling)
        self.vbox.pack_start(frame)

        # --- refinement & options ---
        self.hbox_for_refinement = gtk.HBox()
        self.REFINEbutton = gtk.Button(label="Refine")
        self.RefinementParamsButton = gtk.Button(label="refinement parameters")
        self.REFINEbutton.connect("clicked", self.REFINE)
        self.hbox_for_refinement.add(self.REFINEbutton)
        self.RefinementParamsButton.connect("clicked", self.RefinementParams)
        self.hbox_for_refinement.add(self.RefinementParamsButton)
        self.vbox.add(self.hbox_for_refinement)

        # --- CANCEL button ---
        self.CANCELbutton = gtk.Button(label="CANCEL")
        self.CANCELbutton.connect("clicked", self.CANCEL)
        self.vbox.add(self.CANCELbutton)

        self.window.add(self.vbox)
        self.window.show_all()

    def CANCEL(self, widget):
        self.window.destroy()

    def ChangeXtal(self, widget, data=None):
        self.index = self.index + data
        if self.index < 0:
            self.index = 0
        if self.index >= len(self.Todo):
            self.index = 0
        self.cb.set_active(self.index)

    def ChooseXtal(self, widget):
        self.xtalID = str(widget.get_active_text())
        for n, item in enumerate(self.Todo):
            if str(item[0]) == self.xtalID:
                self.index = n
        self.merge_ligand_button.set_sensitive(True)
        self.place_ligand_here_button.set_sensitive(True)

        self.ligand_site_name_value.set_label("-")
        self.ligand_inspect_confidence_value.set_label("-")
        self.ligand_inspect_interesting_value.set_label("-")
        self.ligand_inspect_comment_value.set_label("-")

        self.refresh_site_combobox()
        self.db_dict_mainTable = {}
        self.db_dict_panddaTable = {}
        if str(self.Todo[self.index][0]) is not None:
            self.compoundID = str(self.Todo[self.index][1])
            self.refinement_outcome = str(self.Todo[self.index][5])
            self.update_RefinementOutcome_radiobutton()
        if (
            self.xtalID not in self.siteDict
        ):  # i.e. we are not working with a PanDDA model
            self.ligand_confidence = str(self.Todo[self.index][6])
            self.update_LigandConfidence_radiobutton()

        self.RefreshData()

    def update_RefinementOutcome_radiobutton(self):
        # updating dataset outcome radiobuttons
        current_stage = 0
        for i, entry in enumerate(self.experiment_stage):
            if entry[1].split()[0] == self.refinement_outcome.split()[0]:
                current_stage = i
                break
        for i, button in enumerate(self.experiment_stage_button_list):
            if i == current_stage:
                button.set_active(True)
                break

    def update_LigandConfidence_radiobutton(self):
        # updating ligand confidence radiobuttons
        current_stage = 0
        for i, entry in enumerate(self.ligand_confidence_category):
            print("--->", entry, self.ligand_confidence)
            try:
                if entry.split()[0] == self.ligand_confidence.split()[0]:
                    current_stage = i
                    break
            except IndexError:
                pass
        for i, button in enumerate(self.ligand_confidence_button_list):
            if i == current_stage:
                button.set_active(True)
                break

    def refresh_site_combobox(self):
        # reset self.pandda_index
        self.pandda_index = -1
        # clear CB first, 100 is sort of arbitrary since it's unlikely there will ever
        # be 100 sites
        for n in range(-1, 100):
            self.cb_site.remove_text(0)
        self.site_index = "0"
        self.event_index = "0"
        # only repopulate if site exists
        if self.xtalID in self.siteDict:
            for item in sorted(self.siteDict[self.xtalID]):
                self.cb_site.append_text("site: %s - event: %s" % (item[5], item[6]))

    def ChangeSite(self, widget, data=None):
        if self.xtalID in self.siteDict:
            self.pandda_index = self.pandda_index + data
            if self.pandda_index < 0:
                self.pandda_index = 0
            if self.pandda_index >= len(self.siteDict[self.xtalID]):
                self.pandda_index = 0
            self.cb_site.set_active(self.pandda_index)

    def ChooseSite(self, widget):
        tmp = str(widget.get_active_text())
        self.site_index = tmp.split()[1]
        self.event_index = tmp.split()[4]
        for n, item in enumerate(self.siteDict[self.xtalID]):
            if item[5] == self.site_index and item[6] == self.event_index:
                self.pandda_index = n
        self.RefreshSiteData()

    def RefreshSiteData(self):
        if self.pandda_index == -1:
            self.merge_ligand_button.set_sensitive(True)
            self.place_ligand_here_button.set_sensitive(True)
        else:
            self.merge_ligand_button.set_sensitive(False)
            self.place_ligand_here_button.set_sensitive(False)
            # and remove ligand molecule so that there is no temptation to merge it
            if len(coot_utils_XChem.molecule_number_list()) > 0:
                for imol in coot_utils_XChem.molecule_number_list():
                    if self.compoundID + ".pdb" in coot.molecule_name(imol):
                        coot.close_molecule(imol)

        print("pandda index", self.pandda_index)
        self.spider_plot = self.siteDict[self.xtalID][self.pandda_index][4]
        print("new spider plot:", self.spider_plot)
        self.event_map = self.siteDict[self.xtalID][self.pandda_index][0]
        print("new event map:", self.event_map)
        self.ligand_confidence = str(self.siteDict[self.xtalID][self.pandda_index][7])
        self.update_LigandConfidence_radiobutton()
        site_x = float(self.siteDict[self.xtalID][self.pandda_index][1])
        site_y = float(self.siteDict[self.xtalID][self.pandda_index][2])
        site_z = float(self.siteDict[self.xtalID][self.pandda_index][3])
        print("new site coordinates:", site_x, site_y, site_z)
        coot.set_rotation_centre(site_x, site_y, site_z)

        self.ligand_site_name_value.set_label(
            str(self.siteDict[self.xtalID][self.pandda_index][8])
        )
        self.ligand_inspect_confidence_value.set_label(
            str(self.siteDict[self.xtalID][self.pandda_index][9])
        )
        self.ligand_inspect_interesting_value.set_label(
            str(self.siteDict[self.xtalID][self.pandda_index][10])
        )
        self.ligand_inspect_comment_value.set_label(
            str(self.siteDict[self.xtalID][self.pandda_index][11])
        )

        self.spider_plot_data = (
            self.db.get_db_pandda_dict_for_sample_and_site_and_event(
                self.xtalID, self.site_index, self.event_index
            )
        )
        print(">>>>> spider plot data", self.spider_plot_data)
        self.ligandIDValue.set_label(self.spider_plot_data["PANDDA_site_ligand_id"])
        try:
            self.ligand_occupancyValue.set_label(
                str(round(float(self.spider_plot_data["PANDDA_site_occupancy"]), 2))
            )
        except ValueError:
            self.ligand_occupancyValue.set_label("-")

        try:
            self.ligand_BaverageValue.set_label(
                str(round(float(self.spider_plot_data["PANDDA_site_B_average"]), 2))
            )
        except ValueError:
            self.ligand_BaverageValue.set_label("-")

        try:
            self.ligand_BratioSurroundingsValue.set_label(
                str(
                    round(
                        float(
                            self.spider_plot_data[
                                "PANDDA_site_B_ratio_residue_surroundings"
                            ]
                        ),
                        2,
                    )
                )
            )
        except ValueError:
            self.ligand_BratioSurroundingsValue.set_label("-")

        try:
            self.ligand_RSCCValue.set_label(
                str(round(float(self.spider_plot_data["PANDDA_site_RSCC"]), 2))
            )
        except ValueError:
            self.ligand_RSCCValue.set_label("-")

        try:
            self.ligand_rmsdValue.set_label(
                str(round(float(self.spider_plot_data["PANDDA_site_rmsd"]), 2))
            )
        except ValueError:
            self.ligand_rmsdValue.set_label("-")

        try:
            self.ligand_RSRValue.set_label(
                str(round(float(self.spider_plot_data["PANDDA_site_RSR"]), 2))
            )
        except ValueError:
            self.ligand_RSRValue.set_label("-")

        try:
            self.ligand_RSZDValue.set_label(
                str(round(float(self.spider_plot_data["PANDDA_site_RSZD"]), 2))
            )
        except ValueError:
            self.ligand_RSZDValue.set_label("-")

        ################################################################################
        # delete old Event MAPs
        if len(coot_utils_XChem.molecule_number_list()) > 0:
            for imol in coot_utils_XChem.molecule_number_list():
                if "map.native.ccp4" in coot.molecule_name(imol):
                    coot.close_molecule(imol)

        ################################################################################
        # Spider plot
        # Note: refinement history was shown instead previously
        if os.path.isfile(self.spider_plot):
            spider_plot_pic = gtk.gdk.pixbuf_new_from_file(self.spider_plot)
        else:
            spider_plot_pic = gtk.gdk.pixbuf_new_from_file(
                os.path.join(
                    os.getenv("XChemExplorer_DIR"),
                    "xce",
                    "image",
                    "NO_SPIDER_PLOT_AVAILABLE.png",
                )
            )
        self.spider_plot_pic = spider_plot_pic.scale_simple(
            190, 190, gtk.gdk.INTERP_BILINEAR
        )
        self.spider_plot_image.set_from_pixbuf(self.spider_plot_pic)

        ################################################################################
        # check for PANDDAs EVENT maps
        if os.path.isfile(self.event_map):
            coot.handle_read_ccp4_map((self.event_map), 0)
            for imol in coot_utils_XChem.molecule_number_list():
                if self.event_map in coot.molecule_name(imol):
                    coot.set_contour_level_in_sigma(imol, 2)
                    coot.set_last_map_colour(0.74, 0.44, 0.02)

    def experiment_stage_button_clicked(self, widget, data=None):
        self.db_dict_mainTable["RefinementOutcome"] = data
        print(
            "==> XCE: setting Refinement Outcome for "
            + self.xtalID
            + " to "
            + str(data)
            + " in mainTable of datasource"
        )
        self.db.update_data_source(self.xtalID, self.db_dict_mainTable)

    def ligand_confidence_button_clicked(self, widget, data=None):
        print("PANDDA_index", self.pandda_index)
        if self.pandda_index == -1:
            self.db_dict_mainTable["RefinementLigandConfidence"] = data
            print(
                "==> XCE: setting Ligand Confidence for "
                + self.xtalID
                + " to "
                + str(data)
                + " in mainTable of datasource"
            )
            self.db.update_data_source(self.xtalID, self.db_dict_mainTable)
            self.Todo[self.index][6] = data
        else:
            self.db_dict_panddaTable["PANDDA_site_confidence"] = data
            print(
                "==> XCE: setting Ligand Confidence for "
                + self.xtalID
                + " (site="
                + str(self.site_index)
                + ", event="
                + str(self.event_index)
                + ") to "
                + str(data)
                + " in panddaTable of datasource"
            )
            self.db.update_site_event_panddaTable(
                self.xtalID, self.site_index, self.event_index, self.db_dict_panddaTable
            )
            self.siteDict[self.xtalID][self.pandda_index][7] = data

    def RefreshData(self):
        # reset spider plot image
        spider_plot_pic = gtk.gdk.pixbuf_new_from_file(
            os.path.join(
                os.getenv("XChemExplorer_DIR"),
                "xce",
                "image",
                "NO_SPIDER_PLOT_AVAILABLE.png",
            )
        )
        self.spider_plot_pic = spider_plot_pic.scale_simple(
            190, 190, gtk.gdk.INTERP_BILINEAR
        )
        self.spider_plot_image.set_from_pixbuf(self.spider_plot_pic)

        # initialize Refinement library
        self.Refine = XChemRefine.RefineOld(
            self.project_directory, self.xtalID, self.compoundID, self.data_source
        )
        self.Serial = self.Refine.GetSerial()
        if self.Serial == 1:
            # i.e. no refinement has been done; data is probably straight out of dimple
            if os.path.isfile(
                os.path.join(self.project_directory, self.xtalID, self.pdb_style)
            ):
                print(
                    "==> XCE: updating quality indicators in data source for "
                    + self.xtalID
                )
                XChemUtils.parse().update_datasource_with_PDBheader(
                    self.xtalID,
                    self.data_source,
                    os.path.join(self.project_directory, self.xtalID, self.pdb_style),
                )
                XChemUtils.parse().update_datasource_with_phenix_validation_summary(
                    self.xtalID, self.data_source, ""
                )  # '' because file does not exist
            elif os.path.isfile(
                os.path.join(self.project_directory, self.xtalID, "dimple.pdb")
            ):
                print(
                    "==> XCE: updating quality indicators in data source for "
                    + self.xtalID
                )
                XChemUtils.parse().update_datasource_with_PDBheader(
                    self.xtalID,
                    self.data_source,
                    os.path.join(self.project_directory, self.xtalID, "dimple.pdb"),
                )
                XChemUtils.parse().update_datasource_with_phenix_validation_summary(
                    self.xtalID, self.data_source, ""
                )  # '' because file does not exist

        # all this information is now updated in the datasource after each refinement
        # cycle
        self.QualityIndicators = self.db.get_db_dict_for_sample(self.xtalID)

        ################################################################################
        # history
        # if the structure was previously refined, try to read the parameters
        if self.Serial > 1:
            self.RefmacParams = self.Refine.ParamsFromPreviousCycle(self.Serial - 1)

        ################################################################################
        # update pdb & maps

        ################################################################################
        # delete old PDB and MAP files
        # - get a list of all molecules which are currently opened in COOT
        # - remove all molecules/ maps before loading a new set
        if len(coot_utils_XChem.molecule_number_list()) > 0:
            for item in coot_utils_XChem.molecule_number_list():
                coot.close_molecule(item)

        ################################################################################
        # read new PDB files
        # read protein molecule after ligand so that this one is the active molecule
        coot.set_nomenclature_errors_on_read("ignore")
        if os.path.isfile(
            os.path.join(self.project_directory, self.xtalID, self.compoundID + ".pdb")
        ):
            imol = coot.handle_read_draw_molecule_with_recentre(
                os.path.join(
                    self.project_directory, self.xtalID, self.compoundID + ".pdb"
                ),
                0,
            )
            self.mol_dict["ligand"] = imol
            coot.read_cif_dictionary(
                os.path.join(
                    self.project_directory, self.xtalID, self.compoundID + ".cif"
                )
            )
        if not os.path.isfile(
            os.path.join(self.project_directory, self.xtalID, self.pdb_style)
        ):
            os.chdir(os.path.join(self.project_directory, self.xtalID))
        if os.path.isfile(
            os.path.join(self.project_directory, self.xtalID, self.pdb_style)
        ):
            os.chdir(os.path.join(self.project_directory, self.xtalID))
            imol = coot.handle_read_draw_molecule_with_recentre(
                os.path.join(self.project_directory, self.xtalID, self.pdb_style), 0
            )
        elif os.path.isfile(
            os.path.join(self.project_directory, self.xtalID, "dimple.pdb")
        ):
            os.chdir(os.path.join(self.project_directory, self.xtalID))
            imol = coot.handle_read_draw_molecule_with_recentre(
                os.path.join(self.project_directory, self.xtalID, "dimple.pdb"), 0
            )
        else:
            self.go_to_next_xtal()
        self.mol_dict["protein"] = imol
        for item in coot_utils_XChem.molecule_number_list():
            if coot.molecule_name(item).endswith(self.pdb_style):
                # master switch to show symmetry molecules
                coot.set_show_symmetry_master(1)
                coot.set_show_symmetry_molecule(item, 1)  # show symm for model

        ################################################################################
        # read fofo maps
        # - read ccp4 map: 0 - 2fofc map, 1 - fofc.map
        # read 2fofc map last so that one can change its contour level
        if os.path.isfile(
            os.path.join(self.project_directory, self.xtalID, "2fofc.map")
        ):
            coot.set_default_initial_contour_level_for_difference_map(3)
            coot.handle_read_ccp4_map(
                os.path.join(self.project_directory, self.xtalID, "fofc.map"), 1
            )
            coot.set_default_initial_contour_level_for_map(1)
            coot.handle_read_ccp4_map(
                os.path.join(self.project_directory, self.xtalID, "2fofc.map"), 0
            )
            coot.set_last_map_colour(0, 0, 1)
        else:
            # try to open mtz file with same name as pdb file
            coot.set_default_initial_contour_level_for_map(1)
            if os.path.isfile(
                os.path.join(self.project_directory, self.xtalID, self.mtz_style)
            ):
                coot.auto_read_make_and_draw_maps(
                    os.path.join(self.project_directory, self.xtalID, self.mtz_style)
                )
            elif os.path.isfile(
                os.path.join(self.project_directory, self.xtalID, "dimple.mtz")
            ):
                coot.auto_read_make_and_draw_maps(
                    os.path.join(self.project_directory, self.xtalID, "dimple.mtz")
                )

        ################################################################################
        # update Quality Indicator table
        try:
            self.RRfreeValue.set_label(
                str(round(float(self.QualityIndicators["RefinementRcryst"]), 3))
                + " / "
                + str(round(float(self.QualityIndicators["RefinementRfree"]), 3))
            )
        except ValueError:
            self.RRfreeValue.set_label("-")

        try:
            self.RRfreeBox.modify_bg(
                gtk.STATE_NORMAL,
                gtk.gdk.color_parse(
                    self.QualityIndicators["RefinementRfreeTraficLight"]
                ),
            )
        except ValueError:
            pass
        self.ResolutionValue.set_label(self.QualityIndicators["RefinementResolution"])
        try:
            self.ResolutionBox.modify_bg(
                gtk.STATE_NORMAL,
                gtk.gdk.color_parse(self.QualityIndicators["RefinementResolutionTL"]),
            )
        except ValueError:
            pass
        self.MolprobityScoreValue.set_label(
            self.QualityIndicators["RefinementMolProbityScore"]
        )
        try:
            self.MolprobityScoreBox.modify_bg(
                gtk.STATE_NORMAL,
                gtk.gdk.color_parse(
                    self.QualityIndicators["RefinementMolProbityScoreTL"]
                ),
            )
        except ValueError:
            pass
        self.RamachandranOutliersValue.set_label(
            self.QualityIndicators["RefinementRamachandranOutliers"]
        )
        try:
            self.RamachandranOutliersBox.modify_bg(
                gtk.STATE_NORMAL,
                gtk.gdk.color_parse(
                    self.QualityIndicators["RefinementRamachandranOutliersTL"]
                ),
            )
        except ValueError:
            pass
        self.RamachandranFavoredValue.set_label(
            self.QualityIndicators["RefinementRamachandranFavored"]
        )
        try:
            self.RamachandranFavoredBox.modify_bg(
                gtk.STATE_NORMAL,
                gtk.gdk.color_parse(
                    self.QualityIndicators["RefinementRamachandranFavoredTL"]
                ),
            )
        except ValueError:
            pass
        self.rmsdBondsValue.set_label(self.QualityIndicators["RefinementRmsdBonds"])
        try:
            self.rmsdBondsBox.modify_bg(
                gtk.STATE_NORMAL,
                gtk.gdk.color_parse(self.QualityIndicators["RefinementRmsdBondsTL"]),
            )
        except ValueError:
            pass
        self.rmsdAnglesValue.set_label(self.QualityIndicators["RefinementRmsdAngles"])
        try:
            self.rmsdAnglesBox.modify_bg(
                gtk.STATE_NORMAL,
                gtk.gdk.color_parse(self.QualityIndicators["RefinementRmsdAnglesTL"]),
            )
        except ValueError:
            pass
        self.MatrixWeightValue.set_label(
            self.QualityIndicators["RefinementMatrixWeight"]
        )

        try:
            pic = gtk.gdk.pixbuf_new_from_file(
                os.path.join(
                    self.project_directory, self.xtalID, self.compoundID + ".png"
                )
            )
        except gobject.GError:
            pic = gtk.gdk.pixbuf_new_from_file(
                os.path.join(
                    os.getenv("XChemExplorer_DIR"),
                    "xce",
                    "image",
                    "NO_COMPOUND_IMAGE_AVAILABLE.png",
                )
            )
        self.pic = pic.scale_simple(190, 190, gtk.gdk.INTERP_BILINEAR)
        self.image.set_from_pixbuf(self.pic)

    def go_to_next_xtal(self):
        self.index += 1
        if self.index >= len(self.Todo):
            self.index = len(self.Todo)
        self.cb.set_active(self.index)

    def REFINE(self, widget):
        #######################################################
        # create folder for new refinement cycle
        os.mkdir(
            os.path.join(
                self.project_directory, self.xtalID, "Refine_" + str(self.Serial)
            )
        )

        #######################################################
        # write PDB file
        # now take protein pdb file and write it to newly create Refine_<serial> folder
        # note: the user has to make sure that the ligand file was merged into main file
        for item in coot_utils_XChem.molecule_number_list():
            if coot.molecule_name(item).endswith(self.pdb_style):
                coot.write_pdb_file(
                    item,
                    os.path.join(
                        self.project_directory,
                        self.xtalID,
                        "Refine_" + str(self.Serial),
                        "in.pdb",
                    ),
                )
                break
            elif coot.molecule_name(item).endswith("dimple.pdb"):
                coot.write_pdb_file(
                    item,
                    os.path.join(
                        self.project_directory,
                        self.xtalID,
                        "Refine_" + str(self.Serial),
                        "in.pdb",
                    ),
                )
                break

        #######################################################
        # run REFMAC
        self.Refine.RunRefmac(
            self.Serial, self.RefmacParams, self.external_software, self.xce_logfile
        )

        self.index += 1
        if self.index >= len(self.Todo):
            self.index = 0
        self.cb.set_active(self.index)

    def RefinementParams(self, widget):
        print("\n==> XCE: changing refinement parameters")
        self.RefmacParams = self.Refine.RefinementParams(self.RefmacParams)

    def set_selection_mode(self, widget):
        self.selection_mode = widget.get_active_text()

    def get_samples_to_look_at(self, widget):
        if self.selection_mode == "":
            self.status_label.set_text("select model stage")
            return
        self.status_label.set_text("checking datasource for samples... ")
        # first remove old samples if present
        if len(self.Todo) != 0:
            for n, item in enumerate(self.Todo):
                self.cb.remove_text(0)
        self.Todo = []
        self.siteDict = {}
        self.Todo, self.siteDict = self.db.get_todoList_for_coot(self.selection_mode)
        self.status_label.set_text("found %s samples" % len(self.Todo))
        # refresh sample CB
        for item in sorted(self.Todo):
            self.cb.append_text("%s" % item[0])
        if self.siteDict == {}:
            self.cb_site.set_sensitive(False)
            self.PREVbuttonSite.set_sensitive(False)
            self.NEXTbuttonSite.set_sensitive(False)
        else:
            self.cb_site.set_sensitive(True)
            self.PREVbuttonSite.set_sensitive(True)
            self.NEXTbuttonSite.set_sensitive(True)

    def update_plot(self, refinement_cycle, Rfree, Rcryst):
        fig = Figure(figsize=(2, 2), dpi=50)
        Plot = fig.add_subplot(111)
        Plot.set_ylim([0, max(Rcryst + Rfree)])
        Plot.set_xlabel("Refinement Cycle", fontsize=12)
        Plot.plot(refinement_cycle, Rfree, label="Rfree", linewidth=2)
        Plot.plot(refinement_cycle, Rcryst, label="Rcryst", linewidth=2)
        Plot.legend(
            bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
            loc=3,
            ncol=2,
            mode="expand",
            borderaxespad=0.0,
            fontsize=12,
        )
        return fig

    def place_ligand_here(self, widget):
        print("===> XCE: moving ligand to pointer")
        print("LIGAND: ", self.mol_dict["ligand"])
        coot_utils_XChem.move_molecule_here(self.mol_dict["ligand"])

    def merge_ligand_into_protein(self, widget):
        print("===> XCE: merge ligand into protein structure")
        coot.merge_molecules_py([self.mol_dict["ligand"]], self.mol_dict["protein"])
        print("===> XCE: deleting ligand molecule")
        coot.close_molecule(self.mol_dict["ligand"])

    def show_molprobity_to_do(self, widget):
        if os.path.isfile(
            os.path.join(
                self.project_directory,
                self.xtalID,
                "Refine_" + str(self.Serial - 1),
                "molprobity_coot.py",
            )
        ):
            print("==> XCE: running MolProbity Summary for", self.xtalID)
            coot.run_script(
                os.path.join(
                    self.project_directory,
                    self.xtalID,
                    "Refine_" + str(self.Serial - 1),
                    "molprobity_coot.py",
                )
            )
        else:
            print(
                "==> XCE: cannot find "
                + os.path.join(
                    self.project_directory,
                    self.xtalID,
                    "Refine_" + str(self.Serial - 1),
                    "molprobity_coot.py",
                )
            )


if __name__ == "__main__":
    GUI().StartGUI()
