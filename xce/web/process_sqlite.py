#!/usr/local/anaconda/sgc_default/envs/sgc_default/bin/python

# Author: Brian Marsden
#
# handover date: 16/12/2016

from rdkit.Chem import Draw
from rdkit import Chem
import shutil
import os
import sys
import getopt
import sqlite3
import csv
import zipfile
import glob


targetID = ""
panddadir = ""


def writeTableRow(row, htmlfile):
    "Write a row of data to the HTML file"
    htmlfile.write("<tr>")
    htmlfile.write("<td>" + row["ModelName"] + "</td>\n")
    htmlfile.write("<td>" + row["CompoundSMILES"] + "</td>\n")
    htmlfile.write(
        "<td><img src='compoundImages/" + row["CompoundCode"] + ".png'></td>\n"
    )
    htmlfile.write("<td>" + row["PANDDA_site_name"] + "</td>\n")
    htmlfile.write("<td>" + row["LigandConfidence"] + "</td>\n")
    htmlfile.write("<td>" + row["ModelStatus"] + "</td>\n")
    htmlfile.write(
        "<td><img src='residueplots/" + row["ModelName"] + ".png' width=150px></td>\n"
    )
    htmlfile.write(
        "<td><div id='"
        + row["ModelName"]
        + "'><a href='icbs/"
        + row["ModelName"]
        + ".html'><img src='mapImages/"
        + row["ModelName"]
        + "_"
        + row["CompoundCode"]
        + "_small.png'></a></div></td>\n"
    )
    #  htmlfile.write("<td><div id='"+row['ModelName']+"'><a href='icbs/"+row['ModelName']+".html'><img src='mapImages/"+row['ModelName']+"_small.png'></a></div></td>\n")
    htmlfile.write("<td>" + row["PANDDA_site_comment"] + "</td>\n")
    #  htmlfile.write("<td>TBD</td>\n")
    #  htmlfile.write("<td>http://www.rcsb.org/pdb/explore/explore.do?structureId=%s</td>\n" %row['Deposition_PDB_ID'])
    htmlfile.write(
        '<td><a target="_blank" href="http://www.rcsb.org/pdb/explore/explore.do?structureId={0!s}">{1!s}</a></td>\n'.format(
            row["Deposition_PDB_ID"], row["Deposition_PDB_ID"]
        )
    )
    htmlfile.write("<td>" + row["DataProcessingResolutionHigh"] + "</td>\n")
    htmlfile.write("<td>" + row["DataProcessingSpaceGroup"] + "</td>\n")
    htmlfile.write("<td>" + row["DataProcessingUnitCell"] + "</td>\n")
    htmlfile.write("<td><a href='pdbs/" + row["ModelName"] + ".pdb'>Save</a></td>\n")
    htmlfile.write("<td><a href='maps/" + row["ModelName"] + ".mtz'>Save</a></td>\n")
    htmlfile.write("<td><a href='maps/" + row["ModelName"] + ".ccp4'>Save</a></td>\n")
    htmlfile.write("</tr>\n")
    return


def writeICBPage(row, panddadir):
    "Write HTML file for ICB visualisation"
    icbhtmlfile = open(panddadir + "/icbs/" + row["ModelName"] + ".html", "w")
    icbhtmlfile.write("<html>\n")
    icbhtmlfile.write("<head>\n")
    icbhtmlfile.write(
        '<meta http-equiv="Content-type" content="text/html; charset=utf-8">\n'
    )
    icbhtmlfile.write(
        '<meta name="viewport" content="width=device-width,initial-scale=1">\n'
    )
    icbhtmlfile.write('<script src="http://molsoft.com/lib/acticm.js"> </script>\n')
    icbhtmlfile.write("<title>" + row["ModelName"] + "</title>\n")
    icbhtmlfile.write("</head>\n")
    icbhtmlfile.write("<body>\n")
    icbhtmlfile.write(
        "<h2>" + row["CrystalName"] + " " + row["CompoundCode"] + " event</h2>\n"
    )
    icbhtmlfile.write(
        '<div id="wait"><h3 style="color:red;">Please wait whilst the interactive viewer is loaded!</h3></div>\n'
    )
    icbhtmlfile.write(
        '<div id="con" style="width: 800px; height: 600px; border: 2px solid #ABABAB">\n'
    )
    icbhtmlfile.write(
        '<img id="pdbloader" src="../mapImages/'
        + row["ModelName"]
        + "_"
        + row["CompoundCode"]
        + '_large.png" />\n'
    )
    #  icbhtmlfile.write('<img id="pdbloader" src="../mapImages/'+row['ModelName']+'_large.png" />\n')
    icbhtmlfile.write("</div>\n")
    icbhtmlfile.write('<div id="details">\n')
    icbhtmlfile.write("<table>\n")
    icbhtmlfile.write("<tr><td>Protein: " + row["CrystalName"] + "</td></tr>\n")
    icbhtmlfile.write(
        "<tr><td>Compound ID: "
        + row["CompoundCode"]
        + "<br>Compound SMILES: "
        + row["CompoundSMILES"]
        + "</td>\n"
    )
    icbhtmlfile.write(
        "<td><img src='../compoundImages/" + row["CompoundCode"] + ".png'></td></tr>\n"
    )
    icbhtmlfile.write("</table>\n")
    icbhtmlfile.write("</div>\n")
    icbhtmlfile.write("<script>\n")
    icbhtmlfile.write("function onLoadActiveIcm()\n")
    icbhtmlfile.write("{\n")
    icbhtmlfile.write('  act = new ActiveIcmJS("con");\n')
    icbhtmlfile.write(
        '  act.projectFile = "'
        + row["ModelName"]
        + "_"
        + row["CompoundCode"]
        + '.icb";\n'
    )
    #  icbhtmlfile.write('  act.projectFile = "'+row['ModelName']+'_'+row['CompoundCode']+'.icb.gz";\n')
    icbhtmlfile.write("  act.searchBarVisible = false;\n")
    icbhtmlfile.write("  act.sequenceViewVisibleAuto = false;\n")
    icbhtmlfile.write("  act.tableViewVisibleAuto = false;\n")
    icbhtmlfile.write(
        '  document.getElementById("wait").style.visibility = "hidden";\n'
    )
    icbhtmlfile.write("}\n")
    icbhtmlfile.write("</script>\n")
    icbhtmlfile.write("</body>\n")
    icbhtmlfile.write("</html>\n")
    icbhtmlfile.close()
    return


def main(argv):
    sqlitefile = ""
    # Process command line options
    try:
        opts, args = getopt.getopt(
            argv, "t:s:d:", ["targetID=", "sqlitefile=", "panddadir="]
        )
    except getopt.GetoptError as err:
        print(err)
        print("process.py -t <TargetID> -s <SQLiteFile> -d <PANDDA dir>")
        sys.exit(2)
    if len(opts) < 3:
        print("Missing arguments:")
        print("process.py -t <TargetID> -s <SQLiteFile> -d <PANDDA dir>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("process.py -t <TargetID> -s <SQLiteFile>")
            sys.exit()
        elif opt in ("-t", "--targetID"):
            targetID = arg
        elif opt in ("-s", "--sqlitefile"):
            sqlitefile = arg
        elif opt in ("-d", "--panddadir"):
            panddadir = arg

    # Create directory structure
    if not os.path.exists(panddadir + "/compoundImages"):
        os.makedirs(panddadir + "/compoundImages")
    if not os.path.exists(panddadir + "/icbs"):
        os.makedirs(panddadir + "/icbs")
    if not os.path.exists(panddadir + "/pdbs"):
        os.makedirs(panddadir + "/pdbs")
    if not os.path.exists(panddadir + "/maps"):
        os.makedirs(panddadir + "/maps")
    if not os.path.exists(panddadir + "/residueplots"):
        os.makedirs(panddadir + "/residueplots")
    if not os.path.exists(panddadir + "/mapImages"):
        os.makedirs(panddadir + "/mapImages")

    # Create HTML file and write header
    htmlfile = open(panddadir + "/index.html", "w")
    htmlfile.write("<html>\n")
    htmlfile.write("<head>\n")
    htmlfile.write(
        '<link rel="stylesheet" type="text/css" href="css/jquery.dataTables.min.css">\n'
    )
    htmlfile.write(
        '<meta http-equiv="Content-type" content="text/html; charset=utf-8">\n'
    )
    htmlfile.write(
        '<meta name="viewport" content="width=device-width,initial-scale=1">\n'
    )
    htmlfile.write("<title>" + targetID + " Fragment Hits</title>\n")
    htmlfile.write(
        '<script type="text/javascript" language="javascript" src="js/jquery-1.12.3.min.js">\n'
    )
    htmlfile.write("</script>\n")
    htmlfile.write(
        '<script type="text/javascript" language="javascript" src="js/jquery.dataTables.min.js">\n'
    )
    htmlfile.write("</script>\n")
    htmlfile.write('<script type="text/javascript" class="init">\n')
    htmlfile.write("$(document).ready(function() {\n")
    htmlfile.write("$('#example').DataTable( {\n")
    htmlfile.write("'bautoWidth': false,\n")
    htmlfile.write("'columns': [\n")
    htmlfile.write("{ 'width': '6%' },\n")
    htmlfile.write("{ 'width': '6%' },\n")
    htmlfile.write("{ 'width': '7%' },\n")
    htmlfile.write("{ 'width': '8%' },\n")
    htmlfile.write("{ 'width': '6%' },\n")
    htmlfile.write("{ 'width': '6%' },\n")
    htmlfile.write("{ 'width': '9%' },\n")
    htmlfile.write("{ 'width': '9%' },\n")
    htmlfile.write("{ 'width': '12%' },\n")
    htmlfile.write("{ 'width': '5%' },\n")
    htmlfile.write("{ 'width': '4%' },\n")
    htmlfile.write("{ 'width': '6%' },\n")
    htmlfile.write("{ 'width': '6%' },\n")
    htmlfile.write("{ 'width': '3%' },\n")
    htmlfile.write("{ 'width': '4%' },\n")
    htmlfile.write("{ 'width': '3%' }\n")
    htmlfile.write("]\n")
    htmlfile.write("} )\n")
    htmlfile.write("} );\n")
    htmlfile.write("</script>\n")
    htmlfile.write("</head>\n")
    htmlfile.write("<body>\n")
    htmlfile.write("<H3>Ligand-bound models for " + targetID + "</h3>")
    htmlfile.write(
        """<h4>Interpreting 'Ligand confidence'</h4>
<p><u>4 - High Confidence:</u>  The expected ligand was easily interpretable from clear density, and subsequent refinement was well-behaved.  This ligand can be trusted.
<br><u>3 - Clear density, unexpected ligand:</u>  Density very clearly showed a well-defined ligand, but that ligand was unexpected in that crystal/dataset.  The observed ligand was modelled anyway, because its presence could be explained in some way.
<br><u>2 - Correct ligand, weak density:</u>  Though density was weak, it was possible to model the expected ligand, possibly including other circumstantial evidence (e.g. similar ligand in another model).
<br><u>1 - Low Confidence:</u>  The ligand model is to be treated with scepticism, because the evidence (density, identity, pose) were not convincing.
<h4>Interpreting 'Model status':</h4>
<p><u>6 - Deposited:</u>  The model has been deposited in the PDB.
<br><u>5 - Deposition ready:</u>  The model is fully error-free, in every residue, and is ready for deposition.
<br><u>4 - CompChem ready:</u>  The model is complete and correct in the region of the bound ligand.  There may be remaining small errors elsewhere in the structure, but they are far away and unlikely to be relevant to any computational analysis or compound design.
<h4>Interpreting 'Ligand validation' spider plots:</h4>  Each axis represents one of the values described below; small is better, and large values on any axis implies that further investigation is warranted.
<p><u>Quality (RSCC)</u> reflects the fit of the atoms to the experimental density, and should typically be greater than 0.7.
<br><u>Accuracy (RSZD)</u> measures the amount of difference density that is found around these atoms, and should be below 3.
<br><u>B-factor ratio</u> measures the consistency of the model with surrounding protein, and is calculated from the B factors of respectively the changed atoms and all side-chain atoms within 4&#8491;.  Large values (>3) reflect poor evidence for the model, and intermediate values (1.5+) indicate errors in refinement or modelling; for weakly-binding ligands, systematically large ratios may be justifiable.
<br><u>RMSD</u> compares the positions of all atoms built into event density, with their positions after final refinement, and should be below 1&#8491;.
<br><u>Precision (RSZO/OCC)</u> measures how clear the density is after refinement.  (This is not a quality indicator, but is related to strength of binding but not in a straightforward way.)
<p></p>\n"""
    )
    htmlfile.write("<h4>Download data</h4>\n")
    htmlfile.write("<ul>\n")
    htmlfile.write(
        "<li><a href='pdbs/allPDBs.zip'>Download all PDB model files<a></li>\n"
    )
    htmlfile.write(
        "<li><a href='maps/allEventMaps.zip'>Download all Event Map files<a></li>\n"
    )
    htmlfile.write("</ul>")
    htmlfile.write('<table id="example" class="display" cellspacing="0">\n')
    htmlfile.write("<thead>\n")
    htmlfile.write("<tr>\n")
    htmlfile.write("<th>Model Name</th>\n")
    htmlfile.write("<th>Compound SMILES</th>\n")
    htmlfile.write("<th>Compound Structure</th>\n")
    htmlfile.write("<th>Site Name</th>\n")
    htmlfile.write("<th>Ligand Confidence</th>\n")
    htmlfile.write("<th>Model Status</th>\n")
    htmlfile.write("<th>Ligand Validation</th>\n")
    htmlfile.write("<th>Event Map 3D</th>\n")
    htmlfile.write("<th>Comment</th>\n")
    htmlfile.write("<th>PDB Identifier</th>\n")
    htmlfile.write("<th>Resol</th>\n")
    htmlfile.write("<th>Spacegroup</th>\n")
    htmlfile.write("<th>Cell</th>\n")
    htmlfile.write("<th>PDB</th>\n")
    htmlfile.write("<th>MTZ</th>\n")
    htmlfile.write("<th>Event Map</th>\n")
    htmlfile.write("</tr>\n")
    htmlfile.write("</thead>\n")
    htmlfile.write("<tfoot>\n")
    htmlfile.write("<tr>\n")
    htmlfile.write("<th>Model Name</th>\n")
    htmlfile.write("<th>Compound SMILES</th>\n")
    htmlfile.write("<th>Compound Structure</th>\n")
    htmlfile.write("<th>Site Name</th>\n")
    htmlfile.write("<th>Ligand Confidence</th>\n")
    htmlfile.write("<th>Model Status</th>\n")
    htmlfile.write("<th>Ligand Validation</th>\n")
    htmlfile.write("<th>Event Map 3D</th>\n")
    htmlfile.write("<th>Comment</th>\n")
    htmlfile.write("<th>PDB Identifier</th>\n")
    htmlfile.write("<th>Resol</th>\n")
    htmlfile.write("<th>Spacegroup</th>\n")
    htmlfile.write("<th>Cell</th>\n")
    htmlfile.write("<th>PDB</th>\n")
    htmlfile.write("<th>MTZ</th>\n")
    htmlfile.write("<th>Event Map</th>\n")
    htmlfile.write("</tr>\n")
    htmlfile.write("</tfoot>\n")
    htmlfile.write("<tbody>\n")

    # Now walk through the input data
    with open("foricm.csv", "wb") as f:
        with sqlite3.connect(sqlitefile) as c:
            c.row_factory = sqlite3.Row
            cur = c.cursor()

            #      sql = ( "select p.ID,p.CrystalName,p.PANDDA_site_event_index,p.CrystalName || '_event'|| p.PANDDA_site_event_index "
            #              " as ModelName,m.CompoundCode,m.CompoundSMILES,p.PANDDA_site_name,p.PANDDA_site_confidence "
            #              " as LigandConfidence,p.RefinementOutcome "
            #              " as ModelStatus,p.PANDDA_site_comment,p.PANDDA_site_x,p.PANDDA_site_y,p.PANDDA_site_z, "
            #              "                p.PANDDA_site_spider_plot,m.DataProcessingResolutionHigh,m.DataProcessingSpaceGroup,"
            #              "                m.DataProcessingUnitCell,m.RefinementPDB_latest,m.RefinementMTZ_latest,p.PANDDA_site_event_map "
            #              " from panddaTable as p, mainTable as m "
            #              " where p.CrystalName=m.CrystalName and p.PANDDA_site_ligand_placed='True' and "
            #              "       (LigandConfidence like '1%' or LigandConfidence like '2%' or LigandConfidence like '3%' or LigandConfidence like '4%') "
            #              " order by p.CrystalName,ModelStatus desc,PANDDA_site_event_index"
            #      )

            # query below is without the LigandConfidence being constrained; this is because some older DBs don't have a starting digit
            # here we constrain RefinementOutcome of site
            #      cur.execute("select p.ID,p.CrystalName,p.PANDDA_site_event_index,p.CrystalName || '_event'|| p.PANDDA_site_event_index as ModelName,m.CompoundCode,m.CompoundSMILES,m.Deposition_PDB_ID,p.PANDDA_site_name,p.PANDDA_site_confidence as LigandConfidence,p.RefinementOutcome as ModelStatus,p.PANDDA_site_comment,p.PANDDA_site_x,p.PANDDA_site_y,p.PANDDA_site_z, p.PANDDA_site_spider_plot,m.DataProcessingResolutionHigh,m.DataProcessingSpaceGroup,m.DataProcessingUnitCell,m.RefinementBoundConformation,m.RefinementMTZ_latest,p.PANDDA_site_event_map from panddaTable as p, mainTable as m where p.CrystalName=m.CrystalName and p.PANDDA_site_ligand_placed='True' and (p.RefinementOutcome like '4%' or p.RefinementOutcome like '5%' or p.RefinementOutcome like '6%')  order by p.CrystalName,ModelStatus desc,PANDDA_site_event_index")

            sql = (
                "select p.ID,p.CrystalName,p.PANDDA_site_event_index,p.CrystalName || '_event'|| p.PANDDA_site_event_index "
                " as ModelName,m.CompoundCode,m.CompoundSMILES,m.Deposition_PDB_ID,p.PANDDA_site_name,"
                " p.PANDDA_site_confidence as LigandConfidence,"
                " p.RefinementOutcome as ModelStatus,"
                " p.PANDDA_site_comment,p.PANDDA_site_x,p.PANDDA_site_y,p.PANDDA_site_z, p.PANDDA_site_spider_plot,"
                " m.DataProcessingResolutionHigh,m.DataProcessingSpaceGroup,m.DataProcessingUnitCell,"
                " m.RefinementBoundConformation,m.RefinementMTZ_latest,"
                " p.PANDDA_site_event_map from panddaTable as p, "
                " mainTable as m where p.CrystalName=m.CrystalName and p.PANDDA_site_ligand_placed='True' "
                " and (p.RefinementOutcome like '4%' or p.RefinementOutcome like '5%' or p.RefinementOutcome like '6%')  "
                " and (LigandConfidence like '1%' or LigandConfidence like '2%' or LigandConfidence like '3%' or LigandConfidence like '4%')"
                " order by p.CrystalName,ModelStatus desc,PANDDA_site_event_index"
            )

            sql = (
                "select p.ID,p.CrystalName,p.PANDDA_site_event_index,p.CrystalName || '_event'|| p.PANDDA_site_event_index "
                " as ModelName,m.CompoundCode,m.CompoundSMILES,m.Deposition_PDB_ID,p.PANDDA_site_name,"
                " p.PANDDA_site_confidence as LigandConfidence,"
                " p.RefinementOutcome as ModelStatus,"
                " p.PANDDA_site_comment,p.PANDDA_site_x,p.PANDDA_site_y,p.PANDDA_site_z, p.PANDDA_site_spider_plot,"
                " m.DataProcessingResolutionHigh,m.DataProcessingSpaceGroup,m.DataProcessingUnitCell,"
                " m.RefinementBoundConformation,m.RefinementMTZ_latest,"
                " p.PANDDA_site_event_map from panddaTable as p, "
                " mainTable as m where p.CrystalName=m.CrystalName and p.PANDDA_site_ligand_placed='True' "
                " and (m.RefinementOutcome like '4%' or m.RefinementOutcome like '5%' or m.RefinementOutcome like '6%')  "
                " and (LigandConfidence like '1%' or LigandConfidence like '2%' or LigandConfidence like '3%' or LigandConfidence like '4%')"
                " order by p.CrystalName,ModelStatus desc,PANDDA_site_event_index"
            )

            cur.execute(sql)

            rows = cur.fetchall()
            if not rows:
                print(
                    "==> WARNING: none of your samples seems to be at least CompChem ready (4)"
                )
                return None
            writer = csv.DictWriter(f, fieldnames=list(rows[1].keys()))
            writer.writeheader()
            for row in rows:
                # Make compound structure
                print(row["ModelName"], row["PANDDA_site_spider_plot"])
                compound = Chem.MolFromSmiles(row["CompoundSMILES"].encode("ascii"))
                Draw.MolToFile(
                    compound,
                    panddadir + "/compoundImages/" + row["CompoundCode"] + ".png",
                    (150, 150),
                )
                # Write out table information for event
                eventID = row["ModelName"] + "_" + row["CompoundCode"]
                actID = (row["ModelName"] + row["CompoundCode"]).replace(
                    targetID + "-", ""
                )
                writeTableRow(row, htmlfile)
                writeICBPage(row, panddadir)
                try:
                    shutil.copy(
                        row["RefinementBoundConformation"],
                        panddadir + "/pdbs/" + row["ModelName"] + ".pdb",
                    )
                    shutil.copy(
                        row["RefinementMTZ_latest"],
                        panddadir + "/maps/" + row["ModelName"] + ".mtz",
                    )
                    shutil.copy(
                        row["PANDDA_site_event_map"],
                        panddadir + "/maps/" + row["ModelName"] + ".ccp4",
                    )
                    if row["PANDDA_site_spider_plot"] is not None:
                        shutil.copy(
                            row["PANDDA_site_spider_plot"],
                            panddadir + "/residueplots/" + row["ModelName"] + ".png",
                        )
                except (IOError, TypeError):
                    print(
                        "*** WARNING: cannot find PDB and/or MTZ of "
                        + row["ModelName"]
                        + " ***"
                    )
                    print("PDB bound  :", row["RefinementBoundConformation"])
                    print("MTZ        :", row["RefinementMTZ_latest"])
                    print("event map  :", row["PANDDA_site_event_map"])
                    print("spider plot:", row["PANDDA_site_spider_plot"])
                    pass
                #        shutil.copy(row['RefinementPDB_latest'],panddadir+"/pdbs/"+row['ModelName']+".pdb")
                #        if row['PANDDA_site_spider_plot'] is not None:
                #          shutil.copy(row['PANDDA_site_spider_plot'],panddadir+"/residueplots/"+row['ModelName']+".png")
                # Write row to CSV for ICM
                writer.writerow(dict(row))

    # Conclude HTML
    htmlfile.write("</tbody>\n")
    htmlfile.write("</table>\n")
    htmlfile.write("</body>\n")
    htmlfile.write("</html>\n")
    htmlfile.close()

    # Copy JS & CSS files
    if not os.path.exists(panddadir + "/js"):
        os.makedirs(panddadir + "/js")
    if not os.path.exists(panddadir + "/css"):
        os.makedirs(panddadir + "/css")
    shutil.copy(
        os.path.join(
            os.getenv("XChemExplorer_DIR"), "web/jscss/css/jquery.dataTables.min.css"
        ),
        panddadir + "/css/jquery.dataTables.min.css",
    )
    shutil.copy(
        os.path.join(
            os.getenv("XChemExplorer_DIR"), "web/jscss/js/jquery-1.12.3.min.js"
        ),
        panddadir + "/js/jquery-1.12.3.min.js",
    )
    shutil.copy(
        os.path.join(
            os.getenv("XChemExplorer_DIR"), "web/jscss/js/jquery.dataTables.min.js"
        ),
        panddadir + "/js/jquery.dataTables.min.js",
    )

    # Create zip files
    print("Creating zipfile of PDBs...")
    os.chdir(panddadir + "/pdbs")
    zf = zipfile.ZipFile("allPDBs.zip", "w")
    for pdb in glob.glob("*.pdb"):
        zf.write(pdb)
    zf.close()

    print("Creatig zipfile of event maps...")
    os.chdir("../maps")
    zf = zipfile.ZipFile("allEventMaps.zip", "w")
    for pdb in glob.glob("*.mtz"):
        zf.write(pdb)
    zf.close()

    # change folder permissions
    os.system("chmod -R 775 {0!s}".format(panddadir))

    return


if __name__ == "__main__":
    main(sys.argv[1:])
