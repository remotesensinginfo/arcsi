#! /usr/bin/env python

"""
Module that contains the ARSCI command to query database of Landsat imagery
and create a list of image download commands.
"""

############################################################################
#  arcsigenlandsatdownlst.py
#
#  Copyright 2017 ARCSI.
#
#  ARCSI: 'Atmospheric and Radiometric Correction of Satellite Imagery'
#
#  ARCSI is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ARCSI is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with ARCSI.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Purpose:  A script to query the database created by the arcsisetuplandsatdb.py
#           command and create a list of URLs to download the imagery.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 28/06/2017
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

import argparse
import sqlite3
import sys

import rsgislib

from arcsilib import ARCSI_VERSION
from arcsilib.arcsiexception import ARCSIException


def genLandsatDownloadList(
    dbFile,
    lsPath,
    lsRow,
    outFile,
    outpath,
    sensorID=None,
    spacecraftID=None,
    collection=None,
    cloudCover=None,
    startDate=None,
    endDate=None,
    limit=None,
    multiDwn=False,
    lstCmds=False,
):
    """
    Using sqlite database query and create a list of files to download
    """
    try:
        ggLandsatDBConn = sqlite3.connect(dbFile)
        ggLandsatDBCursor = ggLandsatDBConn.cursor()

        queryVar = [lsPath, lsRow]
        query = "SELECT BASE_URL FROM LANDSAT WHERE WRS_PATH = ? AND WRS_ROW = ?"

        if not sensorID is None:
            query = query + " AND SENSOR_ID = ?"
            queryVar.append(sensorID)

        if not spacecraftID is None:
            query = query + " AND SPACECRAFT_ID = ?"
            queryVar.append(spacecraftID)

        if not collection is None:
            if collection == "PRE":
                collection = "N/A"
            query = query + " AND COLLECTION_CATEGORY = ?"
            queryVar.append(collection)

        if not cloudCover is None:
            query = query + " AND CLOUD_COVER < ?"
            queryVar.append(cloudCover)

        if not startDate is None:
            query = query + " AND date(SENSING_TIME) > date(?)"
            queryVar.append(startDate)

        if not endDate is None:
            query = query + " AND date(SENSING_TIME) < date(?)"
            queryVar.append(endDate)

        if not limit is None:
            query = query + " ORDER BY CLOUD_COVER ASC LIMIT {}".format(limit)

        multiStr = ""
        if multiDwn:
            multiStr = "-m"

        cmdLst = []
        for row in ggLandsatDBCursor.execute(query, queryVar):
            if lstCmds:
                cmdLst.append(
                    "gsutil {} cp -r {} {} ".format(multiStr, row[0], outpath)
                )
            else:
                cmdLst.append(row[0])

        rsgislib.tools.utils.write_list_to_file(cmdLst, outFile)

    except ARCSIException as e:
        print("Error: {}".format(e), file=sys.stderr)
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)


if __name__ == "__main__":
    """
    The command line user interface to ARCSI generate Landsat file download list.
    """
    parser = argparse.ArgumentParser(
        prog="arcsigenlandsatdownlst.py",
        description="""ARSCI command to query 
                                                   database of Landsat imagery""",
        epilog="""A tool to query the sqlite database
                                              with the Google Landsat imagery
                                              to create a list of URLs to download.""",
    )
    # Request the version number.
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s version " + ARCSI_VERSION
    )
    # Define the argument for specifying the input directory to be processed.
    parser.add_argument(
        "-f", "--dbfile", type=str, required=True, help="""Path to the database file."""
    )
    parser.add_argument(
        "-p", "--path", type=str, required=True, help="""Landsat path."""
    )
    parser.add_argument("-r", "--row", type=str, required=True, help="""Landsat row.""")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="""Output file with a list of files to download.""",
    )
    parser.add_argument(
        "--outpath",
        type=str,
        required=True,
        help="""Output path for the landsat files to download to on your system.""",
    )
    parser.add_argument(
        "--sensor",
        type=str,
        choices=["OLI_TIRS", "ETM", "TM", "MSS", "MSS"],
        help="""Specify the landsat sensor you are interested""",
    )
    parser.add_argument(
        "--spacecraft",
        type=str,
        choices=[
            "LANDSAT_8",
            "LANDSAT_7",
            "LANDSAT_5",
            "LANDSAT_4",
            "LANDSAT_3",
            "LANDSAT_2",
            "LANDSAT_1",
        ],
        help="""Specify the landsat spacecraft you are interested""",
    )
    parser.add_argument(
        "--collection",
        type=str,
        choices=["T1", "T2", "RT", "PRE"],
        help="""Specify the landsat collection you are interested. For more information see https://landsat.usgs.gov/landsat-collections""",
    )
    parser.add_argument(
        "--cloudcover",
        type=float,
        help="""Specify an upper limit for acceptable cloud cover.""",
    )
    parser.add_argument(
        "--startdate", type=str, help="""Specify a start date (YYYY-MM-DD)."""
    )
    parser.add_argument(
        "--enddate", type=str, help="""Specify a end date (YYYY-MM-DD)."""
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="""Specify a limit for the number of scenes returned - scenes are sorted by cloud cover""",
    )
    parser.add_argument(
        "--multi",
        action="store_true",
        default=False,
        help="""Adds -m option to the gsutil download command.""",
    )
    parser.add_argument(
        "--lstcmds",
        action="store_true",
        default=False,
        help="""List download commands rather than just list of URLs""",
    )

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    genLandsatDownloadList(
        args.dbfile,
        args.path,
        args.row,
        args.output,
        args.outpath,
        args.sensor,
        args.spacecraft,
        args.collection,
        args.cloudcover,
        args.startdate,
        args.enddate,
        args.limit,
        args.multi,
        args.lstcmds,
    )
