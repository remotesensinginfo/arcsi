#! /usr/bin/env python

"""
Module that contains the ARSCI command to query database of Sentinel-2 imagery
and create a list of image download commands.
"""

############################################################################
#  arcsigensen2downlst.py
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
# Purpose:  A script to query the database created by the arcsisetupsen2db.py
#           command and create a list of URLs to download the imagery.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 26/06/2017
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

import sys
import argparse
from arcsilib import ARCSI_VERSION
from arcsilib.arcsiexception import ARCSIException
import rsgislib


def genSen2DownloadListGoogle(
    dbFile,
    tile,
    outFile,
    outpath,
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
    # Import python sqlite3 module
    import sqlite3

    try:
        ggSen2DBConn = sqlite3.connect(dbFile)
        ggSen2DBCursor = ggSen2DBConn.cursor()

        queryVar = [tile]
        query = "SELECT BASE_URL FROM SEN2 WHERE MGRS_TILE = ?"

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
        for row in ggSen2DBCursor.execute(query, queryVar):
            if lstCmds:
                cmdLst.append("gsutil " + multiStr + " cp -r " + row[0] + " " + outpath)
            else:
                cmdLst.append(row[0])

        rsgislib.tools.utils.write_list_to_file(cmdLst, outFile)

    except ARCSIException as e:
        print("Error: {}".format(e), file=sys.stderr)
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)


if __name__ == "__main__":
    """
    The command line user interface to ARCSI generate Sentinel-2 file download list.
    """
    parser = argparse.ArgumentParser(
        prog="arcsigensen2downlst.py",
        description="""ARSCI command to query 
                                                   database of Sentinel-2 imagery""",
        epilog="""A tool to query the sqlite database
                                              with the Google Sentinel-2 imagery
                                              to create a list of URLs to download.""",
    )
    # Request the version number.
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s version " + ARCSI_VERSION
    )
    # Define the argument for specifying the input directory to be processed.
    parser.add_argument(
        "-f", "--dbfile", type=str, help="""Path to the database file."""
    )
    parser.add_argument(
        "-t",
        "--tile",
        type=str,
        required=True,
        help="""Sentinel-2 tile - note remove the preceeding 'T'.""",
    )
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
        help="""Output path for the sentinel-2 SAFE files to download to on your system.""",
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

    if (args.dbfile == None) or (args.dbfile == ""):
        raise Exception(
            "A database file is required for generating download list from Google."
        )
    genSen2DownloadListGoogle(
        args.dbfile,
        args.tile,
        args.output,
        args.outpath,
        args.cloudcover,
        args.startdate,
        args.enddate,
        args.limit,
        args.multi,
        args.lstcmds,
    )
