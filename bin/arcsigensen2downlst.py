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

# Import updated print function into python 2.7
from __future__ import print_function
# Import updated division operator into python 2.7
from __future__ import division
# Import the python os.path module
import os.path
# Import the python sys module
import sys
# Import the python Argument parser
import argparse
# Import python time module
import time
# Import the arcsi version number
from arcsilib import ARCSI_VERSION
# Import the ARCSI exception class
from arcsilib.arcsiexception import ARCSIException
# Import rsgislib module
import rsgislib
# Import the json module
import json

def genSen2DownloadListGoogle(dbFile, tile, outFile, outpath, cloudCover=None, startDate=None, endDate=None, multiDwn=False, lstCmds=False):
    """
    Using sqlite database query and create a list of files to download
    """
    # Import python sqlite3 module
    import sqlite3

    try:
        ggSen2DBConn = sqlite3.connect(dbFile)
        ggSen2DBCursor = ggSen2DBConn.cursor()
        
        queryVar = [tile]
        query = 'SELECT BASE_URL FROM SEN2 WHERE MGRS_TILE = ?'
        
        if not cloudCover is None:
            query = query + ' AND CLOUD_COVER < ?'
            queryVar.append(cloudCover)
            
        if not startDate is None:
            query = query + ' AND date(SENSING_TIME) > date(?)'
            queryVar.append(startDate)
            
        if not endDate is None:
            query = query + ' AND date(SENSING_TIME) < date(?)'
            queryVar.append(endDate)

        multiStr = ''
        if multiDwn:
            multiStr = '-m'
        
        cmdLst = []
        for row in ggSen2DBCursor.execute(query,  queryVar):
            if lstCmds:
                cmdLst.append("gsutil "+multiStr+" cp -r " + row[0] + " " + outpath)
            else:
                cmdLst.append(row[0])
        
        rsgisUtils = rsgislib.RSGISPyUtils()
        rsgisUtils.writeList2File(cmdLst, outFile)
        
    except ARCSIException as e:
        print("Error: {}".format(e), file=sys.stderr)
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)

def genSen2DownloadListAWS(tile, outFile, outpath, limit=100, cloudCover=None, startDate=None, endDate=None, multiDwn=False, lstCmds=False):
    """
    Using query online AWS database and create a list of files to download
    """
    # Import python requests module
    import requests

    try:
        devSeedURL = "https://api.developmentseed.org/satellites/?limit="+str(limit)+"&search=scene_id:S2*"+tile+"*"
        if startDate is not None:
            devSeedURL = devSeedURL + "&date_from="+startDate
        if startDate is not None:
            devSeedURL = devSeedURL + "&date_to="+endDate
        if cloudCover is not None:
            devSeedURL = devSeedURL + "&cloud_to="+cloudCover
        
        multiStr = ''
        if multiDwn:
            multiStr = ' --threaded '

        r = requests.get(devSeedURL)
        if r.status_code == 200:
            cmdLst = []
            jsonObj = json.loads(r.text)
            nScenesFound = int(jsonObj['meta']['found'])
            if nScenesFound > limit:
                print("WARNING: The scene limit ("+str(limit)+") is smaller than the number of scenes available ("+str(nScenesFound)+")", file=sys.stderr)
            if nScenesFound > 0:
                resultsLst = jsonObj['results']
                for rsult in resultsLst:
                    if lstCmds:
                        cmdLst.append("sentinelhub.aws --product " + rsult["product_id"] + " -f " + outpath + multiStr)
                    else:
                        cmdLst.append(rsult["product_id"])

                rsgisUtils = rsgislib.RSGISPyUtils()
                rsgisUtils.writeList2File(cmdLst, outFile)

        else:
            raise ARCSIException("Did not get response back from the server; try again later. If persists report as bug. Return code "+str(r.status_code))
    except ARCSIException as e:
        print("Error: {}".format(e), file=sys.stderr)
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)



if __name__ == '__main__':
    """
    The command line user interface to ARCSI generate Sentinel-2 file download list.
    """
    parser = argparse.ArgumentParser(prog='arcsigensen2downlst.py',
                                    description='''ARSCI command to query 
                                                   database of Sentinel-2 imagery''',
                                    epilog='''A tool to query the sqlite database
                                              with the Google Sentinel-2 imagery
                                              to create a list of URLs to download.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)
    # Define the argument for specifying the input directory to be processed.
    parser.add_argument("-s", "--source", choices=['AWS','GOOG'], default='AWS', help='''Specify the sensor being processed.''')
    parser.add_argument("-f", "--dbfile", type=str, help='''Path to the database file.''')
    parser.add_argument("-t", "--tile", type=str, required=True, help='''Sentinel-2 tile - note remove the preceeding 'T'.''')
    parser.add_argument("-o", "--output", type=str, required=True, help='''Output file with a list of files to download.''')
    parser.add_argument("--outpath", type=str, help='''Output path for the sentinel-2 SAFE files to download to on your system.''')
    parser.add_argument("--cloudcover", type=float, help='''Specify an upper limit for acceptable cloud cover.''')
    parser.add_argument("--startdate", type=str, help='''Specify a start date (YYYY-MM-DD).''')
    parser.add_argument("--enddate", type=str, help='''Specify a end date (YYYY-MM-DD).''')
    parser.add_argument("--limit", type=int, default=100, help='''Specify the maximum number of files which are retrieved (AWS only).''')
    parser.add_argument("--multi", action='store_true', default=False, help='''Adds -m option to the gsutil download command.''')
    parser.add_argument("--lstcmds", action='store_true', default=False, help='''List download commands rather than just list of URLs''')

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    if args.source == 'GOOG':
        if (args.dbfile == None) or (args.dbfile == ""):
            raise Exception("A database file is required for generating download list from Google.")
        genSen2DownloadListGoogle(args.dbfile, args.tile, args.output, args.outpath, args.cloudcover, args.startdate, args.enddate, args.multi, args.lstcmds)
    elif args.source == 'AWS':
        genSen2DownloadListAWS(args.tile, args.output, args.outpath, args.limit, args.cloudcover, args.startdate, args.enddate, args.multi, args.lstcmds)
    else:
        raise Exception("You must specify whether to search with Google or Amazon.")

