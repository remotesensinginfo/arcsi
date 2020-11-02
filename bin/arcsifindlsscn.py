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
# Import python sqlite3 module
import sqlite3
# Import the arcsi version number
from arcsilib import ARCSI_VERSION
# Import the ARCSI exception class
from arcsilib.arcsiexception import ARCSIException
# Import rsgislib module
import rsgislib


def findLandsatSceneURL(dbFile, lsPath, lsRow, scnDate, sensorID=None, spacecraftID=None, collection=None, outCmd=False, multiDwn=False, outpath=None):
    """
    Using sqlite database query and create a list of files to download
    """
    try:
        ggLandsatDBConn = sqlite3.connect(dbFile)
        ggLandsatDBCursor = ggLandsatDBConn.cursor()
        
        queryVar = [lsPath, lsRow, scnDate]
        query = 'SELECT BASE_URL FROM LANDSAT WHERE WRS_PATH = ? AND WRS_ROW = ? AND date(SENSING_TIME) == date(?)'
        
        if not sensorID is None:
            query = query + ' AND SENSOR_ID = ?'
            queryVar.append(sensorID)

        if not spacecraftID is None:
            query = query + ' AND SPACECRAFT_ID = ?'
            queryVar.append(spacecraftID)

        if not collection is None:
            if collection == 'PRE':
                collection = 'N/A'
            query = query + ' AND COLLECTION_CATEGORY = ?'
            queryVar.append(collection)
        
        

        scn_url = ""
        found_scn = False
        found_mscns = False
        for row in ggLandsatDBCursor.execute(query,  queryVar):
            if not found_scn:
                scn_url = row[0]
                found_scn = True
            else:
                if not found_mscns:
                    print(scn_url)
                print(row[0])
                found_mscns = True
        
        if not found_scn:
            raise Exception("Was not able to find the scene...")
        
        if found_mscns:
            raise Exception("Found multiple scenes... Do not know what to do.")
        
        if outCmd:
            multiStr = ''
            if multiDwn:
                multiStr = '-m'
                
            out_path_str = '.'
            if outpath is not None:
                out_path_str = outpath
            
            cmd = "gsutil {} cp -r {} {}".format(multiStr, scn_url, out_path_str)
            print(cmd)
        else:
            print(scn_url)
        
        
    except ARCSIException as e:
        print("Error: {}".format(e), file=sys.stderr)
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)

if __name__ == '__main__':
    """
    The command line user interface to ARCSI generate Landsat file download list.
    """
    parser = argparse.ArgumentParser(prog='arcsigenlandsatdownlst.py',
                                    description='''ARSCI command to query 
                                                   database of Landsat imagery''',
                                    epilog='''A tool to query the sqlite database
                                              with the Google Landsat imagery
                                              to create a list of URLs to download.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)
    # Define the argument for specifying the input directory to be processed.
    parser.add_argument("-f", "--dbfile", type=str, required=True, help='''Path to the database file.''')
    parser.add_argument("-p", "--path", type=str, required=True, help='''Landsat path.''')
    parser.add_argument("-r", "--row", type=str, required=True, help='''Landsat row.''')
    parser.add_argument("--scndate", type=str, required=True, help='''Specify a scene date (YYYY-MM-DD).''')
    parser.add_argument("--sensor", type=str, choices=['OLI_TIRS', 'ETM', 'TM', 'MSS', 'MSS'], help='''Specify the landsat sensor you are interested''')
    parser.add_argument("--spacecraft", type=str, choices=['LANDSAT_8', 'LANDSAT_7', 'LANDSAT_5', 'LANDSAT_4', 'LANDSAT_3', 'LANDSAT_2', 'LANDSAT_1'], help='''Specify the landsat spacecraft you are interested''')
    parser.add_argument("--collection", type=str,  default='T1', choices=['T1', 'T2', 'RT', 'PRE'], help='''Specify the landsat collection you are interested. For more information see https://landsat.usgs.gov/landsat-collections''')
    parser.add_argument("--multi", action='store_true', default=False, help='''Adds -m option to the gsutil download command.''')
    parser.add_argument("--cmd", action='store_true', default=False, help='''List download command rather than just URL''')
    parser.add_argument("--outpath", type=str, help='''Output path for the landsat files to download to on your system. If not provided then will be defined as the current directory (i.e., '.')''')


    # Call the parser to parse the arguments.
    args = parser.parse_args()

    findLandsatSceneURL(args.dbfile, args.path, args.row, args.scndate, args.sensor, args.spacecraft, args.collection, args.cmd, args.multi, args.outpath)

