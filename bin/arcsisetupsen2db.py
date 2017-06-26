#! /usr/bin/env python

"""
Module that contains the ARSCI command to set up database of Sentinel-2 imagery.
"""

############################################################################
#  arcsisetupsen2db.py
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
# Purpose:  A script to set up a sqlite database with Sentinel-2
#           imagery from Google with links to download the scenes  
#           from the Google bucket.
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
# Import the python glob module
import glob
# Import the python Argument parser
import argparse
# Import the python curl option 
import pycurl
# Import the tempfile python module
import tempfile
# Import python time module
import time
# Import python shutil module
import shutil
# Import python gzip module
import gzip
# Import python sqlite3 module
import sqlite3
# Import the arcsi version number
from arcsilib import ARCSI_VERSION
# Import the ARCSI utilities class
from arcsilib.arcsiutils import ARCSIUtils
# Import the ARCSI exception class
from arcsilib.arcsiexception import ARCSIException

def downloadProgress(download_t, download_d, upload_t, upload_d):
    try:
        frac = float(download_d)/float(download_t)
    except:
        frac = 0
    sys.stdout.write("\r%s %3i%%" % ("Download:", frac*100)  )

def setupSen2DB(dbFile):
    try:
        with tempfile.TemporaryDirectory() as tmpdirname:        
            sen2URL = 'https://storage.googleapis.com/gcp-public-data-sentinel-2/index.csv.gz'
    
            ggCSVSen2FileGZ = os.path.join(tmpdirname, 'index.csv.gz')
            
            fp = open(ggCSVSen2FileGZ, "wb")
            
            curl = pycurl.Curl()
            curl.setopt(pycurl.URL, sen2URL)
            curl.setopt(pycurl.FOLLOWLOCATION, True)
            curl.setopt(pycurl.NOPROGRESS, 0)
            curl.setopt(pycurl.PROGRESSFUNCTION, downloadProgress)
            curl.setopt(pycurl.FOLLOWLOCATION, 1)
            curl.setopt(pycurl.MAXREDIRS, 5)
            curl.setopt(pycurl.CONNECTTIMEOUT, 50)
            curl.setopt(pycurl.TIMEOUT, 2000)
            curl.setopt(pycurl.FTP_RESPONSE_TIMEOUT, 600)
            curl.setopt(pycurl.NOSIGNAL, 1)
            curl.setopt(pycurl.WRITEDATA, fp)
            try:
                print("Downloading image index from Google:")
                print("Start time: " + time.strftime("%c"))
                curl.perform()
                print("\nTotal-time: " + str(curl.getinfo(curl.TOTAL_TIME)))
                print("Download speed: %.2f bytes/second" % (curl.getinfo(curl.SPEED_DOWNLOAD)))
                print("Document size: %d bytes" % (curl.getinfo(curl.SIZE_DOWNLOAD)))
            except:
                raise ARCSIException("Failed to download file from Google.")
            curl.close()
            fp.close()
            sys.stdout.flush()
            
            print("Create and load data into a sqlite db:")   
            ggSen2DBConn = sqlite3.connect(dbFile)
            ggSen2DBConn.execute('''CREATE TABLE sen2 (COUNT PRIMARY KEY, GRANULE_ID text, PRODUCT_ID text, DATATAKE_IDENTIFIER text, MGRS_TILE text, SENSING_TIME text, TOTAL_SIZE text, CLOUD_COVER real, GEOMETRIC_QUALITY_FLAG int1, GENERATION_TIME text, NORTH_LAT real, SOUTH_LAT real, WEST_LON real, EAST_LON real, BASE_URL text)''')
            
            with gzip.open(ggCSVSen2FileGZ,'r') as ggCSVSen2File:
                commitCounter = 0
                keyCount = 0
                committed = False
                first = True
                for line in ggCSVSen2File:
                    if first:
                        first = False
                    else:
                        committed = False
                        line = line.decode().strip()
                        lineComps = line.split(',')
                        #print('Got:', lineComps[7])
                        if lineComps[7].strip() == 'PASSED':
                            lineComps[7] = '1'
                        elif lineComps[7].strip() == 'FAILED':
                            lineComps[7] = '0'
                        sqlcmd = "INSERT INTO sen2 (COUNT, GRANULE_ID, PRODUCT_ID, DATATAKE_IDENTIFIER, MGRS_TILE, SENSING_TIME, TOTAL_SIZE, CLOUD_COVER, GEOMETRIC_QUALITY_FLAG, GENERATION_TIME, NORTH_LAT, SOUTH_LAT, WEST_LON, EAST_LON, BASE_URL) VALUES ("+str(keyCount)+", '"+lineComps[0]+"', '"+lineComps[1]+"', '"+lineComps[2]+"', '"+lineComps[3]+"', '"+lineComps[4]+"', '"+lineComps[5]+"', "+lineComps[6]+", "+lineComps[7]+", '"+lineComps[8]+"', "+lineComps[9]+", "+lineComps[10]+", "+lineComps[11]+", "+lineComps[12]+", '"+lineComps[13]+"')"
                        #print(sqlcmd)
                        ggSen2DBConn.execute(sqlcmd)
                        keyCount = keyCount + 1
                        if commitCounter == 5000:
                            ggSen2DBConn.commit()
                            sys.stdout.write("#")
                            sys.stdout.flush()
                            commitCounter = 0
                            committed = True
                        else:
                            commitCounter = commitCounter + 1
            if not committed:
                ggSen2DBConn.commit()
            ggSen2DBConn.close()
            print("\nFinished loading data\n")
    except ARCSIException as e:
        print("Error: {}".format(e), file=sys.stderr)
    except Exception as e:
        print("Error: {}".format(e), file=sys.stderr)            


if __name__ == '__main__':
    """
    The command line user interface to ARCSI setup Sentinel-2 DB.
    """
    parser = argparse.ArgumentParser(prog='arcsisetupsen2db.py',
                                    description='''ARSCI command to set up 
                                                   database of Sentinel-2 imagery''',
                                    epilog='''A tool to set up a sqlite database
                                              with the Google Sentinel-2 imagery
                                              list.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)
    # Define the argument for specifying the input directory to be processed.
    parser.add_argument("-f", "--dbfile", type=str, required=True,
                        help='''Path to the database file.''')
    # Call the parser to parse the arguments.
    args = parser.parse_args()

    setupSen2DB(args.dbfile)


