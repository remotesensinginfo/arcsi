#! /usr/bin/env python

"""
Module that contains the ARSCI command to check the files within a text file
which is ready for processing in arcsi with the --multi option have only a 
single input image per granule and select higher processing version.
"""

############################################################################
#  arcsichecksen2ver.py
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
# Purpose:  Module that contains the ARSCI command to check the 
#           files within a text file which is ready for processing 
#           in arcsi with the --multi option have only a single input 
#           image per granule and select higher processing version.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 12/11/2017
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
# Import the python Argument parser
import argparse
# Import the arcsi version number
from arcsilib import ARCSI_VERSION
# Import the ARCSI utilities class
from arcsilib.arcsiutils import ARCSIUtils
# Import the ARCSI sensor factory class
from arcsilib.arcsiutils import ARCSISensorFactory
# Import the ARCSI exception class
from arcsilib.arcsiexception import ARCSIException

class Sen2InfoObj(object):

    def __init__(self, hdrFile=None, tileGranuleID=None, processVer=None, generationTime=None, orbitNumber=None, orbitDirection=None):
        self.hdrFile = hdrFile
        self.tileGranuleID = tileGranuleID
        self.processVer = processVer
        self.generationTime = generationTime
        self.orbitNumber = orbitNumber
        self.orbitDirection = orbitDirection


class ARCSICheckSen2FileVersions (object):

    def checkFileVersions(self, inputFile, outputFile):
        """
        """
        arcsiUtils = ARCSIUtils()
        hdrFilesLst = arcsiUtils.readTextFile2List(inputFile)

        granuleLUT = dict()
        for hdrFile in hdrFilesLst:
            print(hdrFile)
            sensorFact = ARCSISensorFactory()
            sensorClass = sensorFact.getSensorClassFromName('sen2', False, None)
            sensorClass.extractHeaderParameters(hdrFile, None)
            tileGranuleID = sensorClass.uniqueTileID
            processVer = sensorClass.processingBaseline
            genTime = sensorClass.generationTime
            orbNum = sensorClass.orbitNumber
            orbDirect = sensorClass.orbitDirection
            scnTileID = tileGranuleID+'_'+orbNum+'_'+orbDirect
            print(scnTileID)
            if scnTileID not in granuleLUT:
                granuleLUT[scnTileID] = list()
            granuleLUT[scnTileID].append(Sen2InfoObj(hdrFile, tileGranuleID, processVer, genTime, orbNum, orbDirect))
            sensorClass = None

        outHdrs = list()
        for scnTileID in granuleLUT:
            if len(granuleLUT[scnTileID]) > 1:
                maxGenTime = 0
                maxHdrHeader = ""
                first = True
                for sen2Info in granuleLUT[scnTileID]:
                    cGenTime = sen2Info.generationTime
                    if first:
                        first = False
                        maxHdrHeader = sen2Info.hdrFile
                        maxGenTime = cGenTime
                    elif cGenTime > maxGenTime:
                        maxHdrHeader = sen2Info.hdrFile
                        maxGenTime = cGenTime
                outHdrs.append(maxHdrHeader)
            else:
                outHdrs.append(granuleLUT[scnTileID][0].hdrFile)

        arcsiUtils.writeList2File(outHdrs, outputFile)

if __name__ == '__main__':
    """
    The command line user interface to ARCSI Data Extraction Tool.
    """
    parser = argparse.ArgumentParser(prog='arcsichecksen2ver.py',
                                    description='''ARCSI command to check multiple versions of the same Sen2 granule are not being submitted for processing.''',
                                    epilog='''A tools to check multiple versions of the same Sen2 granule are not being submitted for processing.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)

    parser.add_argument("-i", "--input", type=str, required=True,
                        help='''Input text file with a list of header files.''')

    parser.add_argument("-o", "--output", type=str, required=True,
                        help='''Output text file with a list of header files.''')

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    arcsiObj = ARCSICheckSen2FileVersions()

    arcsiObj.checkFileVersions(args.input, args.output)
