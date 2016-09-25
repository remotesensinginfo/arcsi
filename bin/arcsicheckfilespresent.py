#! /usr/bin/env python

"""
Module that contains the ARSCI command to check that all expected files are
present following data extraction.
"""

############################################################################
#  arcsicheckfilespresent.py
#
#  Copyright 2016 ARCSI.
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
# Purpose:  A script to check that all expected files are present following
#           data extraction.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 27/02/2016
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

# Import the future functionality (for Python 2)
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
# Import the python os.path module
import os.path
# Import the python sys module
import sys
# Import the python glob module
import glob
# Import the python Argument parser
import argparse
# Import the arcsi version number
from arcsilib import ARCSI_VERSION
# Import os.walk to navigate directory structure.
import os
# Import the ARCSI sensor factory class
from arcsilib.arcsiutils import ARCSISensorFactory
# Import the list of sensors arcsi supports
from arcsilib import ARCSI_SENSORS_LIST
# Import shutil module
import shutil
# Import the list of archive file extensions arcsi supports
from arcsilib import ARCSI_ARCHIVE_EXE_LIST

class ARCSICheckFilesPresent (object):

    def getListOfFiles(self, searchDIR, headerEnding):
        outFiles = []
        for dirName, subdirList, fileList in os.walk(searchDIR):
            for fname in fileList:
                fname = str(fname)
                if fname.endswith(headerEnding):
                    outFiles.append(os.path.abspath(os.path.join(dirName, fname)))
        return outFiles

    def getListOfArchives(self, dirPath):
        archPaths = []
        for exe in ARCSI_ARCHIVE_EXE_LIST:
            fileList = glob.glob(os.path.join(dirPath, "*" + exe))
            archPaths.extend(fileList)
        return archPaths

    def findAssociatedArchive(self, inputDIR, refDIR, archList):
        refDIR = refDIR.replace(inputDIR, "")
        if (refDIR[0] == '/') or (refDIR[0] == '\\'):
            refDIR = refDIR[1:]
        outArch = None
        for arch in archList:
            if arch.count(refDIR) > 0:
                outArch = arch
                break
        return outArch

    def checkExtractedFiles(self, inputDIR, outputFile, headerEnding, sensorStr, archivesDIR):
        inputDIR = os.path.abspath(inputDIR)
        archivesDIR = os.path.abspath(archivesDIR)
        archList = self.getListOfArchives(archivesDIR)
        dirList = glob.glob(os.path.join(inputDIR, "*"))

        sensorFact = ARCSISensorFactory()
        outFileList = open(outputFile, 'w')
        foundErr = False
        for dir in dirList:
            foundErr = False
            if os.path.isdir(dir):
                print(dir)
                hdrList = self.getListOfFiles(dir, headerEnding)
                if len(hdrList) == 0:
                    print("Error: No header present")
                    foundErr = True
                elif len(hdrList) > 1:
                    print("Error: Multiple header files present - don't know what to do with that!")
                    foundErr = True
                else:
                    sensorClass = sensorFact.getSensorClassFromName(sensorStr, False, None)
                    sensorClass.extractHeaderParameters(hdrList[0], "")
                    if not sensorClass.expectedImageDataPresent():
                        print("Error: Images specified in input header file are not present")
                        foundErr = True
                if foundErr:
                    arch = self.findAssociatedArchive(inputDIR, dir, archList)
                    outFileList.write(arch + '\n')
                    shutil.rmtree(dir)

        outFileList.flush()
        outFileList.close()

if __name__ == '__main__':
    """
    The command line user interface to ARCSI.
    """
    parser = argparse.ArgumentParser(prog='arcsicheckfilespresent.py',
                                    description='''ARCSI command to check that all expected files are
                                                   present following data extraction (expecting arcsiextractdata.py
                                                   to have been used to extract data into folder structure).''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)

    parser.add_argument("-i", "--input", type=str, required=True,
                        help='''Input directory containing the folders for the extracted data''')

    parser.add_argument("-o", "--output", type=str, required=True,
                        help='''Output text file listing the archives which have problems''')

    parser.add_argument("-e", "--header", type=str, required=True,
                        help='''The extension / unquie file ending for the input header files.''')

    parser.add_argument("-s", "--sensor", required=True, choices=ARCSI_SENSORS_LIST,
                        help='''Specify the sensor being processed.''')

    parser.add_argument("-a", "--archives", type=str, required=True,
                        help='''Input directory containing the original archives''')

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    arcsiObj = ARCSICheckFilesPresent()

    arcsiObj.checkExtractedFiles(args.input, args.output, args.header, args.sensor, args.archives)

