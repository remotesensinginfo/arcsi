#! /usr/bin/env python

"""
Module that contains the ARSCI command to check that all expected files are
present following data extraction.
"""

############################################################################
#  arcsiarchivesnotextracted.py
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
# Date: 28/02/2016
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

# Import the python os.path module
import os.path
# Import the python glob module
import glob
# Import the python Argument parser
import argparse
# Import the arcsi version number
from arcsilib import ARCSI_VERSION
# Import the list of archive file extensions arcsi supports
from arcsilib import ARCSI_ARCHIVE_EXE_LIST

class ARCSIFindArchivesNotExtract (object):

    def getListOfArchives(self, dirPath):
        archPaths = []
        for exe in ARCSI_ARCHIVE_EXE_LIST:
            fileList = glob.glob(os.path.join(dirPath, "*" + exe))
            archPaths.extend(fileList)
        return archPaths

    def checkExtractedFiles(self, inputDIR, outputFile, archivesDIR):
        inputDIR = os.path.abspath(inputDIR)
        archivesDIR = os.path.abspath(archivesDIR)
        archList = self.getListOfArchives(archivesDIR)
        outFileList = open(outputFile, 'w')
        for arch in archList:
            print("Checking: " + arch)
            dirName = os.path.basename(arch).split(".")[0]
            dirPath = os.path.join(inputDIR, dirName)
            print("\t Testing for: " + dirPath)
            if not os.path.exists(dirPath):
                outFileList.write(arch + "\n")
            elif not os.path.isdir(dirPath):
                outFileList.write(arch + "\n")

        outFileList.flush()
        outFileList.close()

if __name__ == '__main__':
    """
    The command line user interface to ARCSI.
    """
    parser = argparse.ArgumentParser(prog='arcsiarchivesnotextracted.py',
                                    description='''ARCSI command to check that all expected files are
                                                   present following data extraction (expecting arcsiextractdata.py
                                                   to have been used to extract data into folder structure).''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)

    parser.add_argument("-i", "--input", type=str, required=True,
                        help='''Input directory containing the folders for the extracted data''')

    parser.add_argument("-o", "--output", type=str, required=True,
                        help='''Output text file listing the archives which have problems''')

    parser.add_argument("-a", "--archives", type=str, required=True,
                        help='''Input directory containing the original archives''')

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    arcsiObj = ARCSIFindArchivesNotExtract()

    arcsiObj.checkExtractedFiles(args.input, args.output, args.archives)

