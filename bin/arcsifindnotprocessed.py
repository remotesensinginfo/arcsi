#! /usr/bin/env python

"""
Module that contains the ARSCI command to using the files lookup table
to find scene which have yet to be processed to completion.
"""

############################################################################
#  arcsifindnotprocessed.py
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
# Purpose:  A script to build a LUT of the output files names alongside the input
#           file names.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 25/02/2016
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
# Import JSON module
import json

class ARCSIBuildFileNameLUT (object):

    def getListOfFiles(self, searchDIR, fileEnding):
        outFiles = []
        for dirName, subdirList, fileList in os.walk(searchDIR):
            for fname in fileList:
                fname = str(fname)
                if fname.endswith(fileEnding):
                    outFiles.append(os.path.abspath(os.path.join(dirName, fname)))
        return outFiles

    def findUnprocessedFiles(self, inputDIR, fileEnding, outputFile, lutFile, headersPath):
        inputDIR = os.path.abspath(inputDIR)
        fileList = self.getListOfFiles(inputDIR, fileEnding)

        with open(lutFile, 'r') as f:
            jsonStrData = f.read()
        fileLUT = json.loads(jsonStrData)

        if not headersPath is None:
            headersPath = os.path.abspath(headersPath)

        outFile = open(outputFile, 'w+')
        for fileBase in fileLUT:
            print(fileBase)
            found = False
            for file in fileList:
                if file.count(fileBase) > 0:
                    print("\tFound: " + file)
                    found = True
                    break
            if not found:
                headerFile = fileLUT[fileBase]['Header']
                if not headersPath is None:
                    headerFile = os.path.join(headersPath, headerFile)
                outFile.write(headerFile + "\n")
        outFile.flush()
        outFile.close()

if __name__ == '__main__':
    """
    The command line user interface to ARCSI.
    """
    parser = argparse.ArgumentParser(prog='arcsifindnotprocessed.py',
                                    description='''ARCSI command to build LUT of file names vs input files.''',
                                    epilog='''A tools to build LUT of file names vs input headers files and
                                    optionally the archives. Output is a JSON file.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)

    parser.add_argument("-i", "--input", type=str, required=True,
                        help='''Input directory containing the data to be checked''')

    parser.add_argument("-o", "--output", type=str, required=True,
                        help='''Output text file list of header files which have not been processed.''')

    parser.add_argument("-e", "--ending", type=str, required=True,
                        help='''The extension / unquie file ending for the input files.''')

    parser.add_argument("-l", "--lut", type=str,  required=True,
                        help='''Look up table (as generated arcsibuildfilenameslu.py) of the arcsi file names longside input headers''')

    parser.add_argument("-p", "--headerspath", type=str,
                        help='''Path to location for the header files which will be joined into the output.''')

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    arcsiObj = ARCSIBuildFileNameLUT()

    arcsiObj.findUnprocessedFiles(args.input, args.ending, args.output, args.lut, args.headerspath)

