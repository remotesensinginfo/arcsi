#! /usr/bin/env python

"""
Module that contains the ARSCI command to build lists of image header files
for arcsi.py --multi.
"""

############################################################################
#  arcsibuildmultifilelists.py
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
# Purpose:  A script to build lists of input header files which are inputted
#           into arcsi.py --multi
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 08/09/2017
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

import os
import fnmatch
import argparse
from arcsilib import ARCSI_VERSION
import os
from arcsilib.arcsiutils import ARCSISensorFactory
from arcsilib import ARCSI_SENSORS_LIST


class ARCSIBuildMultiFileList(object):
    def getListOfFiles(self, searchDIR, searchStr, depth):
        inDIRCount = searchDIR.count(os.path.sep)
        outFiles = []
        for dirName, subdirList, fileList in os.walk(searchDIR):
            for fname in fileList:
                fname = str(fname)
                if fnmatch.fnmatch(fname, searchStr):
                    outFiles.append(os.path.abspath(os.path.join(dirName, fname)))
            dirLevel = dirName.count(os.path.sep) - inDIRCount
            if dirLevel >= depth:
                subdirList[:] = []
        return outFiles

    def buildCmds(self, inputPath, outputFile, headerSearchStr, searchDepth, sensor):
        """ """
        inputPath = os.path.abspath(inputPath)
        outputFile = os.path.abspath(outputFile)

        headersFilesList = []
        if headerSearchStr.count("*") == 0:
            raise Exception(
                "The search string you have provided does not have any '*' - which is needed for searching."
            )
        headersFilesList = self.getListOfFiles(inputPath, headerSearchStr, searchDepth)

        sensorFact = ARCSISensorFactory()
        sensorClass = sensorFact.getSensorClassFromName(sensor, False, None)

        hdfFileLUT = dict()

        for hdrFile in headersFilesList:
            print(hdrFile)
            sensorClass.extractHeaderParameters(hdrFile, None)
            if sensor == "sen2":
                scnStr = (
                    sensorClass.spacecraftName.replace("-", "")
                    + "_"
                    + sensorClass.acquisitionTime.strftime("%Y%m%d")
                )
            else:
                scnStr = sensor + "_" + sensorClass.acquisitionTime.strftime("%Y%m%d")
            if scnStr not in hdfFileLUT:
                hdfFileLUT[scnStr] = []
            hdfFileLUT[scnStr].append(hdrFile)

        for scnStr in hdfFileLUT:
            rsgislib.tools.utils.write_list_to_file(
                hdfFileLUT[scnStr], outputFile + scnStr + ".txt"
            )
        print("Completed.")


if __name__ == "__main__":
    """
    The command line user interface to ARCSI Data Extraction Tool.
    """
    parser = argparse.ArgumentParser(
        prog="arcsibuildmultifilelists.py",
        description="""ARCSI command to build arcsi.py commands
                                                for a set of input images using the same options.""",
        epilog="""A tools to build arcsi.py commands
                                           for a set of input images using
                                           the same options""",
    )
    # Request the version number.
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s version " + ARCSI_VERSION
    )

    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="""Input directory containing the data to be processed or text file listing paths to header files (don't set --header if text file)""",
    )

    parser.add_argument(
        "-e",
        "--header",
        type=str,
        required=True,
        help="""A \'UNIX\' search string for identifying the image headers. Note, multiple \'*\' can be used for the search string. If no \'*\' is provided then """,
    )

    parser.add_argument(
        "-d",
        "--depth",
        type=int,
        required=True,
        help="""The depth within the directory tree from the input path which should be searched for image header files.""",
    )

    parser.add_argument(
        "-s",
        "--sensor",
        required=True,
        choices=ARCSI_SENSORS_LIST,
        help="""Specify the sensor being processed.""",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="""Output base file path/name.""",
    )

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    arcsiObj = ARCSIBuildMultiFileList()

    arcsiObj.buildCmds(args.input, args.output, args.header, args.depth, args.sensor)
