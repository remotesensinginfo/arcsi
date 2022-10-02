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
# Purpose:  A script which uses the LUT to remove duplicates from the system
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 26/02/2016
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

import os
import sys
import argparse
from arcsilib import ARCSI_VERSION
from arcsilib.arcsiutils import ARCSISensorFactory
from arcsilib.arcsiexception import ARCSIException
import json
import shutil
import random


class ARCSIRemoveDuplicates(object):
    def selectFiles2Keep(
        self,
        dupHdrs,
        headersDIR,
        dupArchs,
        useArchs,
        selectVersion,
        archivesDIR,
        cpArchives2DIR,
    ):
        """ """
        if useArchs and (len(dupHdrs) != len(dupArchs)):
            raise Exception(
                "If using archives there should be the same number of headers and archives."
            )

        selKey = ""
        if selectVersion == "RANDOM":
            selKey = random.choice(list(dupHdrs.keys()))
        elif selectVersion == "LANDSAT":
            sensorFact = ARCSISensorFactory()
            first = True
            selKey = ""
            selKeyGenTime = None
            for baseName in dupHdrs:
                print(baseName)
                filePrefix3 = baseName[:3]
                filePrefix4 = baseName[:4]
                sensor = ""
                if filePrefix3 == "LM1" or filePrefix4 == "LM01":
                    sensor = "ls1"
                elif filePrefix3 == "LM2" or filePrefix4 == "LM02":
                    sensor = "ls2"
                elif filePrefix3 == "LM3" or filePrefix4 == "LM03":
                    sensor = "ls3"
                elif filePrefix3 == "LM4" or filePrefix4 == "LM04":
                    sensor = "ls4mss"
                elif filePrefix3 == "LM5" or filePrefix4 == "LM05":
                    sensor = "ls5mss"
                elif (
                    filePrefix3 == "LT4"
                    or filePrefix4 == "LS04"
                    or filePrefix4 == "LE04"
                    or filePrefix4 == "LT04"
                ):
                    sensor = "ls4tm"
                elif (
                    filePrefix3 == "LT5"
                    or filePrefix4 == "LS05"
                    or filePrefix4 == "LE05"
                    or filePrefix4 == "LT05"
                ):
                    sensor = "ls5tm"
                elif (
                    filePrefix3 == "LE7"
                    or filePrefix4 == "LS07"
                    or filePrefix4 == "LE07"
                    or filePrefix4 == "LT07"
                ):
                    sensor = "ls7"
                elif (
                    filePrefix3 == "LC8"
                    or filePrefix4 == "LS08"
                    or filePrefix4 == "LC08"
                ):
                    sensor = "ls8"
                else:
                    raise ARCSIException(
                        'Sensor was not recognised for file: "' + baseName + '"'
                    )

                hdrFullPath = os.path.join(headersDIR, dupHdrs[baseName])
                sensorClass = sensorFact.getSensorClassFromName(sensor, False, None)
                sensorClass.extractHeaderParameters(hdrFullPath, "")

                cKeyTime = sensorClass.fileDateObj
                if first:
                    selKeyGenTime = cKeyTime
                    selKey = baseName
                    first = False
                elif cKeyTime > selKeyGenTime:
                    selKeyGenTime = cKeyTime
                    selKey = baseName
        else:
            raise Exception("Don't know the selection method.")

        for baseName in dupHdrs:
            if not selKey == baseName:
                fullPath2Del = os.path.join(headersDIR, baseName)
                if os.path.exists(fullPath2Del):
                    print("Deleting: " + fullPath2Del)
                    shutil.rmtree(fullPath2Del)
                if useArchs:
                    archFile = os.path.join(archivesDIR, dupArchs[baseName])
                    archFileMV = os.path.join(cpArchives2DIR, dupArchs[baseName])
                    if os.path.isfile(archFile):
                        print("Moving: " + archFile)
                        shutil.move(archFile, archFileMV)

    def sortDuplicateFiles(
        self, lutFile, headersDIR, archivesDIR, cpArchives2DIR, selectVersion
    ):
        """ """
        headersDIR = os.path.abspath(headersDIR)
        useArchs = False
        if (archivesDIR is not None) and (cpArchives2DIR is not None):
            archivesDIR = os.path.abspath(archivesDIR)
            cpArchives2DIR = os.path.abspath(cpArchives2DIR)
            useArchs = True
        elif (archivesDIR is not None) or (cpArchives2DIR is not None):
            raise Exception(
                "If archives directory is specified then the output directory to which the duplicate archives will be moved must be specified."
            )

        with open(lutFile, "r") as f:
            jsonStrData = f.read()
        fileLUT = json.loads(jsonStrData)

        for arcsiFileName in fileLUT:
            print(arcsiFileName)
            dupArchs = dict()
            dupHdrs = dict()
            dirName = os.path.dirname(fileLUT[arcsiFileName]["Header"])
            dupHdrs[dirName] = fileLUT[arcsiFileName]["Header"]
            if useArchs:
                dupArchs[dirName] = fileLUT[arcsiFileName]["Archive"]
            if "Duplicates" in fileLUT[arcsiFileName]:
                for dup in fileLUT[arcsiFileName]["Duplicates"]:
                    if not headersDIR is None:
                        dirName = os.path.dirname(dup["Header"])
                        dupHdrs[dirName] = dup["Header"]
                        if useArchs:
                            dupArchs[dirName] = dup["Archive"]
            if len(dupHdrs) > 1:
                self.selectFiles2Keep(
                    dupHdrs,
                    headersDIR,
                    dupArchs,
                    useArchs,
                    selectVersion,
                    archivesDIR,
                    cpArchives2DIR,
                )
            print("")


if __name__ == "__main__":
    """
    The command line user interface to ARCSI.
    """
    parser = argparse.ArgumentParser(
        prog="arcsiremoveduplicates.py",
        description="""ARCSI command to remove duplicate files.""",
        epilog="""A tool to remove duplicate files using the LUT
                                    file generated from the arcsibuildfilenameslu.py command.""",
    )
    # Request the version number.
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s version " + ARCSI_VERSION
    )

    parser.add_argument(
        "-l",
        "--lut",
        type=str,
        required=True,
        help="""Look up table (as generated arcsibuildfilenameslu.py) of the arcsi file names longside input headers""",
    )

    parser.add_argument(
        "-w",
        "--workingdir",
        type=str,
        help="""Working directory from which the header file paths are reference from in the LUT""",
    )

    parser.add_argument(
        "-a",
        "--archivedir",
        type=str,
        help="""Directory containing the archvies references in the LUT.""",
    )

    parser.add_argument(
        "-d",
        "--dirout",
        type=str,
        help="""A directory where the archives will be moved to.""",
    )

    parser.add_argument(
        "-s",
        "--select",
        type=str,
        choices=["RANDOM", "LANDSAT"],
        default="RANDOM",
        help="""Specify whether the file kept is selected at random or whether for landsat images the generation time is used to select the file generated more recently.""",
    )

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    if (args.workingdir is None) and (args.archivedir is None):
        print("Either or both the working and archives directories need to be defined")
        sys.exit()

    if not args.archivedir is None:
        if args.dirout is None:
            print(
                "If archives directory is specified then the output directory to which the duplicate archives will be moved must be specified."
            )
            sys.exit()

    arcsiObj = ARCSIRemoveDuplicates()
    arcsiObj.sortDuplicateFiles(
        args.lut, args.workingdir, args.archivedir, args.dirout, args.select
    )
