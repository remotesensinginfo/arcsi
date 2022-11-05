#! /usr/bin/env python

"""
Module that contains the ARSCI command to build the individual
file archive extraction commands using arcsiextractdata.py for
deployment on a HPC system.
"""

############################################################################
#  arcsibuildextractfilecmds.py
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
# Purpose:  A script to build the individual arcsiextractdata.py commands
#           from a directory of files for deployment on a multi-processor
#           system.
#
# Author: Pete Bunting
# Email:  pfb@aber.ac.uk
# Date:   30/01/2016
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
from arcsilib import ARCSI_ARCHIVE_EXE_LIST
import rsgislib.tools.utils

class ARCSIBuildExtractDataCmds(object):
    def runGenCommandsDIR(self, inputDIR, outputFile, outDIR, nofolders):
        outDIR = os.path.abspath(outDIR)
        outFile = open(outputFile, "w")
        for dirName, subdirList, fileList in os.walk(inputDIR):
            for fname in fileList:
                fileExt = os.path.basename(fname).split(".", 1)[-1].lower()
                filePath = os.path.abspath(os.path.join(dirName, fname))
                relFileExe = False
                for ext in ARCSI_ARCHIVE_EXE_LIST:
                    if ext[1:] == fileExt:
                        relFileExe = True
                if relFileExe:
                    print("\t%s" % filePath)
                    command = (
                        'arcsiextractdata.py -f "' + filePath + '"  -o "' + outDIR + '"'
                    )
                    if nofolders:
                        command = command + " --nofolders"
                    command = command + "\n"
                    outFile.write(command)
        outFile.close()

    def runGenCommandsFiles(self, inputFileList, outputFile, outDIR, nofolders):
        outDIR = os.path.abspath(outDIR)

        fileList = rsgislib.tools.utils.read_text_file_no_new_lines2List(inputFileList)
        outFile = open(outputFile, "w")
        for filePath in fileList:
            fileExt = os.path.basename(filePath).split(".", 1)[-1].lower()
            relFileExe = False
            for ext in ARCSI_ARCHIVE_EXE_LIST:
                if ext[1:] == fileExt:
                    relFileExe = True
            if relFileExe:
                print("\t%s" % filePath)
                command = (
                    'arcsiextractdata.py -f "' + filePath + '"  -o "' + outDIR + '"'
                )
                if nofolders:
                    command = command + " --nofolders"
                command = command + "\n"
                outFile.write(command)

        outFile.close()


if __name__ == "__main__":
    """
    The command line user interface to ARCSI build arcsibuildextractfilecmds.py commands.
    """
    parser = argparse.ArgumentParser(
        prog="arcsibuildextractfilecmds.py",
        description="""ARCSI command to build commands for
                                                arcsiextractdata.py.""",
        epilog="""ARCSI command to build commands for
                                              arcsiextractdata.py for extracting
                                              single files on a HPC type system.""",
    )
    # Request the version number.
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s version " + ARCSI_VERSION
    )
    # Define the argument for specifying the input directory to be processed.
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="""Input directory contains archives (tar, tar.gz, bz2, and/or zip).""",
    )
    parser.add_argument(
        "-l",
        "--listfile",
        type=str,
        help="""Input text file containing a list of archive (tar, tar.gz, bz2, and/or zip) files.""",
    )
    # Define the argument for specifying the output directory.
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="""The output shell script the list of commands will be written.""",
    )
    # Define the argument for specifying the output directory for the arcsiextractdata.py command.
    parser.add_argument(
        "-d",
        "--outDIR",
        type=str,
        required=True,
        help="""The output directory to which arcsiextractdata.py will write outputs.""",
    )
    parser.add_argument(
        "--nofolders",
        action="store_true",
        default=False,
        help="""Specifies individual folders should not be
                                created for each archive which is being extracted.""",
    )
    # Call the parser to parse the arguments.
    args = parser.parse_args()

    if (args.input == None) and (args.listfile == None):
        print("An input directory (-i) or file (-l) must be specified.")
        sys.exit()

    arcsiObj = ARCSIBuildExtractDataCmds()
    if not args.input is None:
        arcsiObj.runGenCommandsDIR(args.input, args.output, args.outDIR, args.nofolders)
    elif not args.listfile is None:
        arcsiObj.runGenCommandsFiles(
            args.listfile, args.output, args.outDIR, args.nofolders
        )
