#! /usr/bin/env python

"""
Module that contains the ARSCI command to extract data from archives.
"""

############################################################################
#  arcsiextractdata.py
#
#  Copyright 2013 ARCSI.
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
# Purpose:  A script to unarchive data from tar / tar.gz files into
#           a directory structure with a directory for contents of
#           each archive.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 29/01/2014
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

import argparse
import glob
import os
import sys

import rsgislib.tools.utils

from arcsilib import ARCSI_VERSION


class ARCSIExtractData(object):
    def untargzFiles(self, filelist, outDIR, noFolders):
        tarcommand = "tar -xvzf "
        command = ""
        try:
            processingDIR = ""
            for filename in filelist:
                if noFolders:
                    processingDIR = outDIR
                else:
                    processingDIR = os.path.join(
                        outDIR, os.path.basename(filename).split(".")[0]
                    )
                if not os.path.exists(processingDIR):
                    os.makedirs(processingDIR)
                os.chdir(processingDIR)
                print("Extracting: ", filename)
                print("Output to: ", processingDIR)
                command = tarcommand + '"' + filename + '"'
                os.system(command)
        except Exception as e:
            print("IOError Occurred: " + str(e))

    def untargzFile(self, filename, outDIR, noFolders):
        tarcommand = "tar -xvzf "
        command = ""
        try:
            processingDIR = ""
            if noFolders:
                processingDIR = outDIR
            else:
                processingDIR = os.path.join(
                    outDIR, os.path.basename(filename).split(".")[0]
                )
            if not os.path.exists(processingDIR):
                os.makedirs(processingDIR)
            os.chdir(processingDIR)
            print("Extracting: ", filename)
            print("Output to: ", processingDIR)
            command = tarcommand + '"' + filename + '"'
            os.system(command)

        except Exception as e:
            print("IOError Occurred: " + str(e))

    def untarFiles(self, filelist, outDIR, noFolders):
        tarcommand = "tar -xvf "
        command = ""
        try:
            processingDIR = ""
            for filename in filelist:
                if noFolders:
                    processingDIR = outDIR
                else:
                    processingDIR = os.path.join(
                        outDIR, os.path.basename(filename).split(".")[0]
                    )
                if not os.path.exists(processingDIR):
                    os.makedirs(processingDIR)
                os.chdir(processingDIR)
                print(filename)
                command = tarcommand + '"' + filename + '"'
                os.system(command)
        except Exception as e:
            print("IOError Occurred: " + str(e))

    def untarFile(self, filename, outDIR, noFolders):
        tarcommand = "tar -xvf "
        command = ""
        try:
            processingDIR = ""
            if noFolders:
                processingDIR = outDIR
            else:
                processingDIR = os.path.join(
                    outDIR, os.path.basename(filename).split(".")[0]
                )
            if not os.path.exists(processingDIR):
                os.makedirs(processingDIR)
            os.chdir(processingDIR)
            print(filename)
            command = tarcommand + '"' + filename + '"'
            os.system(command)
        except Exception as e:
            print("IOError Occurred: " + str(e))

    def unZipFiles(self, filelist, outDIR, noFolders):
        zipcommand = "unzip "
        command = ""
        try:
            processingDIR = ""
            for filename in filelist:
                if noFolders:
                    processingDIR = outDIR
                else:
                    processingDIR = os.path.join(
                        outDIR, os.path.basename(filename).split(".")[0]
                    )
                if not os.path.exists(processingDIR):
                    os.makedirs(processingDIR)
                os.chdir(processingDIR)
                print(filename)
                command = zipcommand + '"' + filename + '"'
                os.system(command)
        except Exception as e:
            print("IOError Occurred: " + str(e))

    def unZipFile(self, filename, outDIR, noFolders):
        zipcommand = "unzip "
        command = ""
        try:
            processingDIR = ""
            if noFolders:
                processingDIR = outDIR
            else:
                processingDIR = os.path.join(
                    outDIR, os.path.basename(filename).split(".")[0]
                )
            if not os.path.exists(processingDIR):
                os.makedirs(processingDIR)
            os.chdir(processingDIR)
            print(filename)
            command = zipcommand + '"' + filename + '"'
            os.system(command)
        except Exception as e:
            print("IOError Occurred: " + str(e))

    def untarbzFiles(self, filelist, outDIR, noFolders):
        tarcommand = "tar -xvjf "
        command = ""
        try:
            processingDIR = ""
            for filename in filelist:
                if noFolders:
                    processingDIR = outDIR
                else:
                    processingDIR = os.path.join(
                        outDIR, os.path.basename(filename).split(".")[0]
                    )
                if not os.path.exists(processingDIR):
                    os.makedirs(processingDIR)
                os.chdir(processingDIR)
                print("Extracting: ", filename)
                print("Output to: ", processingDIR)
                command = tarcommand + '"' + filename + '"'
                os.system(command)
        except Exception as e:
            print("IOError Occurred: " + str(e))

    def untarbzFile(self, filename, outDIR, noFolders):
        tarcommand = "tar -xvjf "
        command = ""
        try:
            processingDIR = ""
            if noFolders:
                processingDIR = outDIR
            else:
                processingDIR = os.path.join(
                    outDIR, os.path.basename(filename).split(".")[0]
                )
            if not os.path.exists(processingDIR):
                os.makedirs(processingDIR)
            os.chdir(processingDIR)
            print("Extracting: ", filename)
            print("Output to: ", processingDIR)
            command = tarcommand + '"' + filename + '"'
            os.system(command)

        except Exception as e:
            print("IOError Occurred: " + str(e))

    def run4DIR(self, inputDIR, outputDIR, noFolders):
        inputDIR = os.path.abspath(inputDIR)
        outputDIR = os.path.abspath(outputDIR)

        # First, do the tar.gz files.
        inputFileListTarGz = glob.glob(os.path.join(inputDIR, "*.tar.gz"))
        self.untargzFiles(inputFileListTarGz, outputDIR, noFolders)
        inputFileListTarGz = glob.glob(os.path.join(inputDIR, "*.tgz"))
        self.untargzFiles(inputFileListTarGz, outputDIR, noFolders)
        inputFileListTarGz = glob.glob(os.path.join(inputDIR, "*.TAR.GZ"))
        self.untargzFiles(inputFileListTarGz, outputDIR, noFolders)
        inputFileListTarGz = glob.glob(os.path.join(inputDIR, "*.TGZ"))
        self.untargzFiles(inputFileListTarGz, outputDIR, noFolders)

        # Second, do the tar files.
        inputFileListTar = glob.glob(os.path.join(inputDIR, "*.tar"))
        self.untarFiles(inputFileListTar, outputDIR, noFolders)
        inputFileListTar = glob.glob(os.path.join(inputDIR, "*.TAR"))
        self.untarFiles(inputFileListTar, outputDIR, noFolders)

        # Third, do the zip files.
        inputFileListZip = glob.glob(os.path.join(inputDIR, "*.zip"))
        self.unZipFiles(inputFileListZip, outputDIR, noFolders)
        inputFileListZip = glob.glob(os.path.join(inputDIR, "*.ZIP"))
        self.unZipFiles(inputFileListZip, outputDIR, noFolders)

        # Fourth, do the tar.bz files
        inputFileListTarBz = glob.glob(os.path.join(inputDIR, "*.tar.bz"))
        self.untarbzFiles(inputFileListTarBz, outputDIR, noFolders)
        inputFileListTarBz = glob.glob(os.path.join(inputDIR, "*.TAR.BZ"))
        self.untarbzFiles(inputFileListTarBz, outputDIR, noFolders)
        inputFileListTarBz = glob.glob(os.path.join(inputDIR, "*.tar.bz2"))
        self.untarbzFiles(inputFileListTarBz, outputDIR, noFolders)
        inputFileListTarBz = glob.glob(os.path.join(inputDIR, "*.TAR.BZ2"))
        self.untarbzFiles(inputFileListTarBz, outputDIR, noFolders)

    def run4File(self, inputFile, outputDIR, noFolders):
        inputFile = os.path.abspath(inputFile)
        outputDIR = os.path.abspath(outputDIR)

        fileExt = os.path.basename(inputFile).split(".", 1)[-1].lower()

        if (fileExt == "tar.gz") or (fileExt == "tgz"):
            self.untargzFile(inputFile, outputDIR, noFolders)

        if fileExt == "tar":
            self.untarFile(inputFile, outputDIR, noFolders)

        if fileExt == "zip":
            self.unZipFile(inputFile, outputDIR, noFolders)

        if fileExt == "tar.bz":
            self.untarbzFile(inputFile, outputDIR, noFolders)

        if fileExt == "tar.bz2":
            self.untarbzFile(inputFile, outputDIR, noFolders)

    def run4List(self, inputListFile, outputDIR, noFolders):

        archsList = rsgislib.tools.utils.read_text_file_no_new_lines2List(inputListFile)
        for archFile in archsList:
            self.run4File(archFile, outputDIR, noFolders)


if __name__ == "__main__":
    """
    The command line user interface to ARCSI Data Extraction Tool.
    """
    parser = argparse.ArgumentParser(
        prog="arcsiextractdata.py",
        description="""ARCSI command extract data
                                                   from tar or tar.gz archives.""",
        epilog="""A tools to extract data
                                              from tar or tar.gz archives into
                                              individual directories per image""",
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
        help="""Input directory contains archives (tar, tar.gz, tar.bz, tar.bz2 and/or zip).""",
    )
    # Define the argument for specifying the input file to be processed
    parser.add_argument(
        "-f",
        "--file",
        type=str,
        help="""Input file contains archive (tar, tar.gz, tar.bz, tar.bz2 and/or zip).""",
    )
    # Define the argument for specifying a list of input files to be processed
    parser.add_argument(
        "-l",
        "--list",
        type=str,
        help="""Input file contains archive (tar, tar.gz, tar.bz, tar.bz2 and/or zip).""",
    )
    # Define the argument for specifying the output directory.
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="""The output directory to which all output files are to be written.""",
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

    if (args.input == None) & (args.file == None) & (args.list == None):
        print(
            "Error: An input directory, list as a text file or single archive file must be specified."
        )
        sys.exit()

    arcsiObj = ARCSIExtractData()

    if args.input is not None:
        arcsiObj.run4DIR(args.input, args.output, args.nofolders)

    if args.file is not None:
        arcsiObj.run4File(args.file, args.output, args.nofolders)

    if args.list is not None:
        arcsiObj.run4List(args.list, args.output, args.nofolders)
