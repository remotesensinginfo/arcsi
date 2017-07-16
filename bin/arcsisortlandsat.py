#! /usr/bin/env python

"""
A command to sort the landsat scenes into directories for the specific
sensors (i.e., Landsat 1, Landsat 2 ... Landsat 8 etc).
"""


############################################################################
#  arcsisortlandsat.py
#
#  Copyright 2014 ARCSI.
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
# Purpose:  A script to sort Landsat scenes into sensor specific
#           directories.
#
#           Landsat 1 MSS: LM1
#           Landsat 2 MSS: LM2
#           Landsat 3 MSS: LM3
#           Landsat 4 MSS: LM4
#           Landsat 4 TM:  LS4
#           Landsat 5 MSS: LM5
#           Landsat 5 TM:  LS5
#           Landsat 7 TM:  LS7
#           Landsat 8 TM:  LS8
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 22/07/2014
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
# Import the python os module
import os
# Import the python sys module
import sys
# Import the python Argument parser
import argparse
# Import the ARCSI exception class
from arcsilib.arcsiexception import ARCSIException
# Import the shutil python module.
import shutil
#
import os.path
# Import the arcsi version number
from arcsilib import ARCSI_VERSION

class ARCSISortLandsatData (object):

    def createDIRStruct(self, outputDIR, noDIRStruct, inputFiles):
        if not os.path.exists(outputDIR):
            os.makedirs(outputDIR)
            if not noDIRStruct:
                if inputFiles:
                    os.makedirs(os.path.join(outputDIR, "RAW"))
                os.makedirs(os.path.join(outputDIR, "Inputs"))
                os.makedirs(os.path.join(outputDIR, "Outputs"))
                os.makedirs(os.path.join(outputDIR, "tmp"))


    def moveFile(self, cFileLoc, outFileDIR, userInteract):
        ## CHECK IF FILE EXISTS AT outFileDIR
        ## IF EXISTS and userInteract == True THEN ask users which to use (show file sizes)
        ## ELSE (i.e., userInteract == False) THEN don't overwrite and leave orignal file where it is.
        fileName = os.path.basename(cFileLoc)
        if os.path.exists(os.path.join(outFileDIR,fileName)):
            if userInteract:
                print("1) ", cFileLoc) # DUPLICATE
                print("\t File Size: " + str((os.path.getsize(cFileLoc)/1024)/1024))
                print("2) ", os.path.join(outFileDIR,fileName)) # ORIGINAL
                print("\t File Size: " + str((os.path.getsize(os.path.join(outFileDIR,fileName))/1024)/1024))
                # Print File
                Answer = input("Press 1 to move NEW file or 2 to keep ORIGINAL.\n")
                if Answer == '1': # Move New
                    os.remove(os.path.join(outFileDIR,fileName))
                    shutil.move(cFileLoc, outFileDIR)
                else:
                    pass
            else:
                print("DUPLICATE IGNORING: ", cFileLoc)
        else:
            shutil.move(cFileLoc, outFileDIR)

    def runFiles(self, inputDir, outputDir, noDIRStruct, userInteract):
        inputDir = os.path.abspath(inputDir)
        outputDir = os.path.abspath(outputDir)

        if not os.path.isdir(inputDir):
            raise ARCSIException("The input directory specified does not exist!")
        if not os.path.isdir(outputDir):
            raise ARCSIException("The output directory specified does not exist!")

        #inputFiles = os.listdir(inputDir)

        inputFiles = []

        # Navigate the directory tree
        for dirName, sudirList, fileList in os.walk(inputDir):
            # Append target file to list f files using the absolute filepath
            for fname in fileList:
                filename = os.path.join(dirName, fname)
                inputFiles.append(filename)

        createdLM1DIR = False
        createdLM2DIR = False
        createdLM3DIR = False
        createdLM4DIR = False
        createdLS4DIR = False
        createdLM5DIR = False
        createdLS5DIR = False
        createdLS7DIR = False
        createdLS8DIR = False
        createdLS05DIR = False
        createdLS07DIR = False
        outputFileDIR = ""

        for file in inputFiles:
            #print(file)
            basefilename = os.path.basename(file)
            #print(basefilename)
            filePrefix3 = basefilename[:3]
            filePrefix4 = basefilename[:4]
            #print(filePrefix)
            if '.DS_Store' in file:
                print('Skipping .DS_Store file')
                pass
            elif filePrefix3 == 'LM1' or filePrefix4 == 'LM01':
                outputFileDIR = os.path.join(outputDir, "LM1")
                if (not createdLM1DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, True)
                    createdLM1DIR = True
                # Move file...
                inFile = file
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LM2' or filePrefix4 == 'LM02':
                outputFileDIR = os.path.join(outputDir, "LM2")
                if (not createdLM2DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, True)
                    createdLM2DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LM3' or filePrefix4 == 'LM03':
                outputFileDIR = os.path.join(outputDir, "LM3")
                if (not createdLM3DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, True)
                    createdLM3DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LM4' or filePrefix4 == 'LM04':
                outputFileDIR = os.path.join(outputDir, "LM4")
                if (not createdLM4DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, True)
                    createdLM4DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LM5' or filePrefix4 == 'LM05':
                outputFileDIR = os.path.join(outputDir, "LM5")
                if (not createdLM5DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, True)
                    createdLM5DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LT4' or filePrefix4 == 'LS04' or filePrefix4 == 'LE04' or filePrefix4 == 'LT04':
                outputFileDIR = os.path.join(outputDir, "LS4")
                if (not createdLS4DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, True)
                    createdLS4DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LT5' or filePrefix4 == 'LS05' or filePrefix4 == 'LE05' or filePrefix4 == 'LT05':
                outputFileDIR = os.path.join(outputDir, "LS5")
                if (not createdLS5DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, True)
                    createdLS5DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LE7' or filePrefix4 == 'LS07' or filePrefix4 == 'LE07' or filePrefix4 == 'LT07':
                outputFileDIR = os.path.join(outputDir, "LS7")
                if (not createdLS7DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, True)
                    createdLS7DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LC8' or filePrefix4 == 'LS08' or filePrefix4 == 'LC08':
                outputFileDIR = os.path.join(outputDir, "LS8")
                if (not createdLS8DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, True)
                    createdLS8DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            else:
                print("Sensor was not recognised for file: " + file)

    def runFolders(self, inputDir, outputDir, noDIRStruct, userInteract):
        inputDir = os.path.abspath(inputDir)
        outputDir = os.path.abspath(outputDir)

        if not os.path.isdir(inputDir):
            raise ARCSIException("The input directory specified does not exist!")
        if not os.path.isdir(outputDir):
            raise ARCSIException("The output directory specified does not exist!")

        inputTmpFiles = os.listdir(inputDir)

        inputFiles = []

        # Navigate the directory tree
        for fileName in inputTmpFiles:
            filePath = os.path.join(inputDir, fileName)
            if os.path.isdir(filePath):
                inputFiles.append(filePath)

        createdLM1DIR = False
        createdLM2DIR = False
        createdLM3DIR = False
        createdLM4DIR = False
        createdLS4DIR = False
        createdLM5DIR = False
        createdLS5DIR = False
        createdLS7DIR = False
        createdLS8DIR = False
        createdLS05DIR = False
        createdLS07DIR = False
        outputFileDIR = ""

        for file in inputFiles:
            #print(file)
            basefilename = os.path.basename(file)
            #print(basefilename)
            filePrefix3 = basefilename[:3]
            filePrefix4 = basefilename[:4]
            #print(filePrefix3)
            #print(filePrefix4)
            if '.DS_Store' in file:
                print('Skipping .DS_Store file')
                pass
            elif filePrefix3 == 'LM1' or filePrefix4 == 'LM01':
                outputFileDIR = os.path.join(outputDir, "LM1")
                if (not createdLM1DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, False)
                    createdLM1DIR = True
                # Move file...
                inFile = file
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "Inputs")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LM2' or filePrefix4 == 'LM02':
                outputFileDIR = os.path.join(outputDir, "LM2")
                if (not createdLM2DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, False)
                    createdLM2DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "Inputs")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LM3' or filePrefix4 == 'LM03':
                outputFileDIR = os.path.join(outputDir, "LM3")
                if (not createdLM3DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, False)
                    createdLM3DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "Inputs")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LM4' or filePrefix4 == 'LM04':
                outputFileDIR = os.path.join(outputDir, "LM4")
                if (not createdLM4DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, False)
                    createdLM4DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "Inputs")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LM5' or filePrefix4 == 'LM05':
                outputFileDIR = os.path.join(outputDir, "LM5")
                if (not createdLM5DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, False)
                    createdLM5DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "Inputs")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LT4' or filePrefix4 == 'LS04' or filePrefix4 == 'LE04' or filePrefix4 == 'LT04':
                outputFileDIR = os.path.join(outputDir, "LS4")
                if (not createdLS4DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, False)
                    createdLS4DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "Inputs")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LT5' or filePrefix4 == 'LS05' or filePrefix4 == 'LE05' or filePrefix4 == 'LT05':
                outputFileDIR = os.path.join(outputDir, "LS5")
                if (not createdLS5DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, False)
                    createdLS5DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "Inputs")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LE7' or filePrefix4 == 'LS07' or filePrefix4 == 'LE07' or filePrefix4 == 'LT07':
                outputFileDIR = os.path.join(outputDir, "LS7")
                if (not createdLS7DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, False)
                    createdLS7DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "Inputs")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            elif filePrefix3 == 'LC8' or filePrefix4 == 'LS08' or filePrefix4 == 'LC08':
                outputFileDIR = os.path.join(outputDir, "LS8")
                if (not createdLS8DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct, False)
                    createdLS8DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "Inputs")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                self.moveFile(inFile, outputFileDIR, userInteract)
            else:
                print("Sensor was not recognised for file: " + file)

if __name__ == '__main__':
    """
    The command line user interface to ARCSI Landsat Sort Tool.
    """
    parser = argparse.ArgumentParser(prog='arcsisortlandsat.py',
                                    description='''ARCSI command to sort
                                                   Landsat data into a directory
                                                   structure.''',
                                    epilog='''ARCSI needs to process the Landsat
                                              scenes from the different sensors
                                              independently, this command sorts a
                                              directory of input data into different
                                              directories.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)
    # Define the argument for specifying the input spectral response file.
    parser.add_argument("-i", "--input", type=str, required=True,
                        help='''Input directory containing the input Landsat Scenes.''')
    # Define the argument for specifying input seperator.
    parser.add_argument("-o", "--output", type=str, required=True,
                        help='''The output directory to which the output structure will be written.''')
    parser.add_argument("--nodirstruct", action='store_true', default=False,
                        help='''Specifies that a directory structure should not be built when the new folders are created.''')
    parser.add_argument("--userinteract", action='store_true', default=False,
                        help='''Specifies whether the user should be promoted for decision if two files of same name exist.''')
    parser.add_argument("--inputdirs", action='store_true', default=False,
                        help='''Specifies that the inputs are directories and not archive files.''')

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    arcsiObj = ARCSISortLandsatData()
    try:
        if args.inputdirs:
            arcsiObj.runFolders(args.input, args.output, args.nodirstruct, args.userinteract)
        else:
            arcsiObj.runFiles(args.input, args.output, args.nodirstruct, args.userinteract)
    except ARCSIException as e:
        print("Error: " + str(e))
    except Exception as e:
        print("Error: " + str(e))

