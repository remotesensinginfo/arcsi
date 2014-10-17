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


class ARCSISortLandsatData (object):
    
    def createDIRStruct(self, outputDIR, noDIRStruct):    
        if not os.path.exists(outputDIR):
            os.makedirs(outputDIR)
            if not noDIRStruct:
                os.makedirs(os.path.join(outputDIR, "RAW"))
                os.makedirs(os.path.join(outputDIR, "Inputs"))
                os.makedirs(os.path.join(outputDIR, "Outputs"))
                os.makedirs(os.path.join(outputDIR, "tmp"))
        
    def run(self, inputDir, outputDir, noDIRStruct):
        inputDir = os.path.abspath(inputDir)
        outputDir = os.path.abspath(outputDir)
        
        if not os.path.isdir(inputDir):
            raise ARCSIException("The input directory specified does not exist!")
        if not os.path.isdir(outputDir):
            raise ARCSIException("The output directory specified does not exist!")
        
        inputFiles = os.listdir(inputDir)
        
        createdLM1DIR = False
        createdLM2DIR = False
        createdLM3DIR = False
        createdLM4DIR = False
        createdLS4DIR = False
        createdLM5DIR = False
        createdLS5DIR = False
        createdLS7DIR = False
        createdLS8DIR = False
        outputFileDIR = ""
        
        for file in inputFiles:
            filePrefix = file[:3]
            if filePrefix == 'LM1':
                outputFileDIR = os.path.join(outputDir, "LM1")
                if (not createdLM1DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct)
                    createdLM1DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                shutil.move(inFile, outputFileDIR)
            elif filePrefix == 'LM2':
                outputFileDIR = os.path.join(outputDir, "LM2")
                if (not createdLM2DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct)
                    createdLM2DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                shutil.move(inFile, outputFileDIR)
            elif filePrefix == 'LM3':
                outputFileDIR = os.path.join(outputDir, "LM3")
                if (not createdLM3DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct)
                    createdLM3DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                shutil.move(inFile, outputFileDIR)
            elif filePrefix == 'LM4':
                outputFileDIR = os.path.join(outputDir, "LM4")
                if (not createdLM4DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct)
                    createdLM4DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                shutil.move(inFile, outputFileDIR)
            elif filePrefix == 'LM5':
                outputFileDIR = os.path.join(outputDir, "LM5")
                if (not createdLM5DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct)
                    createdLM5DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                shutil.move(inFile, outputFileDIR)
            elif filePrefix == 'LT4':
                outputFileDIR = os.path.join(outputDir, "LS4")
                if (not createdLS4DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct)
                    createdLS4DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                shutil.move(inFile, outputFileDIR)
            elif filePrefix == 'LT5':
                outputFileDIR = os.path.join(outputDir, "LS5")
                if (not createdLS5DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct)
                    createdLS5DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                shutil.move(inFile, outputFileDIR)
            elif filePrefix == 'LE7':
                outputFileDIR = os.path.join(outputDir, "LS7")
                if (not createdLS7DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct)
                    createdLS7DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                shutil.move(inFile, outputFileDIR)
            elif filePrefix == 'LC8':
                outputFileDIR = os.path.join(outputDir, "LS8")
                if (not createdLS8DIR) and (not os.path.isdir(outputFileDIR)):
                    self.createDIRStruct(outputFileDIR, noDIRStruct)
                    createdLS8DIR = True
                # Move file...
                inFile = os.path.join(inputDir, file)
                if not noDIRStruct:
                    outputFileDIR = os.path.join(outputFileDIR, "RAW")
                print("Moving: " + inFile)
                print("To: " + outputFileDIR)
                shutil.move(inFile, outputFileDIR)
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
    parser.add_argument('-v', '--version', action='version', version='%(prog)s Version 0.9.0')
    # Define the argument for specifying the input spectral response file.
    parser.add_argument("-i", "--input", type=str, 
                        help='''Input directory containing the input Landsat Scenes.''')
    # Define the argument for specifying input seperator.
    parser.add_argument("-o", "--output", type=str,
                        help='''The output directory to which the output structure will be written.''')
    parser.add_argument("--nodirstruct", action='store_true', default=False, 
                        help='''Specifies that a directory structure should not be built when the new folders are created.''')
    # Call the parser to parse the arguments.
    args = parser.parse_args()
    
    if args.input == None:
        print("An input directory was not specified.")
        parser.print_help()
        sys.exit()
    
    if args.output == None:
        print("An output directory was not specified.")
        parser.print_help()
        sys.exit()
    
    arcsiObj = ARCSISortLandsatData()
    try:
        arcsiObj.run(args.input, args.output, args.nodirstruct)
    except ARCSIException as e:
        print("Error: " + str(e))
    except Exception as e:
        print("Error: " + str(e))

