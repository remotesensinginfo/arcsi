#! /usr/bin/env python

"""
Module that contains the ARSCI command to extract data from archives.
"""

############################################################################
#  arcsi.py
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

# Import the python os.path module
import os.path
# Import the python sys module
import sys
# Import the python glob module
import glob
# Import the python Argument parser
import argparse

class ARCSIExtractData (object):
    
    def untargzFiles(self, filelist, outDIR, noFolders):
        tarcommand = 'tar -xvzf '
        command = ''
        try:
            processingDIR = ""
            for filename in filelist:
                if noFolders:
                    processingDIR = outDIR
                else:
                    processingDIR = os.path.join(outDIR, os.path.basename(filename).split(".")[0])
                if not os.path.exists(processingDIR):
                    os.makedirs(processingDIR)
                os.chdir(processingDIR)
                print("Extracting: ", filename)
                print("Output to: ", processingDIR)
                command = tarcommand + filename
                os.system(command)
        except Exception as e:
            print('IOError Occurred: ' + str(e))
            
            
    def untarFiles(self, filelist, outDIR, noFolders):
        tarcommand = 'tar -xvf '
        command = ''
        try:
            processingDIR = ""
            for filename in filelist:
                if noFolders:
                    processingDIR = outDIR
                else:
                    processingDIR = os.path.join(outDIR, os.path.basename(filename).split(".")[0])
                if not os.path.exists(processingDIR):
                    os.makedirs(processingDIR)
                os.chdir(processingDIR)
                print(filename)
                command = tarcommand + filename
                os.system(command)
        except Exception as e:
            print('IOError Occurred: ' + str(e))
            
    
    def run(self, inputDIR, outputDIR, noFolders):
        # First, do the tar.gz files.
        inputFileListTarGz = glob.glob(os.path.join(inputDIR, "*.tar.gz"))
        #print(inputFileListTarGz)
        self.untargzFiles(inputFileListTarGz, outputDIR, noFolders)
        
        # Second, do the tar files.
        inputFileListTar = glob.glob(os.path.join(inputDIR, "*.tar"))
        #print(inputFileListTar)
        self.untarFiles(inputFileListTar, outputDIR, noFolders)
        

if __name__ == '__main__':
    """
    The command line user interface to ARCSI Data Extraction Tool.
    """
    parser = argparse.ArgumentParser(prog='arcsiextractdata',
                                    description='''ARCSI command extract data
                                                   from tar or tar.gz archives.''',
                                    epilog='''A tools to extract data
                                              from tar or tar.gz archives into
                                              individual directories per image''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s Version 0.1a')
    # Define the argument for specifying the input spectral response file.
    parser.add_argument("-i", "--input", type=str, 
                        help='''Input directory contains archives (tar and/or tar.gz).''')
    # Define the argument for specifying input seperator.
    parser.add_argument("-o", "--output", type=str,
                        help='''The output directory to which all output files are to be written.''')
    parser.add_argument("--nofolders", action='store_true', default=False, 
                        help='''Specifies individual folders should not be
                                created for each archive which is being extracted.''')
    # Call the parser to parse the arguments.
    args = parser.parse_args()
    
    if args.input == None:
        print("An input directory was not specified.")
        sys.exit()
    
    if args.output == None:
        print("An output directory was not specified.")
        sys.exit()
    
    arcsiObj = ARCSIExtractData()
    
    arcsiObj.run(args.input, args.output, args.nofolders)