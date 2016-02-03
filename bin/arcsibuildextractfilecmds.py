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

class ARCSIBuildExtractDataCmds (object):

    def runGenCommands(self, inputDIR, outputFile, outDIR, nofolders):
        outDIR = os.path.abspath(outDIR)
        outFile = open(outputFile, 'w')
        for dirName, subdirList, fileList in os.walk(inputDIR):
            for fname in fileList:
                fileExt = os.path.basename(fname).split(".", 1)[-1].lower() 
                filePath = os.path.abspath(os.path.join(dirName, fname))
                if (fileExt == 'tar') or (fileExt == 'tgz') or (fileExt == 'tar.gz') or (fileExt == 'zip'):
                    print('\t%s' % filePath)
                    command = 'arcsiextractdata.py -f \"' + filePath + '\"  -o \"' + outDIR + '\"'
                    if nofolders:
                        command = command + ' --nofolders'
                    command = command + '\n'
                    outFile.write(command)
        outFile.close()


if __name__ == '__main__':
    """
    The command line user interface to ARCSI build arcsiextractdata.py commands.
    """
    parser = argparse.ArgumentParser(prog='arcsibuildextractfilecmds.py',
                                    description='''ARCSI command to build commands for
                                                arcsiextractdata.py.''',
                                    epilog='''ARCSI command to build commands for
                                              arcsiextractdata.py for extracting 
                                              single files on a HPC type system.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)
    # Define the argument for specifying the input directory to be processed.
    parser.add_argument("-i", "--input", type=str, 
                        help='''Input directory contains archives (tar, tar.gz, and/or zip).''')
    # Define the argument for specifying the output directory.
    parser.add_argument("-o", "--output", type=str,
                        help='''The output shell script the list of commands will be written.''')
    # Define the argument for specifying the output directory for the arcsiextractdata.py command.
    parser.add_argument("-d", "--outDIR", type=str,
                        help='''The output directory to which arcsiextractdata.py will write outputs.''')
    parser.add_argument("--nofolders", action='store_true', default=False, 
                        help='''Specifies individual folders should not be
                                created for each archive which is being extracted.''')
    # Call the parser to parse the arguments.
    args = parser.parse_args()
    
    if (args.input == None):
        print("An input directory or file must be specified.")
        parser.print_help()
        sys.exit()
    
    if args.output == None:
        print("An output text file was not specified.")
        parser.print_help()
        sys.exit()
        
    if args.outDIR == None:
        print("An output directory was not specified.")
        parser.print_help()
        sys.exit()
    
    arcsiObj = ARCSIBuildExtractDataCmds()
    arcsiObj.runGenCommands(args.input, args.output, args.outDIR, args.nofolders)
    




