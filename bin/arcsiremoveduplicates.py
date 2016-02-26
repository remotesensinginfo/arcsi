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

# Import the python os.path module
import os.path
# Import the python sys module
import sys
# Import the python Argument parser
import argparse
# Import the arcsi version number
from arcsilib import ARCSI_VERSION
# Import os.walk to navigate directory structure.
import os
# Import JSON module
import json
# Import shutil module
import shutil

class ARCSIRemoveDuplicates (object):
    
    def findBaseDIR(self, filePath):
        baseName = os.path.dirname(filePath)
        found = False
        while not found:
            baseNameTmp = os.path.dirname(baseName)
            if baseNameTmp == baseName:
                found = True
                break
            elif (baseNameTmp == "") or (baseNameTmp == "/"):
                found = True
                break
            else:
                baseName = baseNameTmp
        return baseName
    
    
    def sortDuplicateFiles(self, lutFile, headersDIR, archivesDIR, cpArchives2DIR):
        headersDIR = os.path.abspath(headersDIR)
        archivesDIR = os.path.abspath(archivesDIR)
        cpArchives2DIR = os.path.abspath(cpArchives2DIR)
        if not archivesDIR is None:
            if cpArchives2DIR is None:
                print("If archives directory is specified then the output directory to which the duplicate archives will be moved must be specified.")
                sys.exit()
        with open(lutFile, 'r') as f:
            jsonStrData = f.read()
        fileLUT = json.loads(jsonStrData)

        for fileBase in fileLUT:
            print("Processing: " + fileBase)
            if 'Duplicates' in fileLUT[fileBase]:
                for dup in fileLUT[fileBase]['Duplicates']:
                    if not headersDIR is None:
                        archFile = os.path.join(archivesDIR, dup['Archive'])
                        archOutFile = os.path.join(cpArchives2DIR, dup['Archive'])
                        if os.path.isfile(archFile):
                            print("Moving: " + archFile)
                            shutil.move(archFile, archOutFile)
                    if not headersDIR is None:
                        dirName = self.findBaseDIR(dup['Header'])
                        fullPath2Del = os.path.join(headersDIR, dirName)
                        if os.path.exists(fullPath2Del):
                            print("Deleting: " + fullPath2Del)
                            shutil.rmtree(fullPath2Del)
            print("")

if __name__ == '__main__':
    """
    The command line user interface to ARCSI.
    """
    parser = argparse.ArgumentParser(prog='arcsiremoveduplicates.py',
                                    description='''ARCSI command to remove duplicate files.''',
                                    epilog='''A tool to remove duplicate files using the LUT 
                                    file generated from the arcsibuildfilenameslu.py command.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)
    
    parser.add_argument("-l", "--lut", type=str,  required=True,
                        help='''Look up table (as generated arcsibuildfilenameslu.py) of the arcsi file names longside input headers''')
    
    parser.add_argument("-w", "--workingdir", type=str, 
                        help='''Working directory from which the header file paths are reference from in the LUT''')

    parser.add_argument("-a", "--archivedir", type=str, 
                        help='''Directory containing the archvies references in the LUT.''')
                        
    parser.add_argument("-d", "--dirout", type=str, 
                        help='''A directory where the archives will be moved to.''')
                        
    
    # Call the parser to parse the arguments.
    args = parser.parse_args()
    
    if (args.workingdir is None) and (args.archivedir is None):
        print("Either or both the working and archives directories need to be defined")
        sys.exit()
    
    if not args.archivedir is None:
        if args.dirout is None:
            print("If archives directory is specified then the output directory to which the duplicate archives will be moved must be specified.")
            sys.exit()
            
    arcsiObj = ARCSIRemoveDuplicates()
    arcsiObj.sortDuplicateFiles(args.lut, args.workingdir, args.archivedir, args.dirout)
    
    
