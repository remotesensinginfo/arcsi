#! /usr/bin/env python

"""
Module that contains the ARSCI command to build a LUT for the output file names.
"""

############################################################################
#  arcsibuildfilenameslu.py
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
# Import the list of sensors arcsi supports
from arcsilib import ARCSI_SENSORS_LIST
# Import the list of archive file extensions arcsi supports
from arcsilib import ARCSI_ARCHIVE_EXE_LIST

class ARCSIBuildFileNameLUT (object):
    
    def getListOfFiles(self, searchDIR, headerEnding):
        outFiles = []
        for dirName, subdirList, fileList in os.walk(searchDIR):
            for fname in fileList:
                fname = str(fname)
                if fname.endswith(headerEnding):
                    outFiles.append(os.path.abspath(os.path.join(dirName, fname)))
        return outFiles
    
    def getListOfArchives(self, dirPath):
        archPaths = []
        for exe in ARCSI_ARCHIVE_EXE_LIST:
            fileList = glob.glob(os.path.join(dirPath, "*" + exe))
            archPaths.extend(fileList)
        return archPaths
    
    def buildLookUp(self, inputDIR, headerEnding, outputFile, sensorStr, archivesDIR=None):
        inputDIR = os.path.abspath(inputDIR)
        hdrList = self.getListOfFiles(inputDIR, headerEnding)
        archLUT = dict()
        if not archivesDIR is None:
            archPaths = self.getListOfArchives(archivesDIR)
            for arch in archPaths:
                archBaseName = os.path.basename(arch).split(".")[0]
                archLUT[archBaseName] = arch
        
        fileDict = dict()
        sensorFact = ARCSISensorFactory()
        duplicate = False
        for fileHdr in hdrList:
            duplicate = False
            sensorClass = sensorFact.getSensorClassFromName(sensorStr, False, None)
            sensorClass.extractHeaderParameters(fileHdr, "")
            outBaseName = sensorClass.generateOutputBaseName()
            
            if outBaseName in fileDict:
                duplicate = True
            
            fileHdr = fileHdr.replace(inputDIR, "")
            if (fileHdr[0] == '/') or (fileHdr[0] == '\\'):
                fileHdr = fileHdr[1:]
            tmpList = {"Header" : fileHdr}
            if not archivesDIR is None:
                for baseArchName in archLUT:
                    if fileHdr.count(baseArchName) > 0:
                        tmpList['Archive'] = os.path.basename(archLUT[baseArchName])
                        break
            if not duplicate:
                fileDict[outBaseName] = tmpList
            else:
                if not 'Duplicates' in fileDict[outBaseName]:
                    fileDict[outBaseName]['Duplicates'] = []
                fileDict[outBaseName]['Duplicates'].append(tmpList)
            
        with open(outputFile, 'w') as outfile:
            json.dump(fileDict, outfile, sort_keys=True, indent=4, separators=(',', ': '), ensure_ascii=False)
        
        

if __name__ == '__main__':
    """
    The command line user interface to ARCSI.
    """
    parser = argparse.ArgumentParser(prog='arcsibuildfilenameslu.py',
                                    description='''ARCSI command to build LUT of file names vs input files.''',
                                    epilog='''A tools to build LUT of file names vs input headers files and 
                                    optionally the archives. Output is a JSON file.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)
    
    parser.add_argument("-i", "--input", type=str, required=True, 
                        help='''Input directory containing the data to be processed''')

    parser.add_argument("-o", "--output", type=str, required=True, 
                        help='''Output text file (shell script) with the list of commands.''')
                        
    parser.add_argument("-e", "--header", type=str, required=True, 
                        help='''The extension / unquie file ending for the input header files.''')
    
    parser.add_argument("-s", "--sensor", required=True, choices=ARCSI_SENSORS_LIST,  
                        help='''Specify the sensor being processed.''')
                        
    parser.add_argument("-a", "--archives", type=str, 
                        help='''Input directory containing the original archives''')
    
    # Call the parser to parse the arguments.
    args = parser.parse_args()
    
    arcsiObj = ARCSIBuildFileNameLUT()
    
    arcsiObj.buildLookUp(args.input, args.header, args.output, args.sensor, args.archives)
    
    