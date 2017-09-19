#!/usr/bin/env python

############################################################################
#  arcsisplitsen2granules.py
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
# Purpose:  A tool to split the old Sentinel-2 scenes which contain 
#           multiple granules into single granule packages which can
#           then be processed through ARCSI.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 02/06/2017
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
# Import the glob module.
import glob
# Import python system library
import sys
#Import the python file os module
import os
#Import the python file paths module
import os.path
#Import the python file shutil module
import shutil
# Import the datetime module
import datetime
# Import the python Argument parser
import argparse
# Import the ARCSI exception class
from arcsilib.arcsiexception import ARCSIException
# Import the arcsi version number
from arcsilib import ARCSI_VERSION


def splitSen2Granules(inputDIR, outputDIR):
	"""
	Function to split multi granule sentinel-2 scenes into single granule scenes. 
	"""
	if (not os.path.exists(inputDIR)) or (not os.path.isdir(inputDIR)) or (not os.path.splitext(inputDIR)[1] == '.SAFE'):
		raise ARCSIException("Input directory either does not exist or is not a directory")

	if (not os.path.exists(outputDIR)) or (not os.path.isdir(outputDIR)):
		raise ARCSIException("Output directory either does not exist or is not a directory")

	granulesDIR = os.path.join(inputDIR, 'GRANULE')
	if (not os.path.exists(granulesDIR)) or (not os.path.isdir(granulesDIR)):
		raise ARCSIException("Input directory does not have expected structure as cannot find granules directory.")

	sceneBaseName = os.path.splitext(os.path.basename(inputDIR))[0]

	files2Copy = list()
	dirs2Copy = list()
	sceneFiles = os.listdir(inputDIR)
	for file in sceneFiles:
		if (not file == "GRANULE") and (not file == '.DS_Store'):
			if os.path.isdir(os.path.join(inputDIR, file)):
				dirs2Copy.append(os.path.join(inputDIR, file))
			elif os.path.isfile(os.path.join(inputDIR, file)):
				files2Copy.append(os.path.join(inputDIR, file))
	
	granulesDIRLst = os.listdir(granulesDIR)
	for granuleDIR in granulesDIRLst:
		granuleDIRPath = os.path.join(granulesDIR, granuleDIR)
		if os.path.isdir(granuleDIRPath) and ('S2A_OPER_MSI_L1C' in granuleDIR):
			tileID = granuleDIR.split('_')[9]
			print('Processing: ' + tileID)
			outDIRBase = sceneBaseName+'_'+tileID+'.SAFE'
			outDIR = os.path.join(outputDIR, outDIRBase)
			if os.path.exists(outDIR):
				raise ARCSIException("Output directory ("+outDIR+") already exists.")
			os.makedirs(outDIR)
			for file in files2Copy:
				shutil.copy2(file, os.path.join(outDIR, os.path.basename(file)))
			for dirPath in dirs2Copy:
				shutil.copytree(dirPath, os.path.join(outDIR, os.path.basename(dirPath)))
			outGranuleDIR = os.path.join(outDIR, 'GRANULE')
			os.makedirs(outGranuleDIR)
			shutil.copytree(granuleDIRPath, os.path.join(outGranuleDIR, granuleDIR))




if __name__ == '__main__':
    """
    The command line user interface to ARCSI Split Sentinel-2 Granules Tool.
    """
    parser = argparse.ArgumentParser(prog='arcsisplitsen2granules.py',
                                    description='''ARCSI command for splitting Sentinel-2 scenes''',
                                    epilog='''A tool to split the old Sentinel-2 scenes which contain 
                                              multiple granules into single granule packages which can 
                                              then be processed through ARCSI.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)

    parser.add_argument("-i", "--input", type=str, required=True,
                        help='''Specify the directory of the sentinel-2 scene.''')

    parser.add_argument("-o", "--output", type=str, required=True,
                        help='''The output directory for the split scenes''')


    # Call the parser to parse the arguments.
    args = parser.parse_args()

    splitSen2Granules(args.input, args.output)


