#! /usr/bin/env python

"""
Module that contains the ARSCI command to generate the py6S function call
to define a spectral response.
"""

############################################################################
#  arcsicreatepy6scall.py
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
# Purpose:  A script to take a spectral response text file and to generate
#           the associated py6S function call. This is most likely a utility
#           which will only be used by anyone adding new sensors to ARCSI.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 03/02/2017
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
# Import the python Argument parser
import argparse
# Import the arcsi version number
from arcsilib import ARCSI_VERSION
# import RSGISLib module
import rsgislib

def createPy6SSensorFuncCall(inputFile, outputFile, wvLenInNM):
    """
    * inputFile - 
    * outputFile - 
    """

    rsgisUtils = rsgislib.RSGISPyUtils()
    lineList = rsgisUtils.readTextFile2List(inputFile)

    startWVLen = ''
    endWVLen = ''
    respStr = ''
    first = True
    for line in lineList:
        splitLine = line.split(',')
        if len(splitLine) == 2:
            if first:
                startWVLen = splitLine[0]
                endWVLen = splitLine[0]
                respStr = splitLine[1]
                first = False
            else:
                endWVLen = splitLine[0]
                respStr = respStr + ', ' + splitLine[1]

    if wvLenInNM:
        startWVLen = str(float(startWVLen)/1000)
        endWVLen = str(float(endWVLen)/1000)

    cmd = 's.wavelength = Py6S.Wavelength('+startWVLen+', '+endWVLen+', ['+respStr+'])'
    #print(cmd)
    rsgisUtils.writeList2File([cmd], outputFile)

if __name__ == '__main__':
    """
    The command line user interface to ARCSI.
    """
    parser = argparse.ArgumentParser(prog='arcsicreatepy6scall.py',
                                    description='''ARCSI command create py6S function calls for creating sensor response functions.''')
    # Request the version number.
    parser.add_argument('-v', '--version', action='version', version='%(prog)s version ' + ARCSI_VERSION)

    parser.add_argument("-i", "--input", type=str, required=True,
                        help='''Input file (i.e., output from arcsispecresponsefuncs.py)''')

    parser.add_argument("-o", "--output", type=str, required=True,
                        help='''Output text with the py6S function''')

    parser.add_argument("--nm", action='store_true', default=False,
                    help='''Specify that the input wavelengths are in nm and there need converting''')

    # Call the parser to parse the arguments.
    args = parser.parse_args()

    createPy6SSensorFuncCall(args.input, args.output, args.nm)




