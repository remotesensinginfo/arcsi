"""
Module that contains the ARCSIUtils class.
"""
############################################################################
#  arcsiutils.py
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
# Purpose:  A class with some useful utilites.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 05/07/2013
# Version: 1.0
#
# History:
# Version 1.0 - Created.
#
############################################################################

# Import the ARCSI exception class
from .arcsiexception import ARCSIException
# Import the OS python module
import os

class ARCSIUtils (object):
    """
    A class with useful utilties for the ARCSI System.
    """
    
    def getFileExtension(self, format):
        ext = ".NA"
        if format.lower() == "kea":
            ext = ".kea"
        elif format.lower() == "gtiff":
            ext = ".tif"
        elif format.lower() == "hfa":
            ext = ".img"
        elif format.lower() == "envi":
            ext = ".env"
        else:
            raise ARCSIException("The extension for the format specified is unknown.")
        return ext
        
    def readTextFile(self, file):
        """
        Read a text file into a single string
        removing new lines.
        """
        txtStr = ""
        try:
            dataFile = open(file, 'r')
            for line in dataFile:
                txtStr += line.strip()
            dataFile.close()
        except Exception as e:
            raise e
        return txtStr
    
    def stringTokenizer(self, line, delimiter):
        tokens = list()
        token = str()
        for i in range(len(line)):
            if line[i] == delimiter:
                tokens.append(token)
                token = str()
            else:
                token = token + line[i]
        tokens.append(token)
        return tokens
        
    def readSpectralResponseFunc(self, inFile, seperator, ignoreLines, waveCol, respCol):
        specResp = list()
        try:
            specFile = open(inFile, 'r')
            lineCount = 0
            for line in specFile:
                if lineCount >= ignoreLines:
                    line = line.strip()
                    if line:
                        lineVals = line.split(seperator)
                        if (len(lineVals) <= waveCol) or (len(lineVals) <= respCol):
                            raise ARCSIException("")
                        waveVal = float(lineVals[waveCol].strip())
                        respVal = float(lineVals[respCol].strip())
                        specResp.append([waveVal,respVal])
                lineCount += 1
            specFile.close()
        except ARCSIException as e:
            raise e
        except Exception as e:
            raise e
        return numpy.array(specResp)
        
    
    
    def isSummerOrWinter(self, lat, long, date):
        summerWinter = 0
        if lat < 0:
            # Southern Hemisphere
            if (date.month > 4) & (date.month < 11):
                summerWinter = 2 # Winter
            else:
                summerWinter = 3 # Summer
        else: 
            # Northern Hemisphere
            if (date.month > 3) & (date.month < 10):
                summerWinter = 1 # summer
            else:
                summerWinter = 2 # Winter
        return summerWinter
        
    def getEnvironmentVariable(self, var):
    	outVar = None
    	try:
    		outVar = os.environ[var]
    		print(outVar)
    	except Exception:
    		outVar = None
    	return outVar
    	

