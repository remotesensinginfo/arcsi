"""
Module that contains the ARCSISentinel2Sensor class.
"""
############################################################################
#  arcsisensorsentinel2.py
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
# Purpose:  A class for read the Sentinel-2 sensor header file and applying
#           the pre-processing operations within ARCSI to the Sentinel-2
#           datasets.
#
# Author: Pete Bunting
# Email: pfb@aber.ac.uk
# Date: 09/07/2016
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
# import abstract base class stuff
from .arcsisensor import ARCSIAbstractSensor
# Import the ARCSI exception class
from .arcsiexception import ARCSIException
# Import the ARCSI utilities class
from .arcsiutils import ARCSIUtils
# Import the datetime module
import datetime
# Import the GDAL/OGR spatial reference library
from osgeo import osr
from osgeo import ogr
# Import OS path module for manipulating the file system
import os.path
# Import the RSGISLib Image Calibration Module.
import rsgislib.imagecalibration
# Import the RSGISLib Image Calculation Module
import rsgislib.imagecalc
# Import the collections module
import collections
# Import the py6s module for running 6S from python.
import Py6S
# Import the python maths library
import math
# Import python XML Parser
import xml.etree.ElementTree as ET
# Import the numpy module
import numpy
# Import the GDAL python module
import osgeo.gdal as gdal
# Import the RIOS RAT module
from rios import rat
# Import the glob tool
import glob
# Import subprocess module
import subprocess

class ARCSISen2SpectralBandObj(object):
    """
    This is a class to store the information associated with a Sentinel-2 image band.
    """
    def __init__(self, phyBandName=None, bandID=None, imgRes=0.0, wvLenMin=0.0, wvLenMax=0.0, wvLenCen=0.0, respFuncStep=1.0, respFunc=list()):
        self.phyBandName = phyBandName
        self.bandID = bandID
        self.imgRes = imgRes
        self.wvLenMin = wvLenMin
        self.wvLenMax = wvLenMax
        self.wvLenCen = wvLenCen
        self.respFuncStep = respFuncStep
        self.respFunc = respFunc
        self.respFunc6S = None
        self.wvLenMin6S = 0.0
        self.wvLenMax6S = 0.0

class ARCSISentinel2Sensor (ARCSIAbstractSensor):
    """
    A class which represents the Sentinel-2 sensor to read
    header parameters and apply data processing operations.
    """
    def __init__(self, debugMode, inputImage):
        ARCSIAbstractSensor.__init__(self, debugMode, inputImage)
        self.sensor = "SEN2"
        self.sen2FileBaseDIR = ''
        self.processingLevel = ''
        self.processingBaseline = ''
        self.spacecraftName = ''
        self.dataTakeType = ''
        self.orbitNumber = ''
        self.orbitDirection = ''

        self.sen2ImgB01 = ''
        self.sen2ImgB02 = ''
        self.sen2ImgB03 = ''
        self.sen2ImgB04 = ''
        self.sen2ImgB05 = ''
        self.sen2ImgB06 = ''
        self.sen2ImgB07 = ''
        self.sen2ImgB8A = ''
        self.sen2ImgB08 = ''
        self.sen2ImgB09 = ''
        self.sen2ImgB10 = ''
        self.sen2ImgB11 = ''
        self.sen2ImgB12 = ''
        self.sen2ImgTCI = ''

        # File paths/names for resampled images
        self.sen2ImgB02_20m = None
        self.sen2ImgB03_20m = None
        self.sen2ImgB04_20m = None
        self.sen2ImgB08_20m = None
        self.sen2ImgB05_10m = None
        self.sen2ImgB06_10m = None
        self.sen2ImgB07_10m = None
        self.sen2ImgB8A_10m = None
        self.sen2ImgB11_10m = None
        self.sen2ImgB12_10m = None

        self.inNoDataVal = 0
        self.inSatDataVal = 65535
        self.quantificationVal = 10000

        self.physicalGain_B0 = 0.0
        self.physicalGain_B1 = 0.0
        self.physicalGain_B2 = 0.0
        self.physicalGain_B3 = 0.0
        self.physicalGain_B4 = 0.0
        self.physicalGain_B5 = 0.0
        self.physicalGain_B6 = 0.0
        self.physicalGain_B7 = 0.0
        self.physicalGain_B8 = 0.0
        self.physicalGain_B9 = 0.0
        self.physicalGain_B10 = 0.0
        self.physicalGain_B11 = 0.0
        self.physicalGain_B12 = 0.0

        self.earthSunDist_U = 0.0
        self.esun_B0 = 0.0
        self.esun_B1 = 0.0
        self.esun_B2 = 0.0
        self.esun_B3 = 0.0
        self.esun_B4 = 0.0
        self.esun_B5 = 0.0
        self.esun_B6 = 0.0
        self.esun_B7 = 0.0
        self.esun_B8 = 0.0
        self.esun_B9 = 0.0
        self.esun_B10 = 0.0
        self.esun_B11 = 0.0
        self.esun_B12 = 0.0

        self.specBandInfo = dict()

        self.uniqueTileID = ''

        self.resampleTo20m = False

        self.epsgCode = 0

        self.viewAnglesFmaskImg = None

        self.imgIntScaleFactor = 1000

    ################################
    #### CODE FROM PYTHON FMASK ####
    ################################
    @staticmethod
    def makeValueArray(valuesListNode):
        """
        Take a <Values_List> node from the XML, and return an array of the values contained
        within it. This will be a 2-d numpy array of float32 values (should I pass the dtype in??)
        
        """
        valuesList = valuesListNode.findall('VALUES')
        vals = []
        for valNode in valuesList:
            text = valNode.text
            vals.append([numpy.float32(x) for x in text.strip().split()])
        return numpy.array(vals)

    def buildViewAngleArr(self, viewingAngleNodeList, angleName):
        """
        Build up the named viewing angle array from the various detector strips given as
        separate arrays. I don't really understand this, and may need to re-write it once
        I have worked it out......
        
        The angleName is one of 'Zenith' or 'Azimuth'.
        Returns a dictionary of 2-d arrays, keyed by the bandId string. 
        """
        angleArrDict = {}
        for viewingAngleNode in viewingAngleNodeList:
            bandId = viewingAngleNode.attrib['bandId']
            angleNode = viewingAngleNode.find(angleName)
            angleArr = self.makeValueArray(angleNode.find('Values_List'))
            if bandId not in angleArrDict:
                angleArrDict[bandId] = angleArr
            else:
                mask = (~numpy.isnan(angleArr))
                angleArrDict[bandId][mask] = angleArr[mask]
        return angleArrDict
    ################################

    def extractHeaderParameters(self, inputHeader, wktStr):
        """
        Understands and parses the Sentinel-2 xml header file
        """
        try:
            arcsiUtils = ARCSIUtils()
            inputHeader = os.path.abspath(inputHeader)
            self.headerFileName = os.path.split(inputHeader)[1]
            self.sen2FileBaseDIR = os.path.split(inputHeader)[0]
            
            print("Reading header file")
            tree = ET.parse(inputHeader)
            root = tree.getroot()

            generalInfoTag = root.find('{https://psd-14.sentinel2.eo.esa.int/PSD/User_Product_Level-1C.xsd}General_Info')
            if generalInfoTag == None:
                raise ARCSIException("Cannot open top level section \'General_Info\' - is this really a Sentinel-2 image file?")

            productInfoTag = generalInfoTag.find('Product_Info')
            if productInfoTag == None:
                raise ARCSIException("Cannot open \'General_Info\' section \'Product_Info\' - is this really a Sentinel-2 image file?")

            productType = productInfoTag.find('PRODUCT_TYPE').text.strip()
            if (productType == None) or (productType == ''):
                raise ARCSIException("Cannot find the product type - is this really a Sentinel-2 image file?")
            elif productType != 'S2MSI1C':
                raise ARCSIException("Expecting Sentinel-2 product type \'S2MSI1C\'")
            else:
                print("Found Sentinel-2 image file product \'S2MSI1C\'")

            self.processingLevel = productInfoTag.find('PROCESSING_LEVEL').text.strip()
            self.processingBaseline = productInfoTag.find('PROCESSING_BASELINE').text.strip()

            datatakeTag = productInfoTag.find('Datatake')
            if datatakeTag == None:
                raise ARCSIException("Cannot open \'Product_Info\' section \'Datatake\'")

            acqTimeStr = datatakeTag.find('DATATAKE_SENSING_START').text.strip()
            acqTimeStr = acqTimeStr.replace('Z', '')
            self.acquisitionTime = datetime.datetime.strptime(acqTimeStr, "%Y-%m-%dT%H:%M:%S.%f")

            self.spacecraftName = datatakeTag.find('SPACECRAFT_NAME').text.strip()
            self.dataTakeType = datatakeTag.find('DATATAKE_TYPE').text.strip()
            self.orbitNumber = datatakeTag.find('SENSING_ORBIT_NUMBER').text.strip()
            self.orbitDirection = datatakeTag.find('SENSING_ORBIT_DIRECTION').text.strip()

            prodURI = productInfoTag.find('PRODUCT_URI').text.strip()
            self.uniqueTileID = prodURI.split('_')[5]


            # Get the input image band file names.
            granuleListTag = productInfoTag.find('Product_Organisation').find('Granule_List')
            granulesTagsLst = list()
            for granuleLstChild in granuleListTag:
                if granuleLstChild.tag == 'Granule':
                    granulesTagsLst.append(granuleLstChild)

            if len(granulesTagsLst) != 1:
                raise ARCSIException("Only expecting a single granule within the file... The input image you have provided is not supported by ARCSI - please report so we can add support.")

            granuleTag = granulesTagsLst[0]
            for granuleChild in granuleTag:
                if granuleChild.tag == 'IMAGE_FILE':
                    imgFile = granuleChild.text.strip()
                    tmpFiles = glob.glob(os.path.join(self.sen2FileBaseDIR, imgFile+'*.jp2'))
                    if len(tmpFiles) == 1:
                        if 'B01' in imgFile:
                            self.sen2ImgB01 = tmpFiles[0]
                        elif 'B02' in imgFile:
                            self.sen2ImgB02 = tmpFiles[0]
                        elif 'B03' in imgFile:
                            self.sen2ImgB03 = tmpFiles[0]
                        elif 'B04' in imgFile:
                            self.sen2ImgB04 = tmpFiles[0]
                        elif 'B05' in imgFile:
                            self.sen2ImgB05 = tmpFiles[0]
                        elif 'B06' in imgFile:
                            self.sen2ImgB06 = tmpFiles[0]
                        elif 'B07' in imgFile:
                            self.sen2ImgB07 = tmpFiles[0]
                        elif 'B8A' in imgFile:
                            self.sen2ImgB8A = tmpFiles[0]
                        elif 'B08' in imgFile:
                            self.sen2ImgB08 = tmpFiles[0]
                        elif 'B09' in imgFile:
                            self.sen2ImgB09 = tmpFiles[0]
                        elif 'B10' in imgFile:
                            self.sen2ImgB10 = tmpFiles[0]
                        elif 'B11' in imgFile:
                            self.sen2ImgB11 = tmpFiles[0]
                        elif 'B12' in imgFile:
                            self.sen2ImgB12 = tmpFiles[0]
                        elif 'TCI' in imgFile:
                            self.sen2ImgTCI = tmpFiles[0]
                        else:
                            raise ARCSIException("Could not associated image file with an expected image band: " + imgFile)
                    else:
                        raise ARCSIException("Could not file image file for: " + imgFile)

            productImgCharTag = generalInfoTag.find('Product_Image_Characteristics')
            if productImgCharTag == None:
                raise ARCSIException("Cannot open \'General_Info\' section \'Product_Image_Characteristics\'")

            specialValsTagsLst = list()
            quantificationValTag = None
            reflectanceConversionTag = None
            spectralInfoListTag = None
            physicalGainsTagsLst = list()
            for prodImgCharChild in productImgCharTag:
                if prodImgCharChild.tag == 'Special_Values':
                    specialValsTagsLst.append(prodImgCharChild)
                elif prodImgCharChild.tag == 'QUANTIFICATION_VALUE':
                    quantificationValTag = prodImgCharChild
                elif prodImgCharChild.tag == 'Reflectance_Conversion':
                    reflectanceConversionTag = prodImgCharChild
                elif prodImgCharChild.tag == 'Spectral_Information_List':
                    spectralInfoListTag = prodImgCharChild
                elif prodImgCharChild.tag == 'PHYSICAL_GAINS':
                    physicalGainsTagsLst.append(prodImgCharChild)

            # Parse out the input no data and saturation value.
            for tag in specialValsTagsLst:
                if tag.find('SPECIAL_VALUE_TEXT').text.strip() == 'NODATA':
                    self.inNoDataVal = arcsiUtils.str2Int(tag.find('SPECIAL_VALUE_INDEX').text.strip())
                elif tag.find('SPECIAL_VALUE_TEXT').text.strip() == 'SATURATED':
                    self.inSatDataVal = arcsiUtils.str2Int(tag.find('SPECIAL_VALUE_INDEX').text.strip())

            # Get the Quantification value.
            self.quantificationVal = arcsiUtils.str2Float(quantificationValTag.text.strip())

            # Get the Physical Gain values.
            for physGainTag in physicalGainsTagsLst:
                if physGainTag.attrib['bandId'] == '0':
                    self.physicalGain_B0 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '1':
                    self.physicalGain_B1 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '2':
                    self.physicalGain_B2 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '3':
                    self.physicalGain_B3 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '4':
                    self.physicalGain_B4 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '5':
                    self.physicalGain_B5 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '6':
                    self.physicalGain_B6 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '7':
                    self.physicalGain_B7 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '8':
                    self.physicalGain_B8 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '9':
                    self.physicalGain_B9 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '10':
                    self.physicalGain_B10 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '11':
                    self.physicalGain_B11 = arcsiUtils.str2Float(physGainTag.text.strip())
                elif physGainTag.attrib['bandId'] == '12':
                    self.physicalGain_B12 = arcsiUtils.str2Float(physGainTag.text.strip())

            # Get Reflectance Conversion info
            self.earthSunDist_U = arcsiUtils.str2Float(reflectanceConversionTag.find('U').text.strip())
            for solarIrrTag in reflectanceConversionTag.find('Solar_Irradiance_List'):
                if solarIrrTag.attrib['bandId'] == '0':
                    self.esun_B0 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '1':
                    self.esun_B1 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '2':
                    self.esun_B2 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '3':
                    self.esun_B3 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '4':
                    self.esun_B4 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '5':
                    self.esun_B5 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '6':
                    self.esun_B6 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '7':
                    self.esun_B7 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '8':
                    self.esun_B8 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '9':
                    self.esun_B9 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '10':
                    self.esun_B10 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '11':
                    self.esun_B11 = arcsiUtils.str2Float(solarIrrTag.text.strip())
                elif solarIrrTag.attrib['bandId'] == '12':
                    self.esun_B12 = arcsiUtils.str2Float(solarIrrTag.text.strip())

            # Get Spectral Info
            for specInfoTag in spectralInfoListTag:
                phyBandName = specInfoTag.attrib['physicalBand']
                specBandObj = ARCSISen2SpectralBandObj()
                specBandObj.phyBandName = phyBandName
                specBandObj.bandID = specInfoTag.attrib['bandId']
                specBandObj.imgRes = arcsiUtils.str2Float(specInfoTag.find('RESOLUTION').text.strip())
                specBandObj.wvLenMin = arcsiUtils.str2Float(specInfoTag.find('Wavelength').find('MIN').text.strip())
                specBandObj.wvLenMax = arcsiUtils.str2Float(specInfoTag.find('Wavelength').find('MAX').text.strip())
                specBandObj.wvLenCen = arcsiUtils.str2Float(specInfoTag.find('Wavelength').find('CENTRAL').text.strip())
                specBandObj.respFuncStep = arcsiUtils.str2Float(specInfoTag.find('Spectral_Response').find('STEP').text.strip())
                specResStrVals = specInfoTag.find('Spectral_Response').find('VALUES').text.strip().split(' ')
                specResVals = list()
                for strVal in specResStrVals:
                    specResVals.append(arcsiUtils.str2Float(strVal))
                specBandObj.respFunc = specResVals
                self.specBandInfo[phyBandName] = specBandObj

            tree = None

            ####### READ GRANULE HEADER FILES ##########
            granuleDIR = os.path.join(self.sen2FileBaseDIR, 'GRANULE')
            granLC1DIRsIn = os.listdir(granuleDIR)
            granLC1DIRs = []
            for dirStr in granLC1DIRsIn:
                if 'L1C' in dirStr:
                    granLC1DIRs.append(dirStr)
            if len(granLC1DIRs) != 1:
                raise ARCSIException("Couldn't find the granule directory")
            granuleHdr = os.path.join(os.path.join(granuleDIR, granLC1DIRs[0]), 'MTD_TL.xml')

            tree = ET.parse(granuleHdr)
            root = tree.getroot()

            # Get Geometric Info Tag
            geometricInfoTag = root.find('{https://psd-12.sentinel2.eo.esa.int/PSD/S2_PDI_Level-1C_Tile_Metadata.xsd}Geometric_Info')

            # Get Tile Geocoding tag.
            tileGeocoding = geometricInfoTag.find('Tile_Geocoding')

            # Get EPSG Code.
            epsgCodeStr = tileGeocoding.find('HORIZONTAL_CS_CODE').text.strip()
            self.epsgCode = arcsiUtils.str2Int(epsgCodeStr.split(':')[1])
            inProj = osr.SpatialReference()
            inProj.ImportFromEPSG(self.epsgCode)
            if self.inWKT is "":
                self.inWKT = inProj.ExportToWkt()

            tlXRes10 = 0.0
            tlYRes10 = 0.0
            xSizeRes10 = 0
            ySizeRes10 = 0
            tlXRes20 = 0.0
            tlYRes20 = 0.0
            xSizeRes20 = 0
            ySizeRes20 = 0
            tlXRes60 = 0.0
            tlYRes60 = 0.0
            xSizeRes60 = 0
            ySizeRes60 = 0
            for tileGeocodeChild in tileGeocoding:
                if tileGeocodeChild.tag == 'Geoposition':
                    if tileGeocodeChild.attrib['resolution'] == '10':
                        tlXRes10 = arcsiUtils.str2Float(tileGeocodeChild.find('ULX').text.strip())
                        tlYRes10 = arcsiUtils.str2Float(tileGeocodeChild.find('ULY').text.strip())
                    elif  tileGeocodeChild.attrib['resolution'] == '20':
                        tlXRes20 = arcsiUtils.str2Float(tileGeocodeChild.find('ULX').text.strip())
                        tlYRes20 = arcsiUtils.str2Float(tileGeocodeChild.find('ULY').text.strip())
                    elif  tileGeocodeChild.attrib['resolution'] == '60':
                        tlXRes60 = arcsiUtils.str2Float(tileGeocodeChild.find('ULX').text.strip())
                        tlYRes60 = arcsiUtils.str2Float(tileGeocodeChild.find('ULY').text.strip())
                    else:
                        raise ARCSIException("Resolution must be 10, 20 or 60m")
                elif tileGeocodeChild.tag == 'Size':
                    if tileGeocodeChild.attrib['resolution'] == '10':
                        xSizeRes10 = arcsiUtils.str2Int(tileGeocodeChild.find('NCOLS').text.strip())
                        ySizeRes10 = arcsiUtils.str2Int(tileGeocodeChild.find('NROWS').text.strip())
                    elif  tileGeocodeChild.attrib['resolution'] == '20':
                        xSizeRes20 = arcsiUtils.str2Int(tileGeocodeChild.find('NCOLS').text.strip())
                        ySizeRes20 = arcsiUtils.str2Int(tileGeocodeChild.find('NROWS').text.strip())
                    elif  tileGeocodeChild.attrib['resolution'] == '60':
                        xSizeRes60 = arcsiUtils.str2Int(tileGeocodeChild.find('NCOLS').text.strip())
                        ySizeRes60 = arcsiUtils.str2Int(tileGeocodeChild.find('NROWS').text.strip())
                    else:
                        raise ARCSIException("Resolution must be 10, 20 or 60m")

            brXRes10 = tlXRes10 + (xSizeRes10 * 10)
            brYRes10 = tlYRes10 - (ySizeRes10 * 10)
            brXRes20 = tlXRes20 + (xSizeRes20 * 20)
            brYRes20 = tlYRes20 - (ySizeRes20 * 20)
            brXRes60 = tlXRes60 + (xSizeRes60 * 60)
            brYRes60 = tlYRes60 - (ySizeRes60 * 60)

            self.xTL = arcsiUtils.getMinVal([tlXRes10, tlXRes20, tlXRes60])
            self.yTL = arcsiUtils.getMaxVal([tlYRes10, tlYRes20, tlYRes60])
            self.xBR = arcsiUtils.getMaxVal([brXRes10, brXRes20, brXRes60])
            self.yBR = arcsiUtils.getMinVal([brYRes10, brYRes20, brYRes60])
            self.xTR = self.xBR
            self.yTR = self.yTL
            self.xBL = self.xTL
            self.yBL = self.yBR

            self.xCentre = self.xTL + ((self.xTR - self.xTL)/2)
            self.yCentre = self.yBR + ((self.yTL - self.yBR)/2)

            self.latCentre, self.lonCentre = arcsiUtils.getLatLong(inProj, self.xCentre, self.yCentre)
            self.latTL, self.lonTL = arcsiUtils.getLatLong(inProj, self.xTL, self.yTL)
            self.latTR, self.lonTR = arcsiUtils.getLatLong(inProj, self.xTR, self.yTR)
            self.latBL, self.lonBL = arcsiUtils.getLatLong(inProj, self.xBL, self.yBL)
            self.latBR, self.lonBR = arcsiUtils.getLatLong(inProj, self.xBR, self.yBR)            

            # Get Tile angles tag.
            tileAngles = geometricInfoTag.find('Tile_Angles')

            sunAngleTag = tileAngles.find('Mean_Sun_Angle')
            self.solarZenith = arcsiUtils.str2Float(sunAngleTag.find('ZENITH_ANGLE').text.strip())
            self.solarAzimuth = arcsiUtils.str2Float(sunAngleTag.find('AZIMUTH_ANGLE').text.strip())

            senZenVals = []
            senAzVals = []
            for meanViewIncAngleTag in tileAngles.find('Mean_Viewing_Incidence_Angle_List'):
                if meanViewIncAngleTag.tag == 'Mean_Viewing_Incidence_Angle':
                    senZenVals.append(arcsiUtils.str2Float(meanViewIncAngleTag.find('ZENITH_ANGLE').text.strip()))
                    senAzVals.append(arcsiUtils.str2Float(meanViewIncAngleTag.find('AZIMUTH_ANGLE').text.strip()))

            self.sensorZenith = arcsiUtils.getMeanVal(senZenVals)
            self.sensorAzimuth = arcsiUtils.getMeanVal(senAzVals)

            ################################
            #### CODE FROM PYTHON FMASK ####
            ################################
            # Upper-left corners of images at different resolutions. As far as I can
            # work out, these coords appear to be the upper left corner of the upper left
            # pixel, i.e. equivalent to GDAL's convention. This also means that they
            # are the same for the different resolutions, which is nice. 
            self.ulxyByRes = {}
            posNodeList = tileGeocoding.findall('Geoposition')
            for posNode in posNodeList:
                res = posNode.attrib['resolution']
                ulx = float(posNode.find('ULX').text)
                uly = float(posNode.find('ULY').text)
                self.ulxyByRes[res] = (ulx, uly)

            # Sun and satellite angles. 
            sunZenithNode = tileAngles.find('Sun_Angles_Grid').find('Zenith')
            self.angleGridXres = float(sunZenithNode.find('COL_STEP').text)
            self.angleGridYres = float(sunZenithNode.find('ROW_STEP').text)
            self.sunZenithGrid = self.makeValueArray(sunZenithNode.find('Values_List'))
            sunAzimuthNode = tileAngles.find('Sun_Angles_Grid').find('Azimuth')
            self.sunAzimuthGrid = self.makeValueArray(sunAzimuthNode.find('Values_List'))
            self.anglesGridShape = self.sunAzimuthGrid.shape
            
            # Now build up the viewing angle per grid cell, from the separate layers
            # given for each detector for each band. Initially I am going to keep
            # the bands separate, just to see how that looks. 
            # The names of things in the XML suggest that these are view angles,
            # but the numbers suggest that they are angles as seen from the pixel's 
            # frame of reference on the ground, i.e. they are in fact what we ultimately want. 
            viewingAngleNodeList = tileAngles.findall('Viewing_Incidence_Angles_Grids')
            self.viewZenithDict = self.buildViewAngleArr(viewingAngleNodeList, 'Zenith')
            self.viewAzimuthDict = self.buildViewAngleArr(viewingAngleNodeList, 'Azimuth')
            
            # Make a guess at the coordinates of the angle grids. These are not given 
            # explicitly in the XML, and don't line up exactly with the other grids, so I am 
            # making a rough estimate. Because the angles don't change rapidly across these 
            # distances, it is not important if I am a bit wrong (although it would be nice
            # to be exactly correct!). 
            (ulx, uly) = self.ulxyByRes["10"]
            self.anglesULXY = (ulx - self.angleGridXres / 2.0, uly + self.angleGridYres / 2.0)
            ################################

            # Get 6S spectral response functions. 
            for specRspBand in self.specBandInfo:
                respFuncSize = len(self.specBandInfo[specRspBand].respFunc)
                tmpWvLenMin = self.specBandInfo[specRspBand].wvLenMin
                tmpWvLenMax = self.specBandInfo[specRspBand].wvLenMax
                numWvLens = math.ceil(tmpWvLenMax - tmpWvLenMin)
                lenDiff = respFuncSize - numWvLens
                if lenDiff > 0:
                    tmpWvLenMax = tmpWvLenMax + lenDiff
                elif lenDiff < 0:
                    tmpWvLenMin = tmpWvLenMin + (lenDiff*(-1))
                wvLensIn = numpy.arange(tmpWvLenMin, tmpWvLenMax, 1)
                olWVLens, olSpecResp = arcsiUtils.resampleSpectralResponseFunc(wvLensIn, self.specBandInfo[specRspBand].respFunc, 2.5, 'linear')
                self.specBandInfo[specRspBand].respFunc6S = olSpecResp
                self.specBandInfo[specRspBand].wvLenMin6S = olWVLens[0]/1000
                self.specBandInfo[specRspBand].wvLenMax6S = olWVLens[-1]/1000

        except Exception as e:
            raise e

    def getSolarIrrStdSolarGeom(self):
        """
        Get Solar Azimuth and Zenith as standard geometry.
        Azimuth: N=0, E=90, S=180, W=270.
        """
        return (self.solarAzimuth, self.solarZenith)

    def getSensorViewGeom(self):
        """
        Get sensor viewing angles
        returns (viewAzimuth, viewZenith)
        """
        return (0.0, 0.0)

    def generateOutputBaseName(self):
        """
        Customises the generic name for the WorldView2 sensor
        """
        outname = self.defaultGenBaseOutFileName()
        outname = outname + '_' + self.uniqueTileID
        return outname

    def expectedImageDataPresent(self):
        imageDataPresent = True

        if not os.path.exists(self.sen2ImgB01):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB02):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB03):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB04):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB05):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB06):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB07):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB8A):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB08):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB09):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB10):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB11):
            imageDataPresent = False
        elif not os.path.exists(self.sen2ImgB12):
            imageDataPresent = False

        return imageDataPresent

    def imgNeedMosaicking(self):
        return False

    def inImgsDiffRes(self):
        return True

    def applyImageDataMask(self, inputHeader, outputPath, outputMaskName, outputImgName, outFormat, outWKTFile):
        raise ARCSIException("Sentinel-2 does not provide any image masks, do not use the MASK option.")

    def mosaicImageTiles(self, outputPath):
        raise ARCSIException("Image data does not need mosaicking")

    def resampleImgRes(self, outputPath, resampleToLowResImg, resampleMethod='cubic'):
        outBaseName = self.generateOutputBaseName()
        if resampleToLowResImg:
            # Resample to 20 m
            self.sen2ImgB02_20m = os.path.join(outputPath, outBaseName+'_B02_20m.kea')
            rsgislib.imageutils.resampleImage2Match(self.sen2ImgB05, self.sen2ImgB02, self.sen2ImgB02_20m, 'KEA', resampleMethod, rsgislib.TYPE_16UINT)
            self.sen2ImgB03_20m = os.path.join(outputPath, outBaseName+'_B03_20m.kea')
            rsgislib.imageutils.resampleImage2Match(self.sen2ImgB05, self.sen2ImgB03, self.sen2ImgB03_20m, 'KEA', resampleMethod, rsgislib.TYPE_16UINT)
            self.sen2ImgB04_20m = os.path.join(outputPath, outBaseName+'_B04_20m.kea')
            rsgislib.imageutils.resampleImage2Match(self.sen2ImgB05, self.sen2ImgB04, self.sen2ImgB04_20m, 'KEA', resampleMethod, rsgislib.TYPE_16UINT)
            self.sen2ImgB08_20m = os.path.join(outputPath, outBaseName+'_B08_20m.kea')
            rsgislib.imageutils.resampleImage2Match(self.sen2ImgB05, self.sen2ImgB08, self.sen2ImgB08_20m, 'KEA', resampleMethod, rsgislib.TYPE_16UINT)
            self.resampleTo20m = True
        else:
            # Resample to 10 m
            self.sen2ImgB05_10m = os.path.join(outputPath, outBaseName+'_B05_10m.kea')
            rsgislib.imageutils.resampleImage2Match(self.sen2ImgB02, self.sen2ImgB05, self.sen2ImgB05_10m, 'KEA', resampleMethod, rsgislib.TYPE_16UINT)
            self.sen2ImgB06_10m = os.path.join(outputPath, outBaseName+'_B06_10m.kea')
            rsgislib.imageutils.resampleImage2Match(self.sen2ImgB02, self.sen2ImgB06, self.sen2ImgB06_10m, 'KEA', resampleMethod, rsgislib.TYPE_16UINT)
            self.sen2ImgB07_10m = os.path.join(outputPath, outBaseName+'_B07_10m.kea')
            rsgislib.imageutils.resampleImage2Match(self.sen2ImgB02, self.sen2ImgB07, self.sen2ImgB07_10m, 'KEA', resampleMethod, rsgislib.TYPE_16UINT)
            self.sen2ImgB8A_10m = os.path.join(outputPath, outBaseName+'_B08A_10m.kea')
            rsgislib.imageutils.resampleImage2Match(self.sen2ImgB02, self.sen2ImgB8A, self.sen2ImgB8A_10m, 'KEA', resampleMethod, rsgislib.TYPE_16UINT)
            self.sen2ImgB11_10m = os.path.join(outputPath, outBaseName+'_B11_10m.kea')
            rsgislib.imageutils.resampleImage2Match(self.sen2ImgB02, self.sen2ImgB11, self.sen2ImgB11_10m, 'KEA', resampleMethod, rsgislib.TYPE_16UINT)
            self.sen2ImgB12_10m = os.path.join(outputPath, outBaseName+'_B12_10m.kea')
            rsgislib.imageutils.resampleImage2Match(self.sen2ImgB02, self.sen2ImgB12, self.sen2ImgB12_10m, 'KEA', resampleMethod, rsgislib.TYPE_16UINT)
            self.resampleTo20m = False

    def generateValidImageDataMask(self, outputPath, outputMaskName, viewAngleImg, outFormat):
        print("Generate valid image mask")
        # Generate the valid image mask
        outputImage = os.path.join(outputPath, outputMaskName)
        inImgBands = []
        if self.resampleTo20m:
            inImgBands.append(self.sen2ImgB02_20m)
            inImgBands.append(self.sen2ImgB03_20m)
            inImgBands.append(self.sen2ImgB04_20m)
            inImgBands.append(self.sen2ImgB05)
            inImgBands.append(self.sen2ImgB06)
            inImgBands.append(self.sen2ImgB07)
            inImgBands.append(self.sen2ImgB08_20m)
            inImgBands.append(self.sen2ImgB8A)
            inImgBands.append(self.sen2ImgB11)
            inImgBands.append(self.sen2ImgB12)
        else: 
            inImgBands.append(self.sen2ImgB02)
            inImgBands.append(self.sen2ImgB03)
            inImgBands.append(self.sen2ImgB04)
            inImgBands.append(self.sen2ImgB05_10m)
            inImgBands.append(self.sen2ImgB06_10m)
            inImgBands.append(self.sen2ImgB07_10m)
            inImgBands.append(self.sen2ImgB08)
            inImgBands.append(self.sen2ImgB8A_10m)
            inImgBands.append(self.sen2ImgB11_10m)
            inImgBands.append(self.sen2ImgB12_10m)
        rsgislib.imageutils.genValidMask(inimages=inImgBands, outimage=outputImage, gdalformat=outFormat, nodata=self.inNoDataVal)

        ################################
        #### CODE FROM PYTHON FMASK ####
        ################################
        tmpViewAngleImg = os.path.join(outputPath, os.path.splitext(os.path.basename(viewAngleImg))[0]+'tmp.kea')
        # Generate the acquasition anagles image.
        # This scale value will convert between DN and radians in output image file, radians = dn * SCALE_TO_RADIANS
        SCALE_TO_RADIANS = 0.01

        drvr = gdal.GetDriverByName('KEA')
        (nrows, ncols) = self.anglesGridShape
        ds = drvr.Create(tmpViewAngleImg, ncols, nrows, 4, gdal.GDT_Int16)
        gt = (self.anglesULXY[0], self.angleGridXres, 0, self.anglesULXY[1], 0.0, -self.angleGridYres)
        ds.SetGeoTransform(gt)
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(int(self.epsgCode))
        ds.SetProjection(sr.ExportToWkt())

        nullValDN = 1000

        # Get a sorted list of the Sentinel-2 band names. Note that sometimes this
        # is an incomplete list of band names, which appears to be due to a bug in 
        # earlier versions of ESA's processing software. I suspect it relates to 
        # Anomaly number 11 in the following page. 
        # https://sentinel.esa.int/web/sentinel/news/-/article/new-processing-baseline-for-sentinel-2-products
        bandNames = sorted(self.viewAzimuthDict.keys())

        # Mean over all bands
        satAzDeg = numpy.array([self.viewAzimuthDict[i] for i in bandNames])
        satAzDegMeanOverBands = satAzDeg.mean(axis=0)

        satZenDeg = numpy.array([self.viewZenithDict[i] for i in bandNames])
        satZenDegMeanOverBands = satZenDeg.mean(axis=0)

        sunAzDeg = self.sunAzimuthGrid
        sunZenDeg = self.sunZenithGrid

        stackDeg = numpy.array([satAzDegMeanOverBands, satZenDegMeanOverBands, sunAzDeg, sunZenDeg])
        stackRadians = numpy.radians(stackDeg)

        stackDN = numpy.round(stackRadians / SCALE_TO_RADIANS).astype(numpy.int16)
        nullmask = numpy.isnan(stackDeg)
        stackDN[nullmask] = nullValDN

        lnames = ['SatelliteAzimuth', 'SatelliteZenith', 'SunAzimuth', 'SunZenith']
        for i in range(ds.RasterCount):
            b = ds.GetRasterBand(i+1)
            b.WriteArray(stackDN[i])
            b.SetNoDataValue(nullValDN)
            b.SetDescription(lnames[i])
        del ds
        ################################
        rsgislib.imageutils.popImageStats(tmpViewAngleImg, usenodataval=True, nodataval=1000, calcpyramids=False)
        rsgislib.imageutils.resampleImage2Match(outputImage, tmpViewAngleImg, viewAngleImg, outFormat, 'nearestneighbour', datatype=None)
        dataset = gdal.Open(viewAngleImg, gdal.GA_Update)
        if not dataset is None:
            dataset.GetRasterBand(1).SetDescription("SatelliteAzimuth")
            dataset.GetRasterBand(2).SetDescription("SatelliteZenith")
            dataset.GetRasterBand(3).SetDescription("SunAzimuth")
            dataset.GetRasterBand(4).SetDescription("SunZenith")
        dataset = None
        drvr.Delete(tmpViewAngleImg)

        self.viewAnglesFmaskImg = viewAngleImg

        return outputImage

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        outputImage = os.path.join(outputPath, outputName)

        s2Band = collections.namedtuple('S2Band', ['bandName', 'fileName', 'bandIndex', 'satVal'])
        bandDefnSeq = list()
        self.inSatDataVal = float(self.inSatDataVal)
        if self.resampleTo20m:
            bandDefnSeq.append(s2Band(bandName="Blue", fileName=self.sen2ImgB02_20m, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="Green", fileName=self.sen2ImgB03_20m, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="Red", fileName=self.sen2ImgB04_20m, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="RE_B05", fileName=self.sen2ImgB05, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="RE_B06", fileName=self.sen2ImgB06, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="RE_B07", fileName=self.sen2ImgB07, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="NIR_B08", fileName=self.sen2ImgB08_20m, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="NIR_B08A", fileName=self.sen2ImgB8A, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="SWIR1", fileName=self.sen2ImgB11, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="SWIR2", fileName=self.sen2ImgB12, bandIndex=1, satVal=self.inSatDataVal))
        else:
            bandDefnSeq.append(s2Band(bandName="Blue", fileName=self.sen2ImgB02, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="Green", fileName=self.sen2ImgB03, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="Red", fileName=self.sen2ImgB04, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="RE_B05", fileName=self.sen2ImgB05_10m, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="RE_B06", fileName=self.sen2ImgB06_10m, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="RE_B07", fileName=self.sen2ImgB07_10m, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="NIR_B08", fileName=self.sen2ImgB08, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="NIR_B08A", fileName=self.sen2ImgB8A_10m, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="SWIR1", fileName=self.sen2ImgB11_10m, bandIndex=1, satVal=self.inSatDataVal))
            bandDefnSeq.append(s2Band(bandName="SWIR2", fileName=self.sen2ImgB12_10m, bandIndex=1, satVal=self.inSatDataVal))

        rsgislib.imagecalibration.saturatedPixelsMask(outputImage, outFormat, bandDefnSeq)
        return outputImage

    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat):
        print("Converting to Radiance")
        outputReflImage = os.path.join(outputPath, outputReflName)
        outputThermalImage = None

        inImgBands = list()
        if self.resampleTo20m:
            inImgBands.append(self.sen2ImgB02_20m)
            inImgBands.append(self.sen2ImgB03_20m)
            inImgBands.append(self.sen2ImgB04_20m)
            inImgBands.append(self.sen2ImgB05)
            inImgBands.append(self.sen2ImgB06)
            inImgBands.append(self.sen2ImgB07)
            inImgBands.append(self.sen2ImgB08_20m)
            inImgBands.append(self.sen2ImgB8A)
            inImgBands.append(self.sen2ImgB11)
            inImgBands.append(self.sen2ImgB12)
        else: 
            inImgBands.append(self.sen2ImgB02)
            inImgBands.append(self.sen2ImgB03)
            inImgBands.append(self.sen2ImgB04)
            inImgBands.append(self.sen2ImgB05_10m)
            inImgBands.append(self.sen2ImgB06_10m)
            inImgBands.append(self.sen2ImgB07_10m)
            inImgBands.append(self.sen2ImgB08)
            inImgBands.append(self.sen2ImgB8A_10m)
            inImgBands.append(self.sen2ImgB11_10m)
            inImgBands.append(self.sen2ImgB12_10m)

        solarIrrVals = list()
        IrrVal = collections.namedtuple('SolarIrradiance', ['irradiance'])
        solarIrrVals.append(IrrVal(irradiance=self.esun_B1))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B2))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B3))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B4))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B5))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B6))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B7))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B8))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B11))
        solarIrrVals.append(IrrVal(irradiance=self.esun_B12))

        rsgislib.imagecalibration.toaRefl2Radiance(inImgBands, outputReflImage, outFormat, rsgislib.TYPE_32FLOAT, self.quantificationVal, self.earthSunDist_U, self.solarZenith, solarIrrVals)
        return outputReflImage, outputThermalImage

    def convertThermalToBrightness(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        raise ARCSIException("Not Implemented")

    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        print("Converting to TOA")
        outputImage = os.path.join(outputPath, outputName)
        inImgBands = list()
        if self.resampleTo20m:
            inImgBands.append(self.sen2ImgB02_20m)
            inImgBands.append(self.sen2ImgB03_20m)
            inImgBands.append(self.sen2ImgB04_20m)
            inImgBands.append(self.sen2ImgB05)
            inImgBands.append(self.sen2ImgB06)
            inImgBands.append(self.sen2ImgB07)
            inImgBands.append(self.sen2ImgB08_20m)
            inImgBands.append(self.sen2ImgB8A)
            inImgBands.append(self.sen2ImgB11)
            inImgBands.append(self.sen2ImgB12)
        else: 
            inImgBands.append(self.sen2ImgB02)
            inImgBands.append(self.sen2ImgB03)
            inImgBands.append(self.sen2ImgB04)
            inImgBands.append(self.sen2ImgB05_10m)
            inImgBands.append(self.sen2ImgB06_10m)
            inImgBands.append(self.sen2ImgB07_10m)
            inImgBands.append(self.sen2ImgB08)
            inImgBands.append(self.sen2ImgB8A_10m)
            inImgBands.append(self.sen2ImgB11_10m)
            inImgBands.append(self.sen2ImgB12_10m)

        rsgislib.imagecalc.calcImageRescale(inImgBands, outputImage, outFormat, rsgislib.TYPE_16UINT, self.inNoDataVal, 0, self.quantificationVal, 0, 0, scaleFactor)
        self.imgIntScaleFactor = scaleFactor
        return outputImage

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor):
        outputImage = os.path.join(outputPath, outputName)
        tmpBaseName = os.path.splitext(outputName)[0]
        tmpBaseDIR = os.path.join(tmpPath, tmpBaseName)

        tmpDIRExisted = True
        if not os.path.exists(tmpBaseDIR):
            os.makedirs(tmpBaseDIR)
            tmpDIRExisted = False

        # Using python-fmask (http://pythonfmask.org)
        from fmask import config
        from fmask import fmask
        from rios import fileinfo

        sen2ImgB01_tmp = os.path.join(tmpBaseDIR, tmpBaseName+'_B01.kea')
        rsgislib.imageutils.resampleImage2Match(inputReflImage, self.sen2ImgB01, sen2ImgB01_tmp, 'KEA', 'nearestneighbour', rsgislib.TYPE_16UINT)
        sen2ImgB09_tmp = os.path.join(tmpBaseDIR, tmpBaseName+'_B09.kea')
        rsgislib.imageutils.resampleImage2Match(inputReflImage, self.sen2ImgB09, sen2ImgB09_tmp, 'KEA', 'nearestneighbour', rsgislib.TYPE_16UINT)
        sen2ImgB10_tmp = os.path.join(tmpBaseDIR, tmpBaseName+'_B10.kea')
        rsgislib.imageutils.resampleImage2Match(inputReflImage, self.sen2ImgB10, sen2ImgB10_tmp, 'KEA', 'nearestneighbour', rsgislib.TYPE_16UINT)

        tmpTOAImg = os.path.join(tmpBaseDIR, tmpBaseName+'_pyfmasktmpTOA.kea')
        if self.imgIntScaleFactor != 10000:
            vrtImgB1B9B10 = os.path.join(tmpBaseDIR, tmpBaseName+'_b01b09b10_60mBands.vrt')
            cmdVRT = 'gdalbuildvrt -separate ' + vrtImgB1B9B10 + ' ' + sen2ImgB01_tmp + ' ' + sen2ImgB09_tmp + ' ' + sen2ImgB10_tmp
            try:
                subprocess.call(cmdVRT, shell=True)
            except OSError as e:
                raise ARCSIException('Could not successful execute command: ' + cmdVRT)
            sen2ImgB1B9B10Rescaled = os.path.join(tmpBaseDIR, tmpBaseName+'_B01B09B10Rescaled.kea')
            rsgislib.imagecalc.imageMath(vrtImgB1B9B10, sen2ImgB1B9B10Rescaled, '(b1/10000)*'+str(self.imgIntScaleFactor), 'KEA', rsgislib.TYPE_16UINT)
            rsgislib.imageutils.stackImageBands([inputReflImage, sen2ImgB1B9B10Rescaled], None, tmpTOAImg, None, 0, 'KEA', rsgislib.TYPE_16UINT)
        else:
            rsgislib.imageutils.stackImageBands([inputReflImage, sen2ImgB01_tmp, sen2ImgB09_tmp, sen2ImgB10_tmp], None, tmpTOAImg, None, 0, 'KEA', rsgislib.TYPE_16UINT)

        fmaskReflImg = os.path.join(tmpBaseDIR, tmpBaseName+'_pyfmaskRefl.kea')
        rsgislib.imageutils.selectImageBands(tmpTOAImg, fmaskReflImg, 'KEA', rsgislib.TYPE_16UINT, [11,1,2,3,4,5,6,7,8,12,13,9,10])

        anglesInfo = config.AnglesFileInfo(self.viewAnglesFmaskImg, 3, self.viewAnglesFmaskImg, 2, self.viewAnglesFmaskImg, 1, self.viewAnglesFmaskImg, 0)        

        fmaskCloudsImg = os.path.join(tmpBaseDIR, tmpBaseName+'_pyfmaskCloudsResult.kea')
        fmaskFilenames = config.FmaskFilenames()
        fmaskFilenames.setTOAReflectanceFile(fmaskReflImg)
        fmaskFilenames.setSaturationMask(inputSatImage)
        fmaskFilenames.setOutputCloudMaskFile(fmaskCloudsImg)
        
        fmaskConfig = config.FmaskConfig(config.FMASK_SENTINEL2)
        fmaskConfig.setAnglesInfo(anglesInfo)
        fmaskConfig.setKeepIntermediates(True)
        fmaskConfig.setVerbose(True)
        fmaskConfig.setTempDir(tmpBaseDIR)
        fmaskConfig.setTOARefScaling(float(self.imgIntScaleFactor))
        fmaskConfig.setMinCloudSize(8)
        
        # Work out a suitable buffer size, in pixels, dependent on the resolution of the input TOA image
        #toaImgInfo = fileinfo.ImageInfo(fmaskReflImg)
        fmaskConfig.setCloudBufferSize(10)
        fmaskConfig.setShadowBufferSize(10)
        
        fmask.doFmask(fmaskFilenames, fmaskConfig)

        rsgislib.imagecalc.imageMath(fmaskCloudsImg, outputImage, '(b1==2)?1:(b1==3)?2:0', outFormat, rsgislib.TYPE_16UINT)
        if outFormat == 'KEA':
            rsgislib.rastergis.populateStats(outputImage, True, True)
            ratDataset = gdal.Open(outputImage, gdal.GA_Update)
            red = rat.readColumn(ratDataset, 'Red')
            green = rat.readColumn(ratDataset, 'Green')
            blue = rat.readColumn(ratDataset, 'Blue')
            ClassName = numpy.empty_like(red, dtype=numpy.dtype('a255'))

            if(red.shape[0] > 0):
                red[1] = 0
                green[1] = 0
                blue[1] = 255
                ClassName[1] = 'Clouds'

                if(red.shape[0] > 1):
                    red[2] = 0
                    green[2] = 255
                    blue[2] = 255
                    ClassName[2] = 'Shadows'

            rat.writeColumn(ratDataset, "Red", red)
            rat.writeColumn(ratDataset, "Green", green)
            rat.writeColumn(ratDataset, "Blue", blue)
            rat.writeColumn(ratDataset, "ClassName", ClassName)
            ratDataset = None
        rsgislib.imageutils.copyProjFromImage(outputImage, inputReflImage)

        if not self.debugMode:
            if not tmpDIRExisted:
                shutil.rmtree(tmpBaseDIR, ignore_errors=True)

        return outputImage

    def createCloudMaskDataArray(self, inImgDataArr):
        return inImgDataArr

    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        sixsCoeffs = numpy.zeros((10, 6), dtype=numpy.float32)
        # Set up 6S model
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.User()
        s.geometry.solar_z = self.solarZenith
        s.geometry.solar_a = self.solarAzimuth
        s.geometry.view_z = self.sensorZenith
        s.geometry.view_a = self.sensorAzimuth
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute)/60.0
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        if useBRDF:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrBRDFFromRadiance(200)
        else:
            s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal

        # Blue
        s.wavelength = Py6S.Wavelength(self.specBandInfo['B2'].wvLenMin6S, self.specBandInfo['B2'].wvLenMax6S, self.specBandInfo['B2'].respFunc6S)
        s.run()
        sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

        # Green
        s.wavelength = Py6S.Wavelength(self.specBandInfo['B3'].wvLenMin6S, self.specBandInfo['B3'].wvLenMax6S, self.specBandInfo['B3'].respFunc6S)
        s.run()
        sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

        # Red
        s.wavelength = Py6S.Wavelength(self.specBandInfo['B4'].wvLenMin6S, self.specBandInfo['B4'].wvLenMax6S, self.specBandInfo['B4'].respFunc6S)
        s.run()
        sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

        # RE B5
        s.wavelength = Py6S.Wavelength(self.specBandInfo['B5'].wvLenMin6S, self.specBandInfo['B5'].wvLenMax6S, self.specBandInfo['B5'].respFunc6S)
        s.run()
        sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[3,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[3,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[3,5] = float(s.outputs.values['environmental_irradiance'])

        # RE B6
        s.wavelength = Py6S.Wavelength(self.specBandInfo['B6'].wvLenMin6S, self.specBandInfo['B6'].wvLenMax6S, self.specBandInfo['B6'].respFunc6S)
        s.run()
        sixsCoeffs[4,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[4,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[4,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[4,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[4,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[4,5] = float(s.outputs.values['environmental_irradiance'])

        # RE B7
        s.wavelength = Py6S.Wavelength(self.specBandInfo['B7'].wvLenMin6S, self.specBandInfo['B7'].wvLenMax6S, self.specBandInfo['B7'].respFunc6S)
        s.run()
        sixsCoeffs[5,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[5,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[5,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[5,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[5,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[5,5] = float(s.outputs.values['environmental_irradiance'])

        # NIR B8
        s.wavelength = Py6S.Wavelength(self.specBandInfo['B8'].wvLenMin6S, self.specBandInfo['B8'].wvLenMax6S, self.specBandInfo['B8'].respFunc6S)
        s.run()
        sixsCoeffs[6,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[6,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[6,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[6,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[6,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[6,5] = float(s.outputs.values['environmental_irradiance'])

        # NIR B8A
        s.wavelength = Py6S.Wavelength(self.specBandInfo['B8A'].wvLenMin6S, self.specBandInfo['B8A'].wvLenMax6S, self.specBandInfo['B8A'].respFunc6S)
        s.run()
        sixsCoeffs[7,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[7,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[7,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[7,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[7,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[7,5] = float(s.outputs.values['environmental_irradiance'])

        # SWIR 1
        s.wavelength = Py6S.Wavelength(self.specBandInfo['B11'].wvLenMin6S, self.specBandInfo['B11'].wvLenMax6S, self.specBandInfo['B11'].respFunc6S)
        s.run()
        sixsCoeffs[8,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[8,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[8,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[8,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[8,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[8,5] = float(s.outputs.values['environmental_irradiance'])

        # SWIR 2
        s.wavelength = Py6S.Wavelength(self.specBandInfo['B12'].wvLenMin6S, self.specBandInfo['B12'].wvLenMax6S, self.specBandInfo['B12'].respFunc6S)
        s.run()
        sixsCoeffs[9,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[9,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[9,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[9,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[9,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[9,5] = float(s.outputs.values['environmental_irradiance'])

        return sixsCoeffs

    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX', 'DirIrr', 'DifIrr', 'EnvIrr'])
        imgBandCoeffs = list()

        sixsCoeffs = self.calc6SCoefficients(aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF)

        imgBandCoeffs.append(Band6S(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
        imgBandCoeffs.append(Band6S(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
        imgBandCoeffs.append(Band6S(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
        imgBandCoeffs.append(Band6S(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
        imgBandCoeffs.append(Band6S(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2]), DirIrr=float(sixsCoeffs[4,3]), DifIrr=float(sixsCoeffs[4,4]), EnvIrr=float(sixsCoeffs[4,5])))
        imgBandCoeffs.append(Band6S(band=6, aX=float(sixsCoeffs[5,0]), bX=float(sixsCoeffs[5,1]), cX=float(sixsCoeffs[5,2]), DirIrr=float(sixsCoeffs[5,3]), DifIrr=float(sixsCoeffs[5,4]), EnvIrr=float(sixsCoeffs[5,5])))
        imgBandCoeffs.append(Band6S(band=7, aX=float(sixsCoeffs[6,0]), bX=float(sixsCoeffs[6,1]), cX=float(sixsCoeffs[6,2]), DirIrr=float(sixsCoeffs[6,3]), DifIrr=float(sixsCoeffs[6,4]), EnvIrr=float(sixsCoeffs[6,5])))
        imgBandCoeffs.append(Band6S(band=8, aX=float(sixsCoeffs[7,0]), bX=float(sixsCoeffs[7,1]), cX=float(sixsCoeffs[7,2]), DirIrr=float(sixsCoeffs[7,3]), DifIrr=float(sixsCoeffs[7,4]), EnvIrr=float(sixsCoeffs[7,5])))
        imgBandCoeffs.append(Band6S(band=9, aX=float(sixsCoeffs[8,0]), bX=float(sixsCoeffs[8,1]), cX=float(sixsCoeffs[8,2]), DirIrr=float(sixsCoeffs[8,3]), DifIrr=float(sixsCoeffs[8,4]), EnvIrr=float(sixsCoeffs[8,5])))
        imgBandCoeffs.append(Band6S(band=10, aX=float(sixsCoeffs[9,0]), bX=float(sixsCoeffs[9,1]), cX=float(sixsCoeffs[9,2]), DirIrr=float(sixsCoeffs[9,3]), DifIrr=float(sixsCoeffs[9,4]), EnvIrr=float(sixsCoeffs[9,5])))

        for band in imgBandCoeffs:
            print(band)
        rsgislib.imagecalibration.apply6SCoeffSingleParam(inputRadImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, imgBandCoeffs)
        return outputImage

    def convertImageToSurfaceReflDEMElevLUT(self, inputRadImage, inputDEMFile, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, scaleFactor, elevCoeffs=None):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevCoeffs is None:
            print("Build an LUT for elevation values.")
            elev6SCoeffsLUT = self.buildElevation6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax)
            print("LUT has been built.")

            elevLUTFeat = collections.namedtuple('ElevLUTFeat', ['Elev', 'Coeffs'])
            Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX', 'DirIrr', 'DifIrr', 'EnvIrr'])

            elevCoeffs = list()
            for elevLUT in elev6SCoeffsLUT:
                imgBandCoeffs = list()
                sixsCoeffs = elevLUT.Coeffs
                elevVal = elevLUT.Elev
                imgBandCoeffs.append(Band6S(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
                imgBandCoeffs.append(Band6S(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
                imgBandCoeffs.append(Band6S(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
                imgBandCoeffs.append(Band6S(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
                imgBandCoeffs.append(Band6S(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2]), DirIrr=float(sixsCoeffs[4,3]), DifIrr=float(sixsCoeffs[4,4]), EnvIrr=float(sixsCoeffs[4,5])))
                imgBandCoeffs.append(Band6S(band=6, aX=float(sixsCoeffs[5,0]), bX=float(sixsCoeffs[5,1]), cX=float(sixsCoeffs[5,2]), DirIrr=float(sixsCoeffs[5,3]), DifIrr=float(sixsCoeffs[5,4]), EnvIrr=float(sixsCoeffs[5,5])))
                imgBandCoeffs.append(Band6S(band=7, aX=float(sixsCoeffs[6,0]), bX=float(sixsCoeffs[6,1]), cX=float(sixsCoeffs[6,2]), DirIrr=float(sixsCoeffs[6,3]), DifIrr=float(sixsCoeffs[6,4]), EnvIrr=float(sixsCoeffs[6,5])))
                imgBandCoeffs.append(Band6S(band=8, aX=float(sixsCoeffs[7,0]), bX=float(sixsCoeffs[7,1]), cX=float(sixsCoeffs[7,2]), DirIrr=float(sixsCoeffs[7,3]), DifIrr=float(sixsCoeffs[7,4]), EnvIrr=float(sixsCoeffs[7,5])))
                imgBandCoeffs.append(Band6S(band=9, aX=float(sixsCoeffs[8,0]), bX=float(sixsCoeffs[8,1]), cX=float(sixsCoeffs[8,2]), DirIrr=float(sixsCoeffs[8,3]), DifIrr=float(sixsCoeffs[8,4]), EnvIrr=float(sixsCoeffs[8,5])))
                imgBandCoeffs.append(Band6S(band=10, aX=float(sixsCoeffs[9,0]), bX=float(sixsCoeffs[9,1]), cX=float(sixsCoeffs[9,2]), DirIrr=float(sixsCoeffs[9,3]), DifIrr=float(sixsCoeffs[9,4]), EnvIrr=float(sixsCoeffs[9,5])))
                elevCoeffs.append(elevLUTFeat(Elev=float(elevVal), Coeffs=imgBandCoeffs))

        rsgislib.imagecalibration.apply6SCoeffElevLUTParam(inputRadImage, inputDEMFile, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, elevCoeffs)
        return outputImage, elevCoeffs

    def convertImageToSurfaceReflAOTDEMElevLUT(self, inputRadImage, inputDEMFile, inputAOTImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax, scaleFactor, elevAOTCoeffs=None):
        print("Converting to Surface Reflectance")
        outputImage = os.path.join(outputPath, outputName)

        if elevAOTCoeffs is None:
            print("Build an LUT for elevation and AOT values.")
            elevAOT6SCoeffsLUT = self.buildElevationAOT6SCoeffLUT(aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax)

            elevLUTFeat = collections.namedtuple('ElevLUTFeat', ['Elev', 'Coeffs'])
            aotLUTFeat = collections.namedtuple('AOTLUTFeat', ['AOT', 'Coeffs'])
            Band6S = collections.namedtuple('Band6SCoeff', ['band', 'aX', 'bX', 'cX', 'DirIrr', 'DifIrr', 'EnvIrr'])

            elevAOTCoeffs = list()
            for elevLUT in elevAOT6SCoeffsLUT:
                elevVal = elevLUT.Elev
                aotLUT = elevLUT.Coeffs
                aot6SCoeffsOut = list()
                for aotFeat in aotLUT:
                    sixsCoeffs = aotFeat.Coeffs
                    aotVal = aotFeat.AOT
                    imgBandCoeffs = list()
                    imgBandCoeffs.append(Band6S(band=1, aX=float(sixsCoeffs[0,0]), bX=float(sixsCoeffs[0,1]), cX=float(sixsCoeffs[0,2]), DirIrr=float(sixsCoeffs[0,3]), DifIrr=float(sixsCoeffs[0,4]), EnvIrr=float(sixsCoeffs[0,5])))
                    imgBandCoeffs.append(Band6S(band=2, aX=float(sixsCoeffs[1,0]), bX=float(sixsCoeffs[1,1]), cX=float(sixsCoeffs[1,2]), DirIrr=float(sixsCoeffs[1,3]), DifIrr=float(sixsCoeffs[1,4]), EnvIrr=float(sixsCoeffs[1,5])))
                    imgBandCoeffs.append(Band6S(band=3, aX=float(sixsCoeffs[2,0]), bX=float(sixsCoeffs[2,1]), cX=float(sixsCoeffs[2,2]), DirIrr=float(sixsCoeffs[2,3]), DifIrr=float(sixsCoeffs[2,4]), EnvIrr=float(sixsCoeffs[2,5])))
                    imgBandCoeffs.append(Band6S(band=4, aX=float(sixsCoeffs[3,0]), bX=float(sixsCoeffs[3,1]), cX=float(sixsCoeffs[3,2]), DirIrr=float(sixsCoeffs[3,3]), DifIrr=float(sixsCoeffs[3,4]), EnvIrr=float(sixsCoeffs[3,5])))
                    imgBandCoeffs.append(Band6S(band=5, aX=float(sixsCoeffs[4,0]), bX=float(sixsCoeffs[4,1]), cX=float(sixsCoeffs[4,2]), DirIrr=float(sixsCoeffs[4,3]), DifIrr=float(sixsCoeffs[4,4]), EnvIrr=float(sixsCoeffs[4,5])))
                    imgBandCoeffs.append(Band6S(band=6, aX=float(sixsCoeffs[5,0]), bX=float(sixsCoeffs[5,1]), cX=float(sixsCoeffs[5,2]), DirIrr=float(sixsCoeffs[5,3]), DifIrr=float(sixsCoeffs[5,4]), EnvIrr=float(sixsCoeffs[5,5])))
                    imgBandCoeffs.append(Band6S(band=7, aX=float(sixsCoeffs[6,0]), bX=float(sixsCoeffs[6,1]), cX=float(sixsCoeffs[6,2]), DirIrr=float(sixsCoeffs[6,3]), DifIrr=float(sixsCoeffs[6,4]), EnvIrr=float(sixsCoeffs[6,5])))
                    imgBandCoeffs.append(Band6S(band=8, aX=float(sixsCoeffs[7,0]), bX=float(sixsCoeffs[7,1]), cX=float(sixsCoeffs[7,2]), DirIrr=float(sixsCoeffs[7,3]), DifIrr=float(sixsCoeffs[7,4]), EnvIrr=float(sixsCoeffs[7,5])))
                    imgBandCoeffs.append(Band6S(band=9, aX=float(sixsCoeffs[8,0]), bX=float(sixsCoeffs[8,1]), cX=float(sixsCoeffs[8,2]), DirIrr=float(sixsCoeffs[8,3]), DifIrr=float(sixsCoeffs[8,4]), EnvIrr=float(sixsCoeffs[8,5])))
                    imgBandCoeffs.append(Band6S(band=10, aX=float(sixsCoeffs[9,0]), bX=float(sixsCoeffs[9,1]), cX=float(sixsCoeffs[9,2]), DirIrr=float(sixsCoeffs[9,3]), DifIrr=float(sixsCoeffs[9,4]), EnvIrr=float(sixsCoeffs[9,5])))
                    aot6SCoeffsOut.append(aotLUTFeat(AOT=float(aotVal), Coeffs=imgBandCoeffs))
                elevAOTCoeffs.append(elevLUTFeat(Elev=float(elevVal), Coeffs=aot6SCoeffsOut))

        rsgislib.imagecalibration.apply6SCoeffElevAOTLUTParam(inputRadImage, inputDEMFile, inputAOTImage, outputImage, outFormat, rsgislib.TYPE_16UINT, scaleFactor, 0, True, elevAOTCoeffs)

    def run6SToOptimiseAODValue(self, aotVal, radBlueVal, predBlueVal, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        print("Testing AOD Val: ", aotVal,)
        s = Py6S.SixS()
        s.atmos_profile = atmosProfile
        s.aero_profile = aeroProfile
        s.ground_reflectance = grdRefl
        s.geometry = Py6S.Geometry.User()
        s.geometry.solar_z = self.solarZenith
        s.geometry.solar_a = self.solarAzimuth
        s.geometry.view_z = self.sensorZenith
        s.geometry.view_a = self.sensorAzimuth
        s.geometry.month = self.acquisitionTime.month
        s.geometry.day = self.acquisitionTime.day
        s.geometry.gmt_decimal_hour = float(self.acquisitionTime.hour) + float(self.acquisitionTime.minute)/60.0
        s.geometry.latitude = self.latCentre
        s.geometry.longitude = self.lonCentre
        s.altitudes = Py6S.Altitudes()
        s.altitudes.set_target_custom_altitude(surfaceAltitude)
        s.altitudes.set_sensor_satellite_level()
        s.atmos_corr = Py6S.AtmosCorr.AtmosCorrLambertianFromRadiance(200)
        s.aot550 = aotVal

        # Band 1 (Blue!)
        s.wavelength = Py6S.Wavelength(self.specBandInfo['B2'].wvLenMin6S, self.specBandInfo['B2'].wvLenMax6S, self.specBandInfo['B2'].respFunc6S)
        s.run()
        aX = float(s.outputs.values['coef_xa'])
        bX = float(s.outputs.values['coef_xb'])
        cX = float(s.outputs.values['coef_xc'])

        tmpVal = (aX*radBlueVal)-bX;
        reflBlueVal = tmpVal/(1.0+cX*tmpVal)
        outDist = math.sqrt(math.pow((reflBlueVal - predBlueVal),2))
        print("\taX: ", aX, " bX: ", bX, " cX: ", cX, "     Dist = ", outDist)
        return outDist

    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        print("Not implemented\n")
        sys.exit()

    def estimateImageToAODUsingDDV(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax):
        print("Not implemented\n")
        sys.exit()

    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, simpleDOS, dosOutRefl):
        raise ARCSIException("Not Implemented")

    def estimateSingleAOTFromDOS(self, radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl):
        try:
            return self.estimateSingleAOTFromDOSBandImpl(radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl, 1)
        except Exception as e:
            raise

    def setBandNames(self, imageFile):
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        if not dataset is None:
            dataset.GetRasterBand(1).SetDescription("Blue")
            dataset.GetRasterBand(2).SetDescription("Green")
            dataset.GetRasterBand(3).SetDescription("Red")
            dataset.GetRasterBand(4).SetDescription("RE_B5")
            dataset.GetRasterBand(5).SetDescription("RE_B6")
            dataset.GetRasterBand(6).SetDescription("RE_B7")
            dataset.GetRasterBand(7).SetDescription("NIR_B8")
            dataset.GetRasterBand(8).SetDescription("NIR_B8A")
            dataset.GetRasterBand(9).SetDescription("SWIR1")
            dataset.GetRasterBand(10).SetDescription("SWIR2")
            dataset = None
        else:
            print("Could not open image to set band names: ", imageFile)

    def cleanLocalFollowProcessing(self):
        if not self.debugMode:
            gdalDriver = gdal.GetDriverByName('KEA')
            if self.resampleTo20m:
                gdalDriver.Delete(self.sen2ImgB02_20m)
                gdalDriver.Delete(self.sen2ImgB03_20m)
                gdalDriver.Delete(self.sen2ImgB04_20m)
                gdalDriver.Delete(self.sen2ImgB08_20m)
            else:
                gdalDriver.Delete(self.sen2ImgB05_10m)
                gdalDriver.Delete(self.sen2ImgB06_10m)
                gdalDriver.Delete(self.sen2ImgB07_10m)
                gdalDriver.Delete(self.sen2ImgB8A_10m)
                gdalDriver.Delete(self.sen2ImgB11_10m)
                gdalDriver.Delete(self.sen2ImgB12_10m)


