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
        self.sen2ImgB08A = ''
        self.sen2ImgB08 = ''
        self.sen2ImgB09 = ''
        self.sen2ImgB10 = ''
        self.sen2ImgB11 = ''
        self.sen2ImgB12 = ''
        self.sen2ImgTCI = ''

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

        self.reflConvert_U = 0.0
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
                    tmpFiles = glob.glob(os.path.join(self.sen2FileBaseDIR, imgFile+'*'))
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
                            self.sen2ImgB08A = tmpFiles[0]
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
            self.reflConvert_U = arcsiUtils.str2Float(reflectanceConversionTag.find('U').text.strip())
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
            epsgCode = arcsiUtils.str2Int(epsgCodeStr.split(':')[1])
            inProj = osr.SpatialReference()
            inProj.ImportFromEPSG(epsgCode)
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
        elif not os.path.exists(self.sen2ImgB08A):
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

    def resampleImgRes(self, outputPath, resampleToLowResImg):
        raise ARCSIException("Image resampling has not yet been implemented...")

    def generateValidImageDataMask(self, outputPath, outputMaskName, viewAngleImg, outFormat):
        raise ARCSIException("Don't know how to create a valid data mask for Sentinel-2")

    def convertImageToRadiance(self, outputPath, outputReflName, outputThermalName, outFormat):
        print("Converting to Radiance")
        raise ARCSIException("Not Implemented")

    def generateImageSaturationMask(self, outputPath, outputName, outFormat):
        print("Generate Saturation Image")
        raise ARCSIException("Not Implemented")

    def generateValidImageDataMask(self, outputPath, outputMaskName, viewAngleImg, outFormat):
        print("Generate valid image mask")
        outputImage = os.path.join(outputPath, outputMaskName)
        return outputImage

    def convertThermalToBrightness(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        raise ARCSIException("Not Implemented")

    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat, scaleFactor):
        print("Converting to TOA")
        raise ARCSIException("Not Implemented")

    def generateCloudMask(self, inputReflImage, inputSatImage, inputThermalImage, inputValidImg, outputPath, outputName, outFormat, tmpPath, scaleFactor):
        raise ARCSIException("Cloud Masking Not Implemented for Sentinel-2.")

    def createCloudMaskDataArray(self, inImgDataArr):
        return inImgDataArr

    def calc6SCoefficients(self, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF):
        sixsCoeffs = numpy.zeros((12, 6), dtype=numpy.float32)
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

        # Band 1
        s.wavelength = Py6S.Wavelength(0.430, 0.4575, [0.015297, 0.19593, 0.511598, 0.587385, 0.699199, 0.792047, 0.993878, 0.990178, 0.956726, 0.723933, 0.051814, 0.000303])
        s.run()
        sixsCoeffs[0,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[0,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[0,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[0,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[0,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[0,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 2
        s.wavelength = Py6S.Wavelength(0.440, 0.535, [0.001206, 0.002623, 0.002076, 0.002224, 0.002377, 0.002856, 0.009028, 0.038955, 0.292197, 0.382418, 0.400158, 0.424686, 0.505323, 0.529543, 0.534656, 0.543691, 0.601967, 0.621092, 0.575863, 0.546131, 0.571684, 0.633236, 0.738396, 0.768325, 0.788363, 0.809151, 0.844983, 0.840111, 0.78694, 0.761923, 0.810031, 0.901671, 1.0, 0.908308, 0.286992, 0.102833, 0.02508, 0.002585, 0.000441])
        s.run()
        sixsCoeffs[1,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[1,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[1,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[1,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[1,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[1,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 3
        s.wavelength = Py6S.Wavelength(0.5375, 0.5825, [0.00084, 0.080665, 0.341374, 0.828036, 0.888565, 0.860271, 0.834035, 0.867734, 0.933938, 1.0, 0.981107, 0.868656, 0.81291, 0.789606, 0.830458, 0.85799, 0.62498, 0.098293, 0.016512])
        s.run()
        sixsCoeffs[2,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[2,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[2,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[2,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[2,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[2,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 4
        s.wavelength = Py6S.Wavelength(0.6475, 0.6825, [0.034529, 0.817746, 0.983869, 0.995449, 0.977215, 0.814166, 0.764864, 0.830828, 0.883581, 0.955931, 0.973219, 0.965712, 0.944811, 0.422967, 0.063172])
        s.run()
        sixsCoeffs[3,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[3,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[3,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[3,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[3,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[3,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 5
        s.wavelength = Py6S.Wavelength(0.695, 0.7125, [0.04126, 0.478496, 1, 0.993239, 0.945953, 0.902399, 0.757197, 0.196706])
        s.run()
        sixsCoeffs[4,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[4,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[4,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[4,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[4,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[4,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 6
        s.wavelength = Py6S.Wavelength(0.7325, 0.7475, [0.085006, 0.920265, 0.934211, 0.981932, 0.993406, 0.962584, 0.506722])
        s.run()
        sixsCoeffs[5,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[5,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[5,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[5,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[5,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[5,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 7
        s.wavelength = Py6S.Wavelength(0.770, 0.7975, [0.014731, 0.199495, 0.898494, 0.994759, 0.964657, 0.846898, 0.777241, 0.800984, 0.757695, 0.536855, 0.077219, 0.003152])
        s.run()
        sixsCoeffs[6,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[6,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[6,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[6,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[6,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[6,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 8
        s.wavelength = Py6S.Wavelength(0.775, 0.9075, [0.019072, 0.056536, 0.203436, 0.450085, 0.81829, 0.960732, 0.985213, 0.93655, 0.941281, 0.962183, 0.959009, 0.945147, 0.945357, 0.937084, 0.900979, 0.86216, 0.801819, 0.755632, 0.708669, 0.690211, 0.682649, 0.67595, 0.660812, 0.65831, 0.685501, 0.720686, 0.776608, 0.78772, 0.776161, 0.759264, 0.720589, 0.69087, 0.649339, 0.627424, 0.604322, 0.591724, 0.581202, 0.580197, 0.589481, 0.596749, 0.605476, 0.613463, 0.637436, 0.659233, 0.659924, 0.615841, 0.526407, 0.49653, 0.529093, 0.537964, 0.326791, 0.14854, 0.033246, 0.007848])
        s.run()
        sixsCoeffs[7,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[7,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[7,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[7,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[7,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[7,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 9
        s.wavelength = Py6S.Wavelength(0.850, 0.880, [0.02471, 0.104944, 0.585731, 0.87843, 0.926043, 0.935962, 0.965458, 0.97988, 0.988474, 0.999626, 0.472189, 0.106955, 0.008819])
        s.run()
        sixsCoeffs[8,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[8,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[8,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[8,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[8,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[8,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 10
        s.wavelength = Py6S.Wavelength(0.9325, 0.9575, [0.018022, 0.408108, 0.873658, 0.983566, 0.996767, 0.998123, 1.0, 0.956408, 0.931094, 0.450443, 0.059807])
        s.run()
        sixsCoeffs[9,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[9,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[9,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[9,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[9,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[9,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 11
        s.wavelength = Py6S.Wavelength(1.540, 1.6825, [7.00E-06, 2.80E-05, 0.000147, 0.00048, 0.000911, 0.001684, 0.005345, 0.012628, 0.039584, 0.07493, 0.182597, 0.330736, 0.647173, 0.815215, 0.88703, 0.891417, 0.916528, 0.935322, 0.951416, 0.956429, 0.96348, 0.96818, 0.975915, 0.979878, 0.981412, 0.980705, 0.982736, 0.987807, 0.993288, 0.990405, 0.980023, 0.972568, 0.966371, 0.96605, 0.973463, 0.983472, 0.995476, 0.998568, 0.998804, 0.99973, 0.999814, 0.99162, 0.969903, 0.953287, 0.938586, 0.928114, 0.82498, 0.641891, 0.32371, 0.163972, 0.046194, 0.019359, 0.006523, 0.003409, 0.001423, 0.000498, 3.40E-05, 1.30E-05])
        s.run()
        sixsCoeffs[10,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[10,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[10,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[10,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[10,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[10,5] = float(s.outputs.values['environmental_irradiance'])

        # Band 12
        s.wavelength = Py6S.Wavelength(2.080, 2.320, [0.002885, 0.006597, 0.00854, 0.010002, 0.013364, 0.017126, 0.027668, 0.040217, 0.073175, 0.11147, 0.203461, 0.284898, 0.408003, 0.476537, 0.543352, 0.568634, 0.598891, 0.621362, 0.663707, 0.696165, 0.741301, 0.772071, 0.809677, 0.828599, 0.851107, 0.854746, 0.859532, 0.863257, 0.869696, 0.878588, 0.889473, 0.896696, 0.904831, 0.905665, 0.904783, 0.903347, 0.901983, 0.904313, 0.908092, 0.91295, 0.921302, 0.927219, 0.934142, 0.937086, 0.937652, 0.942518, 0.942117, 0.938428, 0.933022, 0.921057, 0.908293, 0.908191, 0.922855, 0.919482, 0.924526, 0.931974, 0.946802, 0.954437, 0.962539, 0.966042, 0.96546, 0.963656, 0.957327, 0.953558, 0.951731, 0.952641, 0.960639, 0.968307, 0.982898, 0.990734, 0.998753, 0.999927, 0.993884, 0.983735, 0.958343, 0.938203, 0.905999, 0.881683, 0.84062, 0.809516, 0.749107, 0.688185, 0.566031, 0.474659, 0.342092, 0.263176, 0.16809, 0.124831, 0.082363, 0.062691, 0.042864, 0.034947, 0.027418, 0.023959, 0.016331, 0.007379, 0.002065])
        s.run()
        sixsCoeffs[11,0] = float(s.outputs.values['coef_xa'])
        sixsCoeffs[11,1] = float(s.outputs.values['coef_xb'])
        sixsCoeffs[11,2] = float(s.outputs.values['coef_xc'])
        sixsCoeffs[11,3] = float(s.outputs.values['direct_solar_irradiance'])
        sixsCoeffs[11,4] = float(s.outputs.values['diffuse_solar_irradiance'])
        sixsCoeffs[11,5] = float(s.outputs.values['environmental_irradiance'])

        return sixsCoeffs

    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF, scaleFactor):
        print("Converting to Surface Reflectance")
        raise ARCSIException("Not Implemented")

    def convertImageToSurfaceReflDEMElevLUT(self, inputRadImage, inputDEMFile, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, aotVal, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, scaleFactor, elevCoeffs=None):
        print("Converting to Surface Reflectance")
        raise ARCSIException("Not Implemented")

    def convertImageToSurfaceReflAOTDEMElevLUT(self, inputRadImage, inputDEMFile, inputAOTImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, useBRDF, surfaceAltitudeMin, surfaceAltitudeMax, aotMin, aotMax, scaleFactor, elevAOTCoeffs=None):
        print("Converting to Surface Reflectance")
        raise ARCSIException("Not Implemented")

    def run6SToOptimiseAODValue(self, aotVal, radBlueVal, predBlueVal, aeroProfile, atmosProfile, grdRefl, surfaceAltitude):
        """Used as part of the optimastion for identifying values of AOD"""
        raise ARCSIException("Not Implemented")

    def findDDVTargets(self, inputTOAImage, outputPath, outputName, outFormat, tmpPath):
        print("Not implemented\n")
        sys.exit()

    def estimateImageToAODUsingDDV(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax):
        print("Not implemented\n")
        sys.exit()

    def estimateImageToAODUsingDOS(self, inputRADImage, inputTOAImage, inputDEMFile, shadowMask, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, aotValMin, aotValMax, globalDOS, simpleDOS, dosOutRefl):
        raise ARCSIException("Not Implemented")

    def estimateSingleAOTFromDOS(self, radianceImage, toaImage, inputDEMFile, tmpPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, minAOT, maxAOT, dosOutRefl):
        raise ARCSIException("Not Implemented")

    def setBandNames(self, imageFile):
        dataset = gdal.Open(imageFile, gdal.GA_Update)
        if not dataset is None:
            dataset.GetRasterBand(1).SetDescription("Coastal")
            dataset.GetRasterBand(2).SetDescription("Blue")
            dataset.GetRasterBand(3).SetDescription("Green")
            dataset.GetRasterBand(4).SetDescription("Yellow")
            dataset.GetRasterBand(5).SetDescription("Red")
            dataset.GetRasterBand(6).SetDescription("RedEdge")
            dataset.GetRasterBand(7).SetDescription("NIR1")
            dataset.GetRasterBand(8).SetDescription("NIR2")
            dataset = None
        else:
            print("Could not open image to set band names: ", imageFile)


