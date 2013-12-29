"""
Module that contains the ARCSIAbstractSensor class.
"""
############################################################################
#  arcsisensor.py
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
# Purpose:  An abstract class which provides the base class for each 
#           sensor supported by the ARCSI system.
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

# import abstract base class stuff
from abc import ABCMeta, abstractmethod
# Import the datetime module
import datetime
# Import the ARCSI exception class
from .arcsiexception import ARCSIException


class ARCSIAbstractSensor (object, metaclass=ABCMeta):
    """
    An abstract class which represents a sensor and allows
    the various opperations required to be applied and standard
    variables (e.g., acqusiation date) stored and retrieved.
    """
        
    def __init__(self):
        self.sensor = "NA"
        self.acquisitionTime = datetime.datetime.today()
        self.latTL = 0.0
        self.lonTL = 0.0
        self.latTR = 0.0
        self.lonTR = 0.0
        self.latBL = 0.0
        self.lonBL = 0.0
        self.latBR = 0.0
        self.lonBR = 0.0
        self.latCentre = 0.0
        self.lonCentre = 0.0
        self.xTL = 0.0
        self.yTL = 0.0
        self.xTR = 0.0
        self.yTR = 0.0
        self.xBL = 0.0
        self.yBL = 0.0
        self.xBR = 0.0
        self.yBR = 0.0
        self.xCentre = 0.0
        self.yCentre = 0.0
        self.inWKT = ""
        self.solarZenith = 0.0
        self.solarAzimuth = 0.0
        self.senorZenith = 0.0
        self.senorAzimuth = 0.0
        self.epsgCodes = dict()
        self.epsgCodes["WGS84UTM1N"] = 32601
        self.epsgCodes["WGS84UTM2N"] = 32602
        self.epsgCodes["WGS84UTM3N"] = 32603
        self.epsgCodes["WGS84UTM4N"] = 32604
        self.epsgCodes["WGS84UTM5N"] = 32605
        self.epsgCodes["WGS84UTM6N"] = 32606
        self.epsgCodes["WGS84UTM7N"] = 32607
        self.epsgCodes["WGS84UTM8N"] = 32608
        self.epsgCodes["WGS84UTM9N"] = 32609
        self.epsgCodes["WGS84UTM10N"] = 32610
        self.epsgCodes["WGS84UTM11N"] = 32611
        self.epsgCodes["WGS84UTM12N"] = 32612
        self.epsgCodes["WGS84UTM13N"] = 32613
        self.epsgCodes["WGS84UTM14N"] = 32614
        self.epsgCodes["WGS84UTM15N"] = 32615
        self.epsgCodes["WGS84UTM16N"] = 32616
        self.epsgCodes["WGS84UTM17N"] = 32617
        self.epsgCodes["WGS84UTM18N"] = 32618
        self.epsgCodes["WGS84UTM19N"] = 32619
        self.epsgCodes["WGS84UTM20N"] = 32620
        self.epsgCodes["WGS84UTM21N"] = 32621
        self.epsgCodes["WGS84UTM22N"] = 32622
        self.epsgCodes["WGS84UTM23N"] = 32623
        self.epsgCodes["WGS84UTM24N"] = 32624
        self.epsgCodes["WGS84UTM25N"] = 32625
        self.epsgCodes["WGS84UTM26N"] = 32626
        self.epsgCodes["WGS84UTM27N"] = 32627
        self.epsgCodes["WGS84UTM28N"] = 32628
        self.epsgCodes["WGS84UTM29N"] = 32629
        self.epsgCodes["WGS84UTM30N"] = 32630
        self.epsgCodes["WGS84UTM31N"] = 32631
        self.epsgCodes["WGS84UTM32N"] = 32632
        self.epsgCodes["WGS84UTM33N"] = 32633
        self.epsgCodes["WGS84UTM34N"] = 32634
        self.epsgCodes["WGS84UTM35N"] = 32635
        self.epsgCodes["WGS84UTM36N"] = 32636
        self.epsgCodes["WGS84UTM37N"] = 32637
        self.epsgCodes["WGS84UTM38N"] = 32638
        self.epsgCodes["WGS84UTM39N"] = 32639
        self.epsgCodes["WGS84UTM40N"] = 32640
        self.epsgCodes["WGS84UTM41N"] = 32641
        self.epsgCodes["WGS84UTM42N"] = 32642
        self.epsgCodes["WGS84UTM43N"] = 32643
        self.epsgCodes["WGS84UTM44N"] = 32644
        self.epsgCodes["WGS84UTM45N"] = 32645
        self.epsgCodes["WGS84UTM46N"] = 32646
        self.epsgCodes["WGS84UTM47N"] = 32647
        self.epsgCodes["WGS84UTM48N"] = 32648
        self.epsgCodes["WGS84UTM49N"] = 32649
        self.epsgCodes["WGS84UTM50N"] = 32650
        self.epsgCodes["WGS84UTM51N"] = 32651
        self.epsgCodes["WGS84UTM52N"] = 32652
        self.epsgCodes["WGS84UTM53N"] = 32653
        self.epsgCodes["WGS84UTM54N"] = 32654
        self.epsgCodes["WGS84UTM55N"] = 32655
        self.epsgCodes["WGS84UTM56N"] = 32656
        self.epsgCodes["WGS84UTM57N"] = 32657
        self.epsgCodes["WGS84UTM58N"] = 32658
        self.epsgCodes["WGS84UTM59N"] = 32659
        self.epsgCodes["WGS84UTM60N"] = 32660
        self.epsgCodes["WGS84UTM1S"] = 32701
        self.epsgCodes["WGS84UTM2S"] = 32702
        self.epsgCodes["WGS84UTM3S"] = 32703
        self.epsgCodes["WGS84UTM4S"] = 32704
        self.epsgCodes["WGS84UTM5S"] = 32705
        self.epsgCodes["WGS84UTM6S"] = 32706
        self.epsgCodes["WGS84UTM7S"] = 32707
        self.epsgCodes["WGS84UTM8S"] = 32708
        self.epsgCodes["WGS84UTM9S"] = 32709
        self.epsgCodes["WGS84UTM10S"] = 32710
        self.epsgCodes["WGS84UTM11S"] = 32711
        self.epsgCodes["WGS84UTM12S"] = 32712
        self.epsgCodes["WGS84UTM13S"] = 32713
        self.epsgCodes["WGS84UTM14S"] = 32714
        self.epsgCodes["WGS84UTM15S"] = 32715
        self.epsgCodes["WGS84UTM16S"] = 32716
        self.epsgCodes["WGS84UTM17S"] = 32717
        self.epsgCodes["WGS84UTM18S"] = 32718
        self.epsgCodes["WGS84UTM19S"] = 32719
        self.epsgCodes["WGS84UTM20S"] = 32720
        self.epsgCodes["WGS84UTM21S"] = 32721
        self.epsgCodes["WGS84UTM22S"] = 32722
        self.epsgCodes["WGS84UTM23S"] = 32723
        self.epsgCodes["WGS84UTM24S"] = 32724
        self.epsgCodes["WGS84UTM25S"] = 32725
        self.epsgCodes["WGS84UTM26S"] = 32726
        self.epsgCodes["WGS84UTM27S"] = 32727
        self.epsgCodes["WGS84UTM28S"] = 32728
        self.epsgCodes["WGS84UTM29S"] = 32729
        self.epsgCodes["WGS84UTM30S"] = 32730
        self.epsgCodes["WGS84UTM31S"] = 32731
        self.epsgCodes["WGS84UTM32S"] = 32732
        self.epsgCodes["WGS84UTM33S"] = 32733
        self.epsgCodes["WGS84UTM34S"] = 32734
        self.epsgCodes["WGS84UTM35S"] = 32735
        self.epsgCodes["WGS84UTM36S"] = 32736
        self.epsgCodes["WGS84UTM37S"] = 32737
        self.epsgCodes["WGS84UTM38S"] = 32738
        self.epsgCodes["WGS84UTM39S"] = 32739
        self.epsgCodes["WGS84UTM40S"] = 32740
        self.epsgCodes["WGS84UTM41S"] = 32741
        self.epsgCodes["WGS84UTM42S"] = 32742
        self.epsgCodes["WGS84UTM43S"] = 32743
        self.epsgCodes["WGS84UTM44S"] = 32744
        self.epsgCodes["WGS84UTM45S"] = 32745
        self.epsgCodes["WGS84UTM46S"] = 32746
        self.epsgCodes["WGS84UTM47S"] = 32747
        self.epsgCodes["WGS84UTM48S"] = 32748
        self.epsgCodes["WGS84UTM49S"] = 32749
        self.epsgCodes["WGS84UTM50S"] = 32750
        self.epsgCodes["WGS84UTM51S"] = 32751
        self.epsgCodes["WGS84UTM52S"] = 32752
        self.epsgCodes["WGS84UTM53S"] = 32753
        self.epsgCodes["WGS84UTM54S"] = 32754
        self.epsgCodes["WGS84UTM55S"] = 32755
        self.epsgCodes["WGS84UTM56S"] = 32756
        self.epsgCodes["WGS84UTM57S"] = 32757
        self.epsgCodes["WGS84UTM58S"] = 32758
        self.epsgCodes["WGS84UTM59S"] = 32759
        self.epsgCodes["WGS84UTM60S"] = 32760
        
        
    
    @abstractmethod
    def extractHeaderParameters(self, inputHeader, wktStr): pass
    
    def defaultGenBaseOutFileName(self):
        """
        A function to generate a generic standard file
        base name which will be sensible.
        
        It is expected that individual sensors may override this function.
        """
        date = self.acquisitionTime.strftime("%Y%m%d")
        pos = "lat" + str(round(self.latCentre,)).replace('.', '').replace('-', '') + "lon" + str(round(self.lonCentre,2)).replace('.', '').replace('-', '')
        outname = self.sensor + "_" + date + "_" + pos
        return outname
        
    def generateOutputBaseName(self):
        """
        Provides a default implementation.
        """
        return self.defaultGenBaseOutFileName()
    
    @abstractmethod
    def convertImageToRadiance(self, outputPath, outputName, outFormat): pass
    
    @abstractmethod
    def convertImageToTOARefl(self, inputRadImage, outputPath, outputName, outFormat): pass
    
    @abstractmethod
    def convertImageToSurfaceReflSglParam(self, inputRadImage, outputPath, outputName, outFormat, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotVal, useBRDF): pass

    @abstractmethod
    def estimateImageToAOD(self, inputRADImage, inputTOAImage, outputPath, outputName, outFormat, tmpPath, aeroProfile, atmosProfile, grdRefl, surfaceAltitude, aotValMin, aotValMax): pass
        


