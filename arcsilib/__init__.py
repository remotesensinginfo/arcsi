"""
ARCSI - this file is needed to ensure it can be imported

See other source files for details
"""

from distutils.version import LooseVersion
import os

ARCSI_VERSION_MAJOR = 3
ARCSI_VERSION_MINOR = 1
ARCSI_VERSION_PATCH = 7

ARCSI_VERSION = str(ARCSI_VERSION_MAJOR) + "."  + str(ARCSI_VERSION_MINOR) + "." + str(ARCSI_VERSION_PATCH)
ARCSI_VERSION_OBJ = LooseVersion(ARCSI_VERSION)

ARCSI_COPYRIGHT_YEAR = "2017"
ARCSI_COPYRIGHT_NAMES = "Pete Bunting, Dan Clewley"

ARCSI_SUPPORT_EMAIL = "rsgislib-support@googlegroups.com"

ARCSI_WEBSITE = "http://www.rsgislib.org/arcsi"

ARCSI_SENSORS_LIST = ['ls1', 'ls2', 'ls3', 'ls4mss', 'ls4tm', 'ls5mss', 'ls5tm', 'ls7', 'ls8', 'sen2', 'rapideye', 'wv2', 'wv3', 'spot5', 'spot6', 'spot7', 'pleiades']
ARCSI_PRODUCTS_LIST = ['RAD', 'SATURATE', 'TOA', 'CLOUDS', 'CLEARSKY', 'DDVAOT', 'DOSAOT', 'DOSAOTSGL', 'SREF', 'STDSREF', 'DOS', 'THERMAL', 'TOPOSHADOW', 'FOOTPRINT', 'METADATA', 'SHARP']
ARCSI_ARCHIVE_EXE_LIST = [".tar.gz", ".tgz", ".TAR.GZ", ".TGZ", ".tar", ".TAR", ".zip", ".ZIP", ".tar.bz", ".TAR.BZ", ".tar.bz2", ".TAR.BZ2"]
ARCSI_GDALFORMATS_LIST = ['KEA']

# Get install prefix
install_prefix = __file__[:__file__.find('lib')]

DEFAULT_ARCSI_AEROIMG_PATH = os.path.join(install_prefix, "share","arcsi", "WorldAerosolParams.kea")
DEFAULT_ARCSI_ATMOSIMG_PATH = os.path.join(install_prefix,"share","arcsi", "WorldAtmosphereParams.kea")

# Check files exit - set to none if they don't
if not os.path.isfile(DEFAULT_ARCSI_ATMOSIMG_PATH):
   DEFAULT_ARCSI_ATMOSIMG_PATH = None
if not os.path.isfile(DEFAULT_ARCSI_AEROIMG_PATH):
   DEFAULT_ARCSI_AEROIMG_PATH = None

if os.environ.get('RIOS_DFLT_DRIVER', None) == None:
    os.environ['RIOS_DFLT_DRIVER'] = 'KEA'

