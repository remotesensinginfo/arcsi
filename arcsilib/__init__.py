"""
ARCSI - this file is needed to ensure it can be imported

See other source files for details
"""

from distutils.version import LooseVersion
import os

ARCSI_VERSION_MAJOR = 0
ARCSI_VERSION_MINOR = 14
ARCSI_VERSION_PATCH = 5

ARCSI_VERSION = str(ARCSI_VERSION_MAJOR) + "."  + str(ARCSI_VERSION_MINOR) + "." + str(ARCSI_VERSION_PATCH)
ARCSI_VERSION_OBJ = LooseVersion(ARCSI_VERSION)

ARCSI_COPYRIGHT_YEAR = "2015"
ARCSI_COPYRIGHT_NAMES = "Pete Bunting"

ARCSI_SUPPORT_EMAIL = "rsgislib-support@googlegroups.com"

ARCSI_WEBSITE = "http://www.rsgislib.org/arcsi"

# Get install prefix
install_prefix = __file__[:__file__.find('lib')]

DEFAULT_ARCSI_AEROIMG_PATH = os.path.join(install_prefix, "share","arcsi", "WorldAerosolParams.kea")
DEFAULT_ARCSI_ATMOSIMG_PATH = os.path.join(install_prefix,"share","arcsi", "WorldAtmosphereParams.kea")

# Check files exit - set to none if they don't
if not os.path.isfile(DEFAULT_ARCSI_ATMOSIMG_PATH):
   DEFAULT_ARCSI_ATMOSIMG_PATH = None
if not os.path.isfile(DEFAULT_ARCSI_AEROIMG_PATH):
   DEFAULT_ARCSI_AEROIMG_PATH = None
