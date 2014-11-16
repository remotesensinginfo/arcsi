"""
ARCSI - this file is needed to ensure it can be imported

See other source files for details
"""

from distutils.version import LooseVersion

ARCSI_VERSION_MAJOR = 0
ARCSI_VERSION_MINOR = 13
ARCSI_VERSION_PATCH = 11

ARCSI_VERSION = str(ARCSI_VERSION_MAJOR) + "."  + str(ARCSI_VERSION_MINOR) + "." + str(ARCSI_VERSION_PATCH)
ARCSI_VERSION_OBJ = LooseVersion(ARCSI_VERSION)

ARCSI_COPYRIGHT_YEAR = "2014"
ARCSI_COPYRIGHT_NAMES = "Pete Bunting"

ARCSI_SUPPORT_EMAIL = "rsgislib-support@googlegroups.com"

ARCSI_WEBSITE = "http://www.rsgislib.org/arcsi"