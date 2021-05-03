"""
lipyds
A toolkit for leaflet-based membrane analysis
"""

# Add imports here
from .data import (LIPID_HEADGROUPS,
                   LIPID_CLASSES,
                   LIPID_NUM_UNSATURATIONS,
                   LIPID_SATURATIONS)
from .leafletfinder.leafletfinder import LeafletFinder
from .analysis import AreaPerLipid, LipidEnrichment, LipidFlipFlop

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
