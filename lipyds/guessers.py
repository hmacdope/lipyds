from MDAnalysis.core.universe import Universe

from . import (LIPID_HEADGROUPS,
               LIPID_CLASSES,
               LIPID_NUM_UNSATURATIONS,
               LIPID_SATURATIONS)
from .topologyattrs import (LipidHeadgroup,
                            LipidClass,
                            LipidNumUnsaturations,
                            LipidSaturation)



def guess_lipid_info(universe: Universe,
                     lipid_headgroups=LIPID_HEADGROUPS,
                     lipid_classes=LIPID_CLASSES,
                     lipid_num_unsaturations=LIPID_NUM_UNSATURATIONS,
                     lipid_saturations=LIPID_SATURATIONS):
    classes = [LipidHeadgroup, LipidClass, LipidNumUnsaturations, LipidSaturation]
    defaults = {cls.attrname: "" for cls in classes}
    defaults[LipidNumUnsaturations.attrname] = -1
    values = dict(lipid_headgroups=lipid_headgroups,
                  lipid_classes=lipid_classes,
                  lipid_num_unsaturations=lipid_num_unsaturations,
                  lipid_saturations=lipid_saturations)

    resnames = universe.residues.resnames

    for attrname, default in defaults.items():
        guesses = [values[attrname].get(name, default) for name in resnames]
        if not hasattr(universe, attrname):
            universe.add_TopologyAttr(attrname)
        setattr(universe.residues, attrname, guesses)
