import numpy as np
from MDAnalysis.core.topologyattrs import ResidueAttr

class LipidHeadgroup(ResidueAttr):
    attrname = "lipid_headgroups"
    singular = "lipid_headgroup"

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(nr)], dtype=object)


class LipidClass(ResidueAttr):
    attrname = "lipid_classes"
    singular = "lipid_class"

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(nr)], dtype=object)


class LipidNumUnsaturations(ResidueAttr):
    attrname = "lipid_num_unsaturations"
    singular = "lipid_num_unsaturation"

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.ones(nr, dtype=int) * -1


class LipidSaturation(ResidueAttr):
    attrname = "lipid_saturations"
    singular = "lipid_saturation"

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.array(['' for _ in range(nr)], dtype=object)