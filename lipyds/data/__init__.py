from .lipids import (get_lipid_heads_and_tails,
                    get_lipid_classes,
                    get_n_unsaturations,
                    get_saturations)


LIPID_HEADGROUPS, LIPID_TAIL_UNSATURATIONS = get_lipid_heads_and_tails()
LIPID_CLASSES = get_lipid_classes(LIPID_HEADGROUPS)
LIPID_NUM_UNSATURATIONS = get_n_unsaturations(LIPID_TAIL_UNSATURATIONS)
LIPID_SATURATIONS = get_saturations(LIPID_NUM_UNSATURATIONS)