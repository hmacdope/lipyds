from typing import Dict
from pkg_resources import resource_filename

import yaml

def get_lipid_heads_and_tails():
    head_file = resource_filename(__name__, "headgroups.yml")
    with open(head_file) as file:
        heads = yaml.load(file, Loader=yaml.FullLoader)
    tail_file = resource_filename(__name__, "n_saturations.yml")
    with open(tail_file) as file:
        tails = yaml.load(file, Loader=yaml.FullLoader)
    
    headgroups = heads.get("general", {})
    tailgroups = tails.get("general", {})

    LIPID_HEADGROUPS = {}
    LIPID_TAIL_UNSATURATIONS = {}

    for tail_name, tail in tailgroups.items():
        for head_name, head in headgroups.items():
            lipid = tail_name + head_name
            LIPID_HEADGROUPS[lipid] = head
            LIPID_TAIL_UNSATURATIONS[lipid] = tail
    
    for lipid, tail in tails.get("lipids", {}).items():
        LIPID_TAIL_UNSATURATIONS[lipid] = tail
    
    for lipid, head in heads.get("lipids", {}).items():
        LIPID_HEADGROUPS[lipid] = head
    
    return LIPID_HEADGROUPS, LIPID_TAIL_UNSATURATIONS


def get_lipid_classes(lipid_headgroups: Dict[str, str]) -> Dict[str, str]:
    classes = {"GM1": "GS",
               "GM3": "GS",
               "PIP1": "PI",
               "PIP2": "PI",
               "PIP3": "PI"}
    lipid_classes = {}
    for lipid, headgroup in lipid_headgroups.items():
        lipid_classes[lipid] = classes.get(headgroup, headgroup)
    return lipid_classes


def get_n_unsaturations(lipid_tail_unsaturations):
    n_unsaturations = {}
    for lipid, tail_unsaturations in lipid_tail_unsaturations.items():
        n_unsaturations[lipid] = sum(tail_unsaturations)
    return n_unsaturations


def get_saturations(lipid_n_unsaturations):
    saturations = {}
    for lipid, n_unsaturations in lipid_n_unsaturations.items():
        if n_unsaturations == 0:
            sat = "saturated"
        elif n_unsaturations == 1:
            sat = "monounsaturated"
        elif n_unsaturations > 1:
            sat = "polyunsaturated"
        else:
            continue
        saturations[lipid] = sat
    return saturations