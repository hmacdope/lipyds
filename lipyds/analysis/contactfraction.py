"""
Lipid Depletion-Enrichment Index
================================

Classes
-------

.. autoclass:: LipidEnrichment
    :members:

"""

from typing import Union

import scipy
import numpy as np

from MDAnalysis.core.universe import Universe
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.analysis import distances
from .base import LeafletAnalysisBase
from ..leafletfinder.utils import get_centers_by_residue

class ContactFraction(LeafletAnalysisBase):
    def __init__(self, universe: Union[AtomGroup, Universe], 
                 cutoff: float=12,
                 **kwargs):
        super().__init__(universe, **kwargs)
        self.cutoff = cutoff
        

    def _prepare(self):
        shape = (self.n_leaflets, self.n_unique_ids, self.n_unique_ids, self.n_frames)
        shape2 = (self.n_leaflets, self.n_unique_ids, self.n_frames)
        self.residue_residue_contact_counts = np.zeros(shape, dtype=int)
        self.residue_contact_counts = np.zeros(shape2, dtype=int)
        self.residue_counts = np.zeros(shape2, dtype=int)
        self.total_lipid_counts = np.zeros((self.n_leaflets, self.n_frames), dtype=int)
        self.timewise_cfraction = np.empty(shape)
        self.timewise_cfraction[:] = np.nan

    def _single_frame(self):
        fr_i = self._frame_index
        for lf_i in range(self.n_leaflets):
            ag = self.leaflet_atomgroups[lf_i]
            # coords = get_centers_by_residue(ag)
            ag_ids = getattr(ag, self.group_by_attr)
            ag_resix = ag.resindices
            a, b = distances.capped_distance(ag.positions, ag.positions,
                                             self.cutoff,
                                             box=self.get_box(),
                                             return_distances=False).T
            not_self = a != b
            a = a[not_self]
            b = b[not_self]

            not_same_res = []
            seen = set()
            for res_a, res_b in zip(ag_resix[a], ag_resix[b]):
                pair = (res_a, res_b)
                not_same_res.append(pair not in seen)
                seen.add(pair)
            
            a = a[not_same_res]
            b = b[not_same_res]
            self.total_lipid_counts[lf_i, fr_i] = len(ag.residues)
            res_res_row = self.residue_residue_contact_counts[lf_i]
            for id_i, resid in enumerate(self.unique_ids):
                a_neighbors = ag_ids[a] == resid
                n_contacts = len(set(ag.resindices[b][a_neighbors]))
                self.residue_contact_counts[lf_i, id_i, fr_i] = n_contacts
                self.residue_counts[lf_i, id_i, fr_i] = n_res = sum(ag_ids == resid)
                if not n_res or not n_contacts:
                    continue
                
                for id_j, resid2 in enumerate(self.unique_ids):
                    b_neighbors = ag_ids[b] == resid2
                    b_neighboring_a = a_neighbors & b_neighbors
                    unique_b_neighbors = set(ag.resindices[b][b_neighboring_a])
                    res_res_row[id_i, id_j, fr_i] = len(unique_b_neighbors)

                    local_ratio = len(unique_b_neighbors) / n_contacts
                    global_ratio = sum(ag_ids == resid2) / len(ag_ids)
                    self.timewise_cfraction[lf_i, id_i, id_j, fr_i] = local_ratio / global_ratio

    def _conclude(self):
        self.contact_fraction_by_leaflet = np.ones((self.n_leaflets, self.n_unique_ids, self.n_unique_ids))

        n_specific_contacts = self.residue_residue_contact_counts.mean(axis=-1)
        n_contacts = self.residue_contact_counts.mean(axis=-1)
        n_residues = self.residue_counts.mean(axis=-1)
        n_in_leaflet = self.total_lipid_counts.mean(axis=-1)

        for lf_i in range(self.n_leaflets):
            specific = n_specific_contacts[lf_i]
            contacts = n_contacts[lf_i]
            residues = n_residues[lf_i]
            leaflet = n_in_leaflet[lf_i]
            for i in range(self.n_unique_ids):
                for j in range(self.n_unique_ids):
                    local_ratio = specific[i, j] / contacts[j]
                    global_ratio = residues[i] / leaflet
                    if global_ratio:

                        self.contact_fraction_by_leaflet[lf_i, i, j] = local_ratio/global_ratio


    def collate_as_dataframe(self, ids=None):
        """Convert the results summary into a pandas DataFrame.

        This requires pandas to be installed.
        """

        if not hasattr(self, "contact_fraction_by_leaflet"):
            raise ValueError('Call run() first to get results')
        try:
            import pandas as pd
        except ImportError:
            raise ImportError('pandas is required to use this function '
                              'but is not installed. Please install with '
                              '`conda install pandas` or '
                              '`pip install pandas`.') from None

        dfs = [pd.DataFrame(d, columns=self.unique_ids)
               for d in self.contact_fraction_by_leaflet]
        for df in dfs:
            df["Y"] = self.unique_ids
        
        if ids is not None:
            ids = list(ids) + ["Y"]
            dfs = [df[df.Y.isin(ids)] for df in dfs]
            dfs = [df[ids] for df in dfs]
        for i, df in enumerate(dfs, 1):
            df['Leaflet'] = i
        df = pd.concat(dfs)
        return df

    def average_as_dataframe(self, ids=None):
        """Convert the results summary into a pandas DataFrame.

        This requires pandas to be installed.
        """

        if not hasattr(self, "timewise_cfraction"):
            raise ValueError('Call run() first to get results')
        try:
            import pandas as pd
        except ImportError:
            raise ImportError('pandas is required to use this function '
                              'but is not installed. Please install with '
                              '`conda install pandas` or '
                              '`pip install pandas`.') from None

        means = np.nanmean(self.timewise_cfraction, axis=-1)
        dfs = [pd.DataFrame(d, columns=self.unique_ids)
               for d in means]
        for df in dfs:
            df["Y"] = self.unique_ids
        
        if ids is not None:
            ids = list(ids)
            dfs = [df[df.Y.isin(ids)] for df in dfs]
            rows = []
            for df in dfs:
                df = df.set_index("Y")
                df = df.reindex(ids)
                rows.append(df)
            dfs = [df[ids] for df in rows]
        for i, df in enumerate(dfs, 1):
            df['Leaflet'] = i
        df = pd.concat(dfs)
        return df



class SymmetricContactFraction(ContactFraction):

    def _prepare(self):
        shape = (self.n_leaflets, self.n_unique_ids, self.n_unique_ids, self.n_frames)
        shape2 = (self.n_leaflets, self.n_unique_ids, self.n_frames)
        self.residue_residue_contact_counts = np.zeros(shape, dtype=int)
        self.residue_contact_counts = np.zeros((self.n_leaflets, self.n_frames), dtype=int)
        self.residue_counts = np.zeros(shape2, dtype=int)
        self.total_lipid_counts = np.zeros((self.n_leaflets, self.n_frames), dtype=int)
        self.timewise_cfraction = np.empty(shape)
        self.timewise_cfraction[:] = np.nan

    def _single_frame(self):
        fr_i = self._frame_index
        for lf_i in range(self.n_leaflets):
            ag = self.leaflet_atomgroups[lf_i]
            ag_ids = getattr(ag, self.group_by_attr)
            ag_resix = ag.resindices
            a, b = distances.capped_distance(ag.positions, ag.positions,
                                             self.cutoff,
                                             box=self.get_box(),
                                             return_distances=False).T
            not_self = a != b
            a = a[not_self]
            b = b[not_self]

            not_same_res = []
            seen = set()
            for res_a, res_b in zip(ag_resix[a], ag_resix[b]):
                pair = (res_a, res_b)
                not_same_res.append(pair not in seen)
                seen.add(pair)
            
            a = a[not_same_res]
            b = b[not_same_res]

            self.total_lipid_counts[lf_i, fr_i] = n_tot = len(ag.residues)
            self.residue_contact_counts[lf_i, fr_i] = n_contacts = len(a)

            res_res_row = self.residue_residue_contact_counts[lf_i]
            for id_i, resid in enumerate(self.unique_ids):
                a_neighbors = ag_ids[a] == resid
                
                self.residue_counts[lf_i, id_i, fr_i] = n_res = sum(ag_ids == resid)
                if not n_res or not n_contacts:
                    continue
                
                for id_j, resid2 in enumerate(self.unique_ids[id_i:], id_i):
                    b_neighbors = ag_ids[b] == resid2
                    b_neighboring_a = a_neighbors & b_neighbors
                    res_res_row[id_i, id_j, fr_i] = n_ab = sum(b_neighboring_a)

                    local_ratio = n_ab / n_contacts
                    prob_b = sum(ag_ids == resid2) / n_tot
                    prob_a = n_res / n_tot
                    global_ratio = prob_a * prob_b
                    if global_ratio:
                        self.timewise_cfraction[lf_i, id_i, id_j, fr_i] = local_ratio / global_ratio

    def _conclude(self):
        pass