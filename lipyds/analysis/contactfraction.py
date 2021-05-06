"""
Lipid Depletion-Enrichment Index
================================

Classes
-------

.. autoclass:: LipidEnrichment
    :members:

"""
import functools
from typing import Union
from collections import defaultdict

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
        self.id_to_index = {x: i for i, x in enumerate(self.unique_ids)}
        

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
            # grab slices we're modifying
            res_n_contacts = self.residue_contact_counts[lf_i, :, fr_i]
            res_counts = self.residue_counts[lf_i, :, fr_i]
            res_res_contacts = self.residue_residue_contact_counts[lf_i, :, :, fr_i]
            timewise_cfrac = self.timewise_cfraction[lf_i, :, :, fr_i]

            ag = self.leaflet_atomgroups[lf_i]

            # num residues in leaflet
            self.total_lipid_counts[lf_i, fr_i] = n_lipids = len(ag.residues)

            resids = getattr(ag.residues, self.group_by_attr)
            unique_resids, counts = np.unique(resids, return_counts=True)
            for rid, count in zip(unique_resids, counts):
                res_counts[self.id_to_index[rid]] = count

            # get atom index neighbors
            a, b = distances.capped_distance(ag.positions, ag.positions,
                                             self.cutoff,
                                             box=self.get_box(),
                                             return_distances=False).T
            a_atoms = ag[a]
            b_atoms = ag[b]

            # get residue neighbors
            ab_resindices = set(x for x in zip(a_atoms.resindices, b_atoms.resindices) if x[0] != x[1])
            a_resindices, b_resindices = map(list, zip(*ab_resindices))
            
            a_residues = self.universe.residues[a_resindices]
            b_residues = self.universe.residues[b_resindices]

            a_ids = getattr(a_residues, self.group_by_attr)
            b_ids = getattr(b_residues, self.group_by_attr)

            unique_a, a_counts = np.unique(a_ids, return_counts=True)
            for a_id, count in zip(unique_a, a_counts):
                res_n_contacts[self.id_to_index[a_id]] = count

            # sort into pairs
            pair_counts = defaultdict(int)
            # this map construction looks stupid but is ~10x faster than
            # a for loop
            def count_pairs(a, b):
                pair_counts[(a, b)] += 1
            list(map(count_pairs, a_ids, b_ids))

            n_contacts = sum(res_n_contacts)

            func = functools.partial(self._calculate_timewise_cfraction,
                                     pair_counts=pair_counts,
                                     res_n_contacts=res_n_contacts,
                                     res_counts=res_counts,
                                     res_res_contacts=res_res_contacts,
                                     timewise_cfrac=timewise_cfrac,
                                     n_interactions=len(a_resindices),
                                     n_lipids=n_lipids)

            list(map(func, range(self.n_unique_ids)))

    def _calculate_timewise_cfraction(self, id_i, pair_counts={},
                                      res_n_contacts=[], res_counts=[],
                                      res_res_contacts=[], timewise_cfrac=[],
                                      n_interactions=0, n_lipids=0):
            def row(id_j):
                resi, resj = self.unique_ids[[id_i, id_j]]
                res_res_contacts[id_i, id_j] = n_ab = pair_counts[(resi, resj)]
                local_ratio = n_ab / res_n_contacts[id_i]
                global_ratio = res_counts[id_j] / n_lipids
                timewise_cfrac[id_i, id_j] = local_ratio / global_ratio

            list(map(row, range(self.n_unique_ids)))


    # def _conclude(self):
    #     self.contact_fraction_by_leaflet = np.ones((self.n_leaflets, self.n_unique_ids, self.n_unique_ids))

    #     n_specific_contacts = self.residue_residue_contact_counts.mean(axis=-1)
    #     n_contacts = self.residue_contact_counts.mean(axis=-1)
    #     n_residues = self.residue_counts.mean(axis=-1)
    #     n_in_leaflet = self.total_lipid_counts.mean(axis=-1)

    #     for lf_i in range(self.n_leaflets):
    #         specific = n_specific_contacts[lf_i]
    #         contacts = n_contacts[lf_i]
    #         residues = n_residues[lf_i]
    #         leaflet = n_in_leaflet[lf_i]
    #         for i in range(self.n_unique_ids):
    #             for j in range(self.n_unique_ids):
    #                 local_ratio = specific[i, j] / contacts[j]
    #                 global_ratio = residues[i] / leaflet
    #                 if global_ratio:

    #                     self.contact_fraction_by_leaflet[lf_i, i, j] = local_ratio/global_ratio


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
            df.index = self.unique_ids
        
        if ids is not None:
            ids = list(ids)
            dfs = [df[df.index.isin(ids)] for df in dfs]
            rows = []
            for df in dfs:
                df = df.reindex(ids)
                rows.append(df)
            dfs = [df[ids] for df in rows]
        for i, df in enumerate(dfs, 1):
            df['Leaflet'] = i
        df = pd.concat(dfs)
        return df



class SymmetricContactFraction(ContactFraction):

    def _calculate_timewise_cfraction(self, id_i, pair_counts={},
                                      res_n_contacts=[], res_counts=[],
                                      res_res_contacts=[], timewise_cfrac=[],
                                      n_interactions=0, n_lipids=0):

            ratios = res_counts / n_lipids
            def row(id_j):
                resi, resj = self.unique_ids[[id_i, id_j]]
                n_ab = pair_counts[(resi, resj)]
                res_res_contacts[id_i, id_j] = res_res_contacts[id_j, id_i] = n_ab
                local_ratio = n_ab / n_interactions
                global_ratio = ratios[id_i] * ratios[id_j]
                cfrac = local_ratio / global_ratio
                timewise_cfrac[id_i, id_j] = timewise_cfrac[id_j, id_i] = cfrac

            list(map(row, range(id_i, self.n_unique_ids)))

    

    # def _conclude(self):
    #     self.contact_fraction_by_leaflet = np.ones((self.n_leaflets, self.n_unique_ids, self.n_unique_ids))

    #     n_specific_contacts = self.residue_residue_contact_counts.mean(axis=-1)
    #     n_contacts = self.residue_contact_counts.mean(axis=-1)
    #     n_residues = self.residue_counts.mean(axis=-1)
    #     n_in_leaflet = self.total_lipid_counts.mean(axis=-1)

    #     for lf_i in range(self.n_leaflets):
    #         specific = n_specific_contacts[lf_i]
    #         contacts = n_contacts[lf_i]
    #         residues = n_residues[lf_i]
    #         leaflet = n_in_leaflet[lf_i]
    #         for i in range(self.n_unique_ids):
    #             for j in range(self.n_unique_ids):
    #                 local_ratio = specific[i, j] / contacts[j]
    #                 global_ratio = residues[i] / leaflet
    #                 if global_ratio:

    #                     self.contact_fraction_by_leaflet[lf_i, i, j] = local_ratio/global_ratio