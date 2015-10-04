from collections import defaultdict
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

from glycresoft_sqlalchemy.data_model import (Hypothesis, Protein, GlycopeptideMatch,
                                              Glycan, object_session)


class GlycoproteinMicroheterogeneitySummary(object):
    def __init__(self, protein, filter_fn):
        self.protein = protein
        self.site_list = protein.glycosylation_sites
        self.site_to_glycan_map = {}
        self.filter_fn = filter_fn
        self.build_site_to_glycan_map()

    def build_site_to_glycan_map(self):
        site_map = {}
        for site in self.site_list:
            species = defaultdict(list)
            scan_id_ranges = set()
            for start, end, glycan_composition, volume, ms2_score, scan_id_range in self.filter_fn(object_session(
                    self.protein).query(GlycopeptideMatch.start_position,
                                        GlycopeptideMatch.end_position,
                                        GlycopeptideMatch.glycan_composition_str,
                                        GlycopeptideMatch.volume,
                                        GlycopeptideMatch.ms2_score,
                                        GlycopeptideMatch.scan_id_range).filter(
                    GlycopeptideMatch.spans(site),
                    GlycopeptideMatch.protein_id == self.protein.id)):
                scan_id_range = tuple(scan_id_range)
                if scan_id_range in scan_id_ranges:
                    continue
                scan_id_ranges.add(scan_id_range)
                species[glycan_composition].append(volume)
            site_map[site] = {k: sum(v) for k, v in species.items()}
        self.site_to_glycan_map = site_map
        return site_map

    def glycan_site_specific_abundance_plot(self):
        axes = {}
        bar_width = 0.3
        for site, abundances in self.site_to_glycan_map.items():
            if len(abundances) == 0:
                continue
            compositions, volumes = zip(*abundances.items())
            fig, ax = plt.subplots(1)
            indices = np.arange(len(compositions))
            ax.bar(indices + bar_width, volumes, bar_width, alpha=0.45)
            ax.set_ylabel("Relative Intensity")
            ax.set_xticks(indices + bar_width * 1.5)
            ax.set_xticklabels(compositions, rotation=90, ha='center')
            ax.set_title("%s at %d" % (self.protein.name, site))
            axes[site] = ax

        return axes

    def glycan_peptide_abundance_plot(self):
        '''
        Plot the relative abundance of each glycoform for each peptide as a bar plot.
        '''
        peptide_bundles = []
        current_peptide = None
        current_bundle = []
        for row in self.filter_fn(self.protein.glycopeptide_matches.filter(
                GlycopeptideMatch.ms2_score > self.score_threshold).group_by(
                GlycopeptideMatch.base_peptide_sequence,
                GlycopeptideMatch.glycan_composition_str)):
            if row[0] != current_peptide:
                if current_peptide is not None:
                    peptide_bundles.append((current_peptide, current_bundle))
                current_bundle = []
                current_peptide = row[0]
            current_bundle.append((row[1], row[2]))
        peptide_bundles.append((current_peptide, current_bundle))
        axes = []
        bar_width = 0.3

        for bundle in peptide_bundles:
            fig, ax = plt.subplots(1)
            peptide, glycans = bundle
            indices = np.arange(len(glycans))
            glycans, volumes = zip(*glycans)
            ax.bar(indices + bar_width, volumes, bar_width, alpha=0.45)
            ax.set_ylabel("Relative Intensity")
            ax.set_xticks(indices + bar_width * 1.5)
            ax.set_xticklabels(glycans, rotation=90, ha='center')
            ax.set_title(peptide)
            axes.append(ax)
        return axes

    def __repr__(self):
        return "<GlycoproteinMicroheterogeneitySummary {} {}>".format(
            self.protein.name, {s: len(v) for s, v in self.site_to_glycan_map.items()})

    def __iter__(self):
        for site, composition_volumes in self.site_to_glycan_map.items():
            yield site, composition_volumes
