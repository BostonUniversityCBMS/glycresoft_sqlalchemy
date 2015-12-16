from collections import defaultdict

from matplotlib import pyplot as plt
import numpy as np

from .colors import (
    NGlycanCompositionColorizer,
    NGlycanCompositionOrderer,
    GlycanLabelTransformer,
    allset, _null_color_chooser)

from glycresoft_sqlalchemy.data_model import (GlycopeptideMatch, object_session)

from glypy.composition.glycan_composition import FrozenGlycanComposition


class GlycoproteinMicroheterogeneitySummary(object):
    def __init__(self, protein, filter_fn, color_chooser=_null_color_chooser, order_chooser=NGlycanCompositionOrderer):
        self.protein = protein
        self.site_list = protein.glycosylation_sites
        self.order_chooser = order_chooser
        self.color_chooser = color_chooser
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

    def glycan_site_specific_abundance_plot(self, legend=True, alpha=0.65):
        axes = {}
        bar_width = 0.3
        for site, abundances in self.site_to_glycan_map.items():
            if len(abundances) == 0:
                continue
            compositions = sorted(abundances, cmp=self.order_chooser)
            volumes = [abundances[c] for c in compositions]
            colors = map(self.color_chooser, compositions)
            xtick_labeler = GlycanLabelTransformer(compositions, self.order_chooser)

            fig, ax = plt.subplots(1)
            indices = np.arange(len(compositions))

            try:
                include_classes = set(map(self.color_chooser.classify, compositions))
            except:
                include_classes = allset()

            ax.bar(indices + bar_width, volumes, bar_width, alpha=alpha, color=colors)
            ax.set_ylabel("Relative Intensity")
            ax.set_xticks(indices + bar_width * 1.5)
            ax.set_xlabel(xtick_labeler.label_key)
            ax.set_xticklabels(tuple(xtick_labeler), rotation=90, ha='center')
            ax.xaxis.set_ticks_position('none')

            ax.set_title("%s at %d" % (self.protein.name, site))
            if legend:
                handles = self.color_chooser.make_legend(
                    include_classes, alpha=alpha)
                if handles:
                    ax.legend(
                        handles=handles,
                        bbox_to_anchor=(1.20, 1.0))
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
            ax.bar(indices + bar_width, volumes, bar_width, alpha=0.5)
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


class GlycoformAbundancePlot(object):
    def __init__(self, query, color_chooser=_null_color_chooser, order_chooser=NGlycanCompositionOrderer):
        self.query = query
        self.abundance_map = defaultdict(float)
        self.color_chooser = color_chooser
        self.order_chooser = order_chooser
        self._total_abundances()

    def _total_abundances(self):
        for item in self.query:
            composition = item.glycan_composition
            composition = FrozenGlycanComposition(composition).serialize()
            abundance = item.total_volume
            self.abundance_map[composition] += abundance

    def abundance_bar_plot(self, legend=True, alpha=0.65, bar_width=0.5, log=False):
        keys = sorted(self.abundance_map, cmp=self.order_chooser)
        indices = np.arange(len(keys)) * 2
        colors = map(self.color_chooser, keys)
        volumes = [self.abundance_map[k] for k in keys]

        try:
            include_classes = set(map(self.color_chooser.classify, keys))
        except:
            include_classes = allset()

        if log:
            volumes = np.log(volumes)

        xtick_labeler = GlycanLabelTransformer(keys, self.order_chooser)

        fig, ax = plt.subplots(1)
        bars = ax.bar(indices + bar_width, volumes, bar_width, color=colors, alpha=alpha, lw=0)
        if log:
            ax.set_ylabel(r"$\log(Relative Intensity)$")
        else:
            ax.set_ylabel("Relative Intensity")
        ax.set_xticks(indices + bar_width * 1.5)

        font_size = max((200. / (len(indices) / 2.)), 3)

        ax.set_xlabel(xtick_labeler.label_key)
        ax.set_xticklabels(tuple(xtick_labeler), rotation=90, ha='center', size=font_size)
        ax.xaxis.set_ticks_position('none')
        ax.set_title("Glycan Composition Total Abundances")
        if legend:
            handles = self.color_chooser.make_legend(
                include_classes, alpha=alpha)
            if handles:
                ax.legend(
                    handles=handles,
                    bbox_to_anchor=(1.20, 1.0))
        return ax
