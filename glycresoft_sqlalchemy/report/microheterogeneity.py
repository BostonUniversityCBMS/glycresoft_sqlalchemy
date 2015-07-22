from glycresoft_sqlalchemy.data_model import Hypothesis, Protein, GlycopeptideMatch, Glycan


class GlycoproteinMicroheterogeneitySummary(object):
    def __init__(self, protein, score_threshold=0.2):
        self.protein = protein
        self.site_list = protein.glycosylation_sites
        self.site_to_glycan_map = {}
        self.score_threshold = score_threshold

    def build_site_to_glycan_map(self):
        site_map = {}
        for site in self.site_list:
            species = set()
            for match in self.protein.glycopeptide_matches.filter(
                    GlycopeptideMatch.spans(site),
                    GlycopeptideMatch.ms2_score > self.score_threshold):
                species.add(match.glycan_composition_str)

            site_map[site] = species
        self.site_to_glycan_map = site_map
        return site_map
