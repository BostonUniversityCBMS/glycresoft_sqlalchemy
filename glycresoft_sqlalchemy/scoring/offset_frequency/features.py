from glycresoft_sqlalchemy.data_model import GlycopeptideMatch
import logging
from glycresoft_sqlalchemy.structure import composition
from glycresoft_sqlalchemy.structure import modification
from .peak_relations import (
    MassOffsetFeature, search_features_on_spectrum,
    feature_function_estimator, FittedFeature)
from .utils import chain_iterable
from glypy import MonosaccharideResidue, monosaccharides
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder.ms2 import residue_counter

logger = logging.getLogger("offset_frequency")


def get_link_features(database_path, hypothesis_id, score_threshold=0.7):
    '''
    Given a hypothesis id for a hypothesis with high scoring GlycopeptideMatch hits,
    infer what the base linking features are. These features will include separate entries
    for different modification states.
    '''
    def filter_fn(q):
        return q.filter(GlycopeptideMatch.ms2_score >= score_threshold)

    job = residue_counter.ResidueCounter(database_path, hypothesis_id, filter_fn, n_processes=1)
    job.source = GlycopeptideMatch
    counts = job.start()
    link_features = []
    for building_block in counts:
        f = MassOffsetFeature(building_block.neutral_mass, name=repr(building_block), feature_type='link')
        link_features.append(f)
    return link_features


def fit_features(training_set, features, kinds=('b', 'y')):
    fitted_features = []
    for feat in features:
        logger.info("Fitting %r", feat)
        for kind in kinds:
            logger.info("for %s", kind)
            u, v, rels = feature_function_estimator(training_set, feat, kind)
            fit_feat = FittedFeature(feat, kind, u, v, rels)
            logger.info("fit: %r", fit_feat)
            fitted_features.append(fit_feat)
    return fitted_features


def specialize_features(fitted_features):
    specialized = []
    for feat in fitted_features:
        for charge_pairs_intensity_ratio, count in feat.charge_intensity_ratio().items():
            if count < 3:
                continue
            charge_pairs, intensity_ratio = charge_pairs_intensity_ratio
            from_charge, to_charge = charge_pairs
            offset = feat.feature.offset
            name = feat.feature.name + " %d->%d / %d" % (from_charge, to_charge, intensity_ratio)
            f = MassOffsetFeature(
                offset=offset, name=name, tolerance=2e-5,
                from_charge=from_charge, to_charge=to_charge,
                intensity_ratio=intensity_ratio)
            specialized.append(f)
    return specialized


def peak_probability_distribution(peak_list, features):
    for peak in peak_list:
        satisfied_features = search_features_on_spectrum(peak, peak_list, features)
        distribution = {}
        



shifts = [
    MassOffsetFeature(
        -composition.Composition("NH3").mass, name='Neutral-Loss-Ammonia',
        tolerance=2e-5, feature_type='neutral_loss'),
    MassOffsetFeature(
        -composition.Composition("H2O").mass, name='Neutral-Loss-Water',
        tolerance=2e-5, feature_type='neutral_loss'),
    MassOffsetFeature(
        0.0, name="Charge-Increased Species", from_charge=1,
        to_charge=2, tolerance=2e-5, feature_type='support_relation'),
    MassOffsetFeature(
        0.0, name="Charge-Decreased Species", from_charge=2,
        to_charge=1, tolerance=2e-5, feature_type='support_relation'),
]

for monosacch_name in ['Hex', 'HexNAc', 'Fuc', 'NeuAc']:
    m = MonosaccharideResidue.from_monosaccharide(monosaccharides[monosacch_name])
    shifts.append(
        MassOffsetFeature(
            m.mass(), name=monosacch_name, tolerance=2e-5, feature_type='glycosylation'))
