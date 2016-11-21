import logging
try:
    logger = logging.getLogger("peak_grouping")
    logging.basicConfig(level='DEBUG')
except Exception, e:
    logging.exception("Logger could not be initialized", exc_info=e)
    raise e

from sqlalchemy.ext.baked import bakery

from glycresoft_sqlalchemy.data_model import (
    Decon2LSPeakGroup,
    PeakGroupMatch, TempPeakGroupMatch, JointPeakGroupMatch,
    func)

TDecon2LSPeakGroup = Decon2LSPeakGroup.__table__
T_TempPeakGroupMatch = TempPeakGroupMatch.__table__
TPeakGroupMatch = PeakGroupMatch.__table__
T_JointPeakGroupMatch = JointPeakGroupMatch.__table__


query_oven = bakery()


def ppm_error(x, y):
    return (x - y) / y


def expanding_window(series, threshold_gap_size=100):
    windows = []
    current_window = []
    last_item = 0
    for item in series:
        if item - last_item > threshold_gap_size:
            if len(current_window) > 0:
                windows.append(current_window)
            current_window = []
        current_window.append(item)
        last_item = item
    if len(current_window) > 0:
        windows.append(current_window)
    return windows


def centroid_scan_error_regression(session, minimum_abundance=250,
                                   source_model=Decon2LSPeakGroup, filter_fn=lambda x: x):
    '''
    alpha, beta = centroid_scan_error_regression(session)
    r = session.execute(select([func.abs(
            source_model.centroid_scan_estimate - (
                alpha + source_model.weighted_monoisotopic_mass * beta))]))
    r.fetchmany(200)
    '''
    mean_centroid_scan_estimate = float(filter_fn(session.query(func.avg(source_model.centroid_scan_estimate)).filter(
        source_model.total_volume > minimum_abundance)).first()[0])
    mean_weighted_mass = float(filter_fn(session.query(func.avg(source_model.weighted_monoisotopic_mass)).filter(
        source_model.total_volume > minimum_abundance)).first()[0])
    beta = float(
        filter_fn(session.query(
            (func.sum(
                (source_model.weighted_monoisotopic_mass - mean_weighted_mass) *
                source_model.centroid_scan_estimate)) / func.sum(
                ((source_model.weighted_monoisotopic_mass - mean_weighted_mass) *
                 (source_model.weighted_monoisotopic_mass - mean_weighted_mass))
            )).filter(source_model.total_volume > minimum_abundance)).first()[0])
    alpha = mean_centroid_scan_estimate - beta * mean_weighted_mass
    return alpha, beta


def expected_a_peak_regression(session, source_model=Decon2LSPeakGroup, filter_fn=lambda x: x):
    '''
    alpha, beta = expected_a_peak_regression(session)
    r = session.execute(select([func.abs(
            source_model.average_a_to_a_plus_2_ratio - (
                alpha + source_model.weighted_monoisotopic_mass * beta))]))
    r.fetchmany(200)
    '''
    mean_weighted_mass = float(filter_fn(session.query(func.avg(
        source_model.weighted_monoisotopic_mass))).first()[0])
    mean_a_to_a_plus_2_ratio = float(filter_fn(session.query(func.avg(
        source_model.average_a_to_a_plus_2_ratio))).first()[0])
    beta = float(
        filter_fn(
            session.query((func.sum(
                (source_model.weighted_monoisotopic_mass - mean_weighted_mass) *
                source_model.average_a_to_a_plus_2_ratio)) / func.sum(
                (source_model.weighted_monoisotopic_mass - mean_weighted_mass) *
                (source_model.weighted_monoisotopic_mass - mean_weighted_mass)))).first()[0])
    alpha = mean_a_to_a_plus_2_ratio - beta * mean_weighted_mass
    return alpha, beta
