from sqlalchemy import or_
from sqlalchemy.orm import aliased

from glycresoft_sqlalchemy.data_model import DatabaseManager, TandemScan, Peak
from glypy import MonosaccharideResidue


def ppm_error_expr(expr, value, ppm_tolerance):
    span = value * ppm_tolerance
    return expr.between(value - span, value + span)

oxonium_ions = map(MonosaccharideResidue.from_iupac_lite, ["HexNAc", "Hex", "Fuc", "NeuAc", "NeuGc"])


def expression_filter(tolerance):
    return [ppm_error_expr(Peak.neutral_mass, m.mass(), tolerance) for m in oxonium_ions]


def detect_oxonium_ions(session, tolerance=2e-5, minimum_intensity=1000,
                        select=TandemScan, sample_run_id=None):
    filter_expression = [ppm_error_expr(Peak.neutral_mass, m.mass(), tolerance) for m in oxonium_ions]
    q = session.query(select).join(Peak, TandemScan.id == Peak.scan_id).filter(
        or_(*filter_expression), Peak.intensity >= minimum_intensity).group_by(TandemScan.id)
    if sample_run_id is not None:
        q = q.filter(TandemScan.sample_run_id == sample_run_id)
    return q


def detect_monosaccharide_losses(session, tolerance=2e-5, minimum_intensity=1000,
                                 select=TandemScan, sample_run_id=None):
    OtherPeak = aliased(Peak)
    expr = (Peak.neutral_mass - OtherPeak.neutral_mass)
    filter_expression = [ppm_error_expr(expr, m.mass(), tolerance) for m in oxonium_ions]
    query = session.query(select).join(Peak, TandemScan.id == Peak.scan_id).join(
        OtherPeak, TandemScan.id == OtherPeak.scan_id)
    q = query.filter(or_(*filter_expression), Peak.intensity >= minimum_intensity).group_by(TandemScan.id)
    if sample_run_id is not None:
        q = q.filter(TandemScan.sample_run_id == sample_run_id)
    return q
