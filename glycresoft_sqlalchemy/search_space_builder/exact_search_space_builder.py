import csv
import os
import re
import datetime
import multiprocessing
import logging

import functools

from glycresoft_ms2_classification.structure.modification import ModificationTable
from glycresoft_ms2_classification.structure.sequence import Sequence
from glycresoft_ms2_classification.structure.stub_glycopeptides import StubGlycopeptide
from glycresoft_ms2_classification.structure import constants
from glycresoft_ms2_classification.proteomics import get_enzyme, msdigest_xml_parser

from .search_space_builder import MS1GlycopeptideResult, get_monosaccharide_identities, parse_site_file, MS1ResultsFile
from .. import data_model as model
from ..data_model import PipelineModule

TheoreticalGlycopeptide = model.TheoreticalGlycopeptide
Protein = model.Protein

logger = logging.getLogger("search_space_builder")


def generate_fragments(seq, ms1_result):
    seq_mod = seq.get_sequence()
    fragments = zip(*map(seq.break_at, range(1, len(seq))))
    b_type = fragments[0]
    b_ions = []
    b_ions_hexnac = []
    for b in b_type:
        for fm in b:
            key = fm.get_fragment_name()
            if key == ("b1" or re.search(r'b1\+', key)) and constants.EXCLUDE_B1:
                # B1 Ions aren't actually seen in reality, but are an artefact of the generation process
                # so do not include them in the output
                continue
            mass = fm.get_mass()
            golden_pairs = fm.golden_pairs
            if "HexNAc" in key:
                b_ions_hexnac.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                b_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    y_type = fragments[1]  # seq.get_fragments('Y')
    y_ions = []
    y_ions_hexnac = []
    for y in y_type:
        for fm in y:
            key = fm.get_fragment_name()
            mass = fm.get_mass()
            golden_pairs = fm.golden_pairs
            if "HexNAc" in key:
                y_ions_hexnac.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                y_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    pep_stubs = StubGlycopeptide(
        ms1_result.base_peptide_sequence,
        ms1_result.peptide_modifications,
        ms1_result.count_glycosylation_sites,
        ms1_result.glycan_composition_str)

    stub_ions = pep_stubs.get_stubs()
    oxonium_ions = pep_stubs.get_oxonium_ions()

    theoretical_glycopeptide = TheoreticalGlycopeptide(
        ms1_score=ms1_result.ms1_score,
        observed_mass=ms1_result.observed_mass,
        calculated_mass=ms1_result.calculated_mass,
        ppm_error=ms1_result.ppm_error,
        volume=ms1_result.volume,
        count_glycosylation_sites=ms1_result.count_glycosylation_sites,
        count_missed_cleavages=ms1_result.count_missed_cleavages,
        start_position=ms1_result.start_position,
        end_position=ms1_result.end_position,
        base_peptide_sequence=ms1_result.base_peptide_sequence,
        modified_peptide_sequence=seq_mod,
        peptide_modifications=ms1_result.peptide_modifications,
        glycopeptide_sequence=seq_mod + ms1_result.glycan_composition_str,
        sequence_length=len(seq),
        glycan_composition_str=ms1_result.glycan_composition_str,
        bare_b_ions=b_ions,
        bare_y_ions=y_ions,
        oxonium_ions=oxonium_ions,
        stub_ions=stub_ions,
        glycosylated_b_ions=b_ions_hexnac,
        glycosylated_y_ions=y_ions_hexnac,
        protein_id=ms1_result.protein_id
        )
    return theoretical_glycopeptide


def from_sequence(ms1_result, protein_map):
    if len(ms1_result.base_peptide_sequence) == 0:
        return None
    if ms1_result.modified_peptide_sequence is None:
        ms1_result.modified_peptide_sequence = ms1_result.base_peptide_sequence
    seq = Sequence(ms1_result.modified_peptide_sequence)
    seq.glycan = ''
    product = generate_fragments(seq, ms1_result)
    if not isinstance(product.protein_id, int):
        product.protein_id = protein_map[product.protein_id]
    return product


class ExactSearchSpaceBuilder(PipelineModule):
    manager_type = model.DatabaseManager

    def __init__(self, ms1_results_file, db_file_name, hypothesis_id,
                 enzyme=None,
                 constant_modifications=None,
                 variable_modifications=None, n_processes=4, **kwargs):
        self.db_file_name = db_file_name

        self.hypothesis_id = hypothesis_id
        self.manager = self.manager_type(db_file_name)
        self.session = self.manager.session()

        self.ms1_results_file = ms1_results_file

        tag = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")

        self.hypothesis = self.session.query(model.Hypothesis).get(hypothesis_id)

        self.ms1_results_reader = MS1ResultsFile(self.ms1_results_file)

        self.n_processes = n_processes

        self.monosaccharide_identities = self.ms1_results_reader.monosaccharide_identities
        enzyme = map(get_enzyme, enzyme)

        self.hypothesis.parameters = {
            "monosaccharide_identities": self.monosaccharide_identities,
            "enzyme": enzyme,
            "constant_modification_list": constant_modifications,
            "variable_modification_list": variable_modifications,
            "ms1_output_file": ms1_results_file,
            "enzyme": enzyme,
            "tag": tag,
            "enable_partial_hexnac_match": constants.PARTIAL_HEXNAC_LOSS
        }
        self.session.add(self.hypothesis)

        self.session.commit()

    def prepare_task_fn(self):
        protein_map = dict(self.session.query(Protein.name, Protein.id).filter(Protein.hypothesis_id == self.hypothesis.id))
        return functools.partial(from_sequence, protein_map=protein_map)

    def run(self):
        session = self.session
        cntr = 0
        last = 0
        commit_interval = 1000
        task_fn = self.prepare_task_fn()
        if self.n_processes > 1:
            pool = multiprocessing.Pool(self.n_processes)
            for theoretical in pool.imap_unordered(task_fn, self.ms1_results_reader, chunksize=500):
                if theoretical is None:
                    continue
                session.add(theoretical)
                cntr += 1
                if cntr > last + commit_interval == 0:
                    session.commit()
                    logger.info("%d records handled", cntr)
                    last = cntr

        for ms1_result in self.ms1_results_reader:
            theoretical = task_fn(ms1_result)
            if theoretical is None:
                continue
            session.add(theoretical)
            cntr += 1
            if cntr > last + commit_interval == 0:
                session.commit()
                logger.info("%d records handled", cntr)
                last = cntr
        session.commit()
        session.close()
