import itertools
import re
import logging
try:
    logger = logging.getLogger("enrich_peptides")
except:
    pass
from glycresoft_sqlalchemy.data_model import InformedPeptide, Protein, PipelineModule, make_transient

try:
    range = xrange
except:
    pass


def flatten(iterable):
    return tuple(itertools.chain.from_iterable(iterable))


def edit_distance(query_seq, target_seq):
    previous = range(len(target_seq) + 1)
    for i, new_pos in enumerate((query_seq)):
        current = [i + 1]
        for j, prev_pos in enumerate((target_seq)):
            insertions = previous[j + 1] + 1
            deletions = current[j] + 1
            substitutions = previous[j] + (not new_pos == prev_pos)
            current.append(min(insertions, deletions, substitutions))
        previous = current
    return previous[-1]


def substring_edit_distance(query_seq, target_seq, max_distance=0):
    # assume target_seq is longer than query_seq
    windows = len(target_seq) - len(query_seq) + max_distance
    window_size = len(query_seq)
    for i in range(windows):
        dist = edit_distance(query_seq, target_seq[i:i + window_size])
        if dist <= max_distance:
            return i, i + window_size, dist
    return False


def fast_exact_match(query_seq, target_seq, **kwargs):
    match = re.search(query_seq, target_seq)
    if match:
        return match.start(), match.end(), 0
    return False


class EnrichDistinctPeptides(PipelineModule):
    def __init__(self, database_path, hypothesis_id, target_proteins=None, max_distance=0):
        self.manager = self.manager_type(database_path)
        self.target_proteins = target_proteins
        self.hypothesis_id = hypothesis_id
        self.max_distance = max_distance

    def stream_distinct_peptides(self, protein):
        session = self.manager.session()
        for i, in session.query(InformedPeptide.id).join(Protein).filter(
                Protein.hypothesis_id == self.hypothesis_id).group_by(
                InformedPeptide.modified_peptide_sequence,
                ~InformedPeptide.modified_peptide_sequence.in_(
                   session.query(InformedPeptide.modified_peptide_sequence).filter(
                    InformedPeptide.protein_id == protein.id))
                ).group_by(InformedPeptide.modified_peptide_sequence):
            yield session.query(InformedPeptide).get(i)

    def run(self):
        session = self.manager.session()
        protein_ids = self.target_proteins
        hypothesis_id = self.hypothesis_id

        if protein_ids is None:
            protein_ids = flatten(session.query(Protein.id).filter(Protein.hypothesis_id == hypothesis_id))
        elif isinstance(protein_ids[0], basestring):
            protein_ids = flatten(session.query(Protein.id).filter(Protein.name == name).first()
                                  for name in protein_ids)
        if self.max_distance == 0:
            decider_fn = fast_exact_match
        else:
            decider_fn = substring_edit_distance

        for protein_id in protein_ids:
            target_protein = session.query(Protein).get(protein_id)
            i = 0
            for peptide in self.stream_distinct_peptides(target_protein):
                match = decider_fn(peptide.base_peptide_sequence, target_protein.protein_sequence)
                if match is not False:
                    start, end, distance = match
                    make_transient(peptide)
                    peptide.id = None
                    peptide.protein_id = target_protein.id
                    peptide.start_position = start
                    peptide.end_position = end
                    session.add(peptide)
                i += 1
                if i % 1000 == 0:
                    logger.info("%d peptides handled for %r", i, target_protein)
                    session.commit()
            session.commit()
            ids = session.query(InformedPeptide.id).filter(
                InformedPeptide.protein_id == Protein.id,
                Protein.hypothesis_id == self.hypothesis_id).group_by(
                InformedPeptide.modified_peptide_sequence,
                InformedPeptide.peptide_modifications,
                InformedPeptide.glycosylation_sites,
                InformedPeptide.protein_id).order_by(InformedPeptide.peptide_score.desc())

            q = session.query(InformedPeptide.id).filter(
                InformedPeptide.protein_id == Protein.id,
                Protein.hypothesis_id == self.hypothesis_id,
                ~InformedPeptide.id.in_(ids.correlate(None)))
            conn = session.connection()
            conn.execute(InformedPeptide.__table__.delete(
                InformedPeptide.__table__.c.id.in_(q.selectable)))
            session.commit()
