from glycresoft_sqlalchemy.data_model import DatabaseManager, AminoAcidComposition
from glycresoft_sqlalchemy.structure import residue, sequence_composition

from glycresoft_sqlalchemy.utils.database_utils import temp_table
from sqlalchemy import select


def generate_all_compositions(session, blocks=None, size=5):
    if blocks is None:
        blocks = sequence_composition.AminoAcidSequenceBuildingBlock.get_all_common_residues()
    T_AminoAcidComposition = AminoAcidComposition.__table__
    accumulator = []
    for composition in sequence_composition.all_compositions(blocks, size):
        d = {
            "composition": str(composition),
            "mass": composition.mass,
            "size": len(composition),
            "count_n_glycosylation": composition.maybe_n_glycosylation()
        }
        accumulator.append(d)
        if len(d) > 1000:
            session.execute(T_AminoAcidComposition.insert(), accumulator)
            accumulator = []
            print "Inserting!"

    session.execute(T_AminoAcidComposition.insert(), accumulator)
    accumulator = []
    return session
