from glycresoft_sqlalchemy.data_model import AminoAcidComposition
from glycresoft_sqlalchemy.structure import sequence_composition


Composition = sequence_composition.Composition
SequenceComposition = sequence_composition.SequenceComposition


def generate_all_compositions(session, blocks=None, size=5, terminal_composition=Composition()):
    if blocks is None:
        blocks = sequence_composition.AminoAcidSequenceBuildingBlock.get_all_sequencing_residues()
    T_AminoAcidComposition = AminoAcidComposition.__table__
    accumulator = []
    for composition in sequence_composition.all_compositions(blocks, size):
        composition.composition_offset = terminal_composition
        d = {
            "composition": str(composition),
            "mass": composition.mass,
            "size": sum(composition.values()),
            "count_n_glycosylation": composition.maybe_n_glycosylation()
        }
        accumulator.append(d)
        if len(d) > 1000:
            session.execute(T_AminoAcidComposition.insert(), accumulator)
            accumulator = []

    session.execute(T_AminoAcidComposition.insert(), accumulator)
    accumulator = []
    return session
