from glycresoft_sqlalchemy.data_model import AminoAcidComposition, PipelineModule
from glycresoft_sqlalchemy.search_space_builder.glycopeptide_builder.ms1 import generate_sequence_compositions
from glycresoft_sqlalchemy.structure import sequence_composition, modification, residue

Residue = residue.Residue
Modification = modification.Modification
ModificationBuildingBlock = sequence_composition.ModificationBuildingBlock
AminoAcidSequenceBuildingBlock = sequence_composition.AminoAcidSequenceBuildingBlock
SequenceComposition = sequence_composition.SequenceComposition


class DegenerateSegmentBlock(object):
    def __init__(self, composition):
        self.composition = composition
        self.neutral_mass = composition.mass
        self._string = str(composition)
        self._hash = hash(self.composition.keys()[0])

    def _reset(self):
        self._string = str(self.composition)
        self._hash = hash(self.composition.keys()[0])
        self.neutral_mass = self.composition.mass

    def __eq__(self, other):
        try:
            return self.composition == other.composition
        except:
            return str(self) == str(other)

    def __ne__(self, other):
        try:
            return self.composition != other.composition
        except:
            return str(self) != str(other)

    def clone(self):
        return self.__class__(self.composition.clone())

    def __hash__(self):
        return self._hash

    def __repr__(self):
        return self._string

    def orderings(self):
        for seq in self.composition.orderings():
            yield seq

    def __contains__(self, item):
        return item in self.composition

    def add_modification(self, modification):
        self.composition[ModificationBuildingBlock(modification)] += 1
        self._reset()
        return self


class DegenerateSegmentBlockBuilder(PipelineModule):
    def __init__(self, database_path, building_blocks=None, max_size=5):
        if building_blocks is None:
            building_blocks = AminoAcidSequenceBuildingBlock.get_all_sequencing_residues()
            building_blocks.append(
                AminoAcidSequenceBuildingBlock(
                    Residue(symbol="C"), modifications=[Modification("Carbamidomethyl")]))
            building_blocks.append(
                AminoAcidSequenceBuildingBlock(
                    Residue(symbol="N"), modifications=[Modification("HexNAc")]))

        self.manager = self.manager_type(database_path)
        self.building_blocks = building_blocks
        self.max_size = max_size

    def run(self):
        generate_sequence_compositions.generate_all_compositions()


def get_segments(session, size=2):
    return [
        DegenerateSegmentBlock(a.to_sequence_composition()) for a in session.query(AminoAcidComposition).filter(
            AminoAcidComposition.size == size)
    ]
