import datetime
from functools import partial
from glycresoft_sqlalchemy.data_model import MS1GlycanHypothesis, Hypothesis, TheoreticalGlycanComposition, PipelineModule, Hierarchy
from glycresoft_sqlalchemy.search_space_builder.glycan_builder import constrained_combinatorics, glycomedb_utils, registry

from glypy.composition import glycan_composition, composition_transform
from glypy.structure.monosaccharide import ReducedEnd


def composition_from_text_file(path, **kwargs):
    with open(path) as f:
        for line in f:
            line = line.replace("\n", "").replace("\r", "")
            yield glycan_composition.parse(line), []  # No Motif Information


@registry.composition_source_type.register("txt")
class TextGlycanCompositionHypothesisBuilder(PipelineModule):
    HypothesisType = MS1GlycanHypothesis

    def __init__(self, database_path, text_file_path, reduction=None, derivatization=None,
                 hypothesis_id=None, **kwargs):
        self.manager = self.manager_type(database_path)
        self.text_file_path = text_file_path
        self.reduction = reduction
        self.derivatization = derivatization
        self.options = kwargs
        self.hypothesis_id = hypothesis_id
        self.composition_source = text_file_path

    def run(self):
        session = self.manager.session()
        hypothesis = None
        if self.hypothesis_id is None:
            hypothesis = self.HypothesisType(name=self.options.get(
                "hypothesis_name",
                "glycan-hypothesis-%s" % datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")))

            session.add(hypothesis)
            session.commit()

            hypothesis_id = hypothesis.id

        else:
            hypothesis = session.query(Hypothesis).get(self.hypothesis_id)
            hypothesis_id = hypothesis.id

        if self.derivatization is not None:
            derivatization = partial(composition_transform, substituent=self.derivatization)
        else:
            def derivatization(composition):
                pass

        if self.reduction is not None:
            def reduction(composition):
                composition.reducing_end = reduction_type[self.reduction]()
        else:
            def reduction(composition):
                pass

        i = 0
        monosaccharide_identities = set()
        for composition, motifs in composition_from_text_file(self.composition_source, **self.options):
            reduction(composition)
            derivatization(composition)

            theoretical_composition = TheoreticalGlycanComposition(
                composition=composition.serialize(),
                calculated_mass=composition.mass(),
                derivatization=self.derivatization,
                hypothesis_id=hypothesis_id
            )
            # include motifs here
            monosaccharide_identities |= set(composition)

            session.add(theoretical_composition)
            i += 1
            if i % 1000 == 0:
                session.commit()
        monosaccharide_identities = map(str, monosaccharide_identities)
        hypothesis.parameters["monosaccharide_identities"] = monosaccharide_identities
        session.add(hypothesis)
        session.commit()     



class GlycanCompositionHypothesisBuilder(PipelineModule):
    HypothesisType = MS1GlycanHypothesis

    def __init__(self, database_path, composition_source, composition_source_type='txt',
                 reduction=None, derivatization=None, hypothesis_id=None, **kwargs):
        self.manager = self.manager_type(database_path)
        self.composition_source = composition_source
        self.composition_source_type = composition_source_type
        self.options = kwargs
        self.reduction = reduction
        self.derivatization = derivatization
        self.hypothesis_id = hypothesis_id

    def run(self):
        task_fn = registry.composition_source_type[self.composition_source_type]
        task = task_fn(self.database_path, self.composition_source, self.reduction, self.derivatization, hypothesis_id=self.hypothesis_id, **self.options)
        task.start()

    def old_run(self):
        loader = composition_source_type_map[self.composition_source_type]
        session = self.manager.session()
        hypothesis = None
        if self.hypothesis_id is None:
            hypothesis = self.HypothesisType(name=self.options.get(
                "hypothesis_name",
                "glycan-hypothesis-%s" % datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d-%H%M%S")))

            session.add(hypothesis)
            session.commit()

            hypothesis_id = hypothesis.id

        else:
            hypothesis = session.query(Hypothesis).get(self.hypothesis_id)
            hypothesis_id = hypothesis.id

        if self.derivatization is not None:
            derivatization = partial(composition_transform, substituent=self.derivatization)
        else:
            def derivatization(composition):
                pass

        if self.reduction is not None:
            def reduction(composition):
                composition.reducing_end = reduction_map[self.reduction]()
        else:
            def reduction(composition):
                pass

        i = 0
        monosaccharide_identities = set()
        for composition, motifs in loader(self.composition_source, **self.options):
            reduction(composition)
            derivatization(composition)

            theoretical_composition = TheoreticalGlycanComposition(
                composition=composition.serialize(),
                calculated_mass=composition.mass(),
                derivatization=self.derivatization,
                hypothesis_id=hypothesis_id
            )
            # include motifs here
            monosaccharide_identities |= set(composition)

            session.add(theoretical_composition)
            i += 1
            if i % 1000 == 0:
                session.commit()
        monosaccharide_identities = map(str, monosaccharide_identities)
        hypothesis.parameters["monosaccharide_identities"] = monosaccharide_identities
        session.add(hypothesis)
        session.commit()
