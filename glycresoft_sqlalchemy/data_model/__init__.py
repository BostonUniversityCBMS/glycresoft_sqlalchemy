from sqlalchemy import func
from sqlalchemy import Numeric
from sqlalchemy.orm.session import object_session
from sqlalchemy.orm import make_transient

Numeric.asdecimal = False

from .base import Hierarchy, Namespace, Base, slurp

from .pipeline_module import PipelineModule, PipelineException

from .generic import (
    MutableList, MutableDict, Taxon, HasTaxonomy,
    ReferenceDatabase, ReferenceAccessionNumber, HasReferenceAccessionNumber,
    _TemplateNumberStore)

from .connection import (
    DatabaseManager, ConnectionManager, SQLiteConnectionManager, session,
    instance_from_uri)

from .hypothesis import (
    Hypothesis, MS1GlycanHypothesis, MS1GlycopeptideHypothesis,
    ExactMS1GlycopeptideHypothesis, MS2GlycanHypothesis, MS2GlycopeptideHypothesis,
    ExactMS2GlycopeptideHypothesis)

from .hypothesis_sample_match import (
    HypothesisSampleMatch, MS1GlycanHypothesisSampleMatch, MS1GlycopeptideHypothesisSampleMatch,
    ExactMS1GlycopeptideHypothesisSampleMatch, MS2GlycanHypothesisSampleMatch,
    ExactMS2GlycopeptideHypothesis, HypothesisSampleMatchToSample, hypothesis_sample_match_type,
    hypothesis_sample_match_root, MS2GlycopeptideHypothesisSampleMatch)


from .sequence_model.peptide import(
    PeptideBase, GlycopeptideBase, NaivePeptide, InformedPeptide,
    TheoreticalGlycopeptideComposition, InformedTheoreticalGlycopeptideComposition,
    TheoreticalGlycopeptide, HasMS1Information, HasTheoreticalFragments, Protein)


from .sequence_model.fragment import(
    TheoreticalPeptideProductIon, TheoreticalGlycopeptideStubIon,
    HasTheoreticalGlycopeptideProductIons)


from .search_result_model.sequence_identification import (
    GlycopeptideMatch, GlycopeptideSpectrumMatch, GlycopeptideSpectrumMatchScore,
    GlycanStructureMatch, GlycanSpectrumMatch, SpectrumMatchBase)

from .search_result_model.peak_grouping import (
    PeakGroupMatchBase, PeakGroupMatchType, PeakGroupMatch, JointPeakGroupMatch,
    TempPeakGroupMatch, PeakGroupMatchToJointPeakGroupMatch, TheoreticalCompositionMap)

from .observed_ions import (
    SampleRun, BUPIDDeconvolutedLCMSMSSampleRun, Decon2LSLCMSSampleRun,
    ScanBase, MSScan, TandemScan, Peak, Decon2LSPeak, Decon2LSPeakGroup,
    Decon2LSPeakToPeakGroupMap, PeakGroupDatabase, MSMSSqlDB, Decon2LSPeakToPeakGroupMap,
    HasPeakChromatogramData)

from .glycomics import (
    GlycanBase, with_glycan_composition, has_glycan_composition, has_glycan_composition_listener,
    TheoreticalGlycanComposition, TheoreticalGlycanCombination, TheoreticalGlycanStructure,
    StructureMotif, TheoreticalGlycanStructureToMotifTable, TheoreticalGlycanCompositionToMotifTable,
    TheoreticalGlycanCombinationTheoreticalGlycanComposition, MassShift, glycoct_parser, FrozenGlycanComposition)

from .json_type import (
    tryjson, clean_dict, new_alchemy_encoder, JSONType)

from .sequence_model.sequencing import (
    SequenceBuildingBlock, SequenceSegment, AminoAcidComposition)

from .scoring_model.peak_group_scoring import PeakGroupScoringModel

from glycresoft_sqlalchemy.utils.database_utils import get_or_create

# Ensure that the PeptideBase relationship is manifest.
# The attribute on Protein only seems to appear after the
# relationship Expression is evaluated at least once. This
# appears to be due to the use of @declared_attr in the PeptideBase
# Mixin.
assert str(GlycopeptideMatch.protein.expression)
assert hasattr(Protein, "glycopeptide_matches")
