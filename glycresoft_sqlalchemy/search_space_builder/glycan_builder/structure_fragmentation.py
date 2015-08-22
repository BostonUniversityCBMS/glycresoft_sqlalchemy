import functools
import multiprocessing
import logging

from glypy.io import glycoct

from glycresoft_sqlalchemy.data_model import (
    PipelineModule, TheoreticalGlycanComposition, TheoreticalGlycanStructure, MS2GlycanHypothesis,
    )


logger = logging.getLogger("glycan_structure_fragmentation")


def fragment_glycan(glycan_structure, database_manager, kind=("B", "Y"), max_cleavages=1):
    structure = glycoct.loads(glycan_structure.glycoct).next()
    fragments = {}
    for f in structure.fragments(kind=kind, max_cleavages=max_cleavages):
        fragments[f.name] = f
    glycan_structure._fragments = fragments
    return glycan_structure
