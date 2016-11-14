# Deprecated

This program is deprecated in favor of https://github.com/mobiusklein/glycresoft_desktop, which wraps the application server in https://github.com/mobiusklein/glycresoft_app to run https://github.com/mobiusklein/glycan_profiling

## GlycReSoft
This is a complete re-write of the GlycReSoft glycomics and glycoproteomics pipeline in `Python`. It is suffixed `_sqlalchemy` to reflect the fact that the majority of its data model is implemented using the SQLAlchemy database toolkit.

The pipeline was redesigned to be platform-agnostic.

### Features
 - [Combinatorial Naive Hypotheses](#combinatorial-naive-hypothesis)
 - Data-Driven Hypotheses
 - LC-MS Profiling/Database Search
 - LC-MS/MS Sequence Database Search

### Combinatorial Naive Hypotheses
We implement an algorithm for constructing combinatorial hypotheses from simple composition rules for glycans and combinations for *N-Linked* glycopeptides.

We also provide an interface to glycomics databases such as Glycome-DB and GlyTouCan for constructing large naive hypotheses using the [glypy](https://github.com/mobiusklein/glypy) library.

### Data-Driven Hypotheses
We provide a framework for building refined glycopeptide hypotheses from proteomics and/or glycomics experimental results.

Proteomics results in `mzIdentML 1.1` format provide a definition for modified peptide sequences and proteins.
