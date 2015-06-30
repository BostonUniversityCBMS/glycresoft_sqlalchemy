from glycresoft_sqlalchemy.data_model import DatabaseManager, Hypothesis, Protein, TheoreticalGlycopeptide, GlycopeptideMatch
from glycresoft_sqlalchemy.search_space_builder import naive_glycopeptide_hypothesis


def taskmain(database_path, hypothesis_name, protein_file, site_list_file,
             glycan_file, glycan_file_type, constant_modifications, variable_modifications,
             enzyme, max_missed_cleavages=1, output_path=None, n_processes=4):
    task = naive_glycopeptide_hypothesis.NaiveGlycopeptideHypothesisBuilder(
        database_path=database_path,
        hypothesis_name=hypothesis_name,
        protein_file=protein_file,
        site_list_file=site_list_file,
        glycan_file=glycan_file,
        glycan_file_type=glycan_file_type,
        constant_modifications=constant_modifications,
        variable_modifications=variable_modifications,
        enzyme=enzyme,
        max_missed_cleavages=max_missed_cleavages,
        n_processes=n_processes
        )
    hypothesis_id = task.start()
    if task.status != 0:
        raise task.error
    task = naive_glycopeptide_hypothesis.NaiveGlycopeptideHypothesiMS1LegacyCSV(
        database_path=database_path,
        hypothesis_id=hypothesis_id,
        output_path=output_path
        )
    task.start()
    if task.status != 0:
        raise task.error
    return hypothesis_id
