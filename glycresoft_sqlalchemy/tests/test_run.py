import logging
from glycresoft_sqlalchemy.search_space_builder import pooling_search_space_builder, pooling_make_decoys
from glycresoft_sqlalchemy.matching import matching
from glycresoft_sqlalchemy.scoring import target_decoy

logging.basicConfig(level=logging.DEBUG, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")


def test_main(ms1_results_path, digest_path, site_list_path, observed_ions_path,
         observed_ions_type='bupid_yaml', ms1_tolerance=10e-6, ms2_tolerance=20e-6,
         output_path=None, decoy_type=0):
    digest = pooling_search_space_builder.parse_digest(digest_path)
    builder = pooling_search_space_builder.PoolingTheoreticalSearchSpaceBuilder(
        ms1_results_path, output_path, site_list=site_list_path, n_processes=6,
        **digest.__dict__)
    builder.start()
    builder.session.commit()

    decoy_maker = pooling_make_decoys.PoolingDecoySearchSpaceBuilder(
        builder.db_file_name, n_processes=6, decoy_type=decoy_type)
    decoy_maker.start()
    decoy_maker.session.commit()

    matcher = matching.IonMatching(builder.db_file_name, 1, observed_ions_path,
                                   observed_ions_type, None, ms1_tolerance,
                                   ms2_tolerance, n_processes=8)
    matcher.start()
    matcher = matching.IonMatching(builder.db_file_name, 2, observed_ions_path,
                                   observed_ions_type, None, ms1_tolerance,
                                   ms2_tolerance, n_processes=8)
    matcher.start()

    matcher.session.commit()
    tda = target_decoy.TargetDecoyAnalyzer(builder.db_file_name, 1, 2)
    tda.start()


if __name__ == '__main__':
    test_main("datafiles/Pure_ResultsOf.csv", "datafiles/KK-Keratin-type1-prospector.xml",
         "datafiles/sitelist.txt", "datafiles/20140918_01.db",
         output_path="datafiles/Pure_ResultsOf.preserve_sequons_reverse.db",
         decoy_type=0)
