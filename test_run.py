import logging
from glycresoft_sqlalchemy.search_space_builder import pooling_search_space_builder, make_decoys
from glycresoft_sqlalchemy.matching import matching
import summarize
logging.basicConfig(level=logging.DEBUG, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")


def main(ms1_results_path, digest_path, site_list_path, observed_ions_path,
         observed_ions_type='bupid_yaml', ms1_tolerance=1e-6, ms2_tolerance=2e-6):
    digest = pooling_search_space_builder.parse_digest(digest_path)
    builder = pooling_search_space_builder.TheoreticalSearchSpace(
        ms1_results_path, site_list=site_list_path, n_processes=6,
        **digest.__dict__)
    builder.run()
    builder.session.close()

    decoy_maker = make_decoys.DecoySearchSpaceBuilder(builder.db_file_name, n_processes=8)
    decoy_maker.run()
    decoy_maker.session.close()

    matcher = matching.IonMatching(builder.db_file_name, observed_ions_path,
                                   observed_ions_type, None, ms1_tolerance,
                                   ms2_tolerance, 8)
    matcher.run()
    matcher.session.close()

    summarize.main(builder.db_file_name)


if __name__ == '__main__':
    main("datafiles/Pure_ResultsOf.csv", "datafiles/KK-Keratin-type1-prospector.xml",
         "datafiles/sitelist.txt", "datafiles/20140918_01.yaml.db")
