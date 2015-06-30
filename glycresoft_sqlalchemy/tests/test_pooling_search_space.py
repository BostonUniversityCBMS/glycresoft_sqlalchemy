import logging
logging.basicConfig(level=logging.DEBUG, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")

from glycresoft_sqlalchemy.search_space_builder import pooling_search_space_builder


def test_main():
    digest = pooling_search_space_builder.parse_digest("./datafiles/KK-Keratin-type1-prospector.xml")
    s = pooling_search_space_builder.PoolingTheoreticalSearchSpace("./datafiles/ResultOf20140918_01_isos.csv", site_list="./datafiles/sitelist.txt",
                                                    n_processes=6, **digest.__dict__)
    s.start()


if __name__ == '__main__':
    test_main()
