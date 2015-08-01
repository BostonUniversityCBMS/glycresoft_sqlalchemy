import logging
logging.basicConfig(level=logging.DEBUG, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")
from glycresoft_sqlalchemy.search_space_builder import pooling_search_space_builder


def test_main():
    s = pooling_search_space_builder.PoolingTheoreticalSearchSpaceBuilder.from_hypothesis(
        "datafiles/naive_glycopeptide.db", 1, 4)
    s.manager.connection_manager.connect_args['timeout'] = 10
    s.start()

if __name__ == '__main__':
    test_main()
