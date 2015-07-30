import logging
logging.basicConfig(level=logging.DEBUG, filemode='w',
                    format="%(asctime)s - %(name)s:%(funcName)s:%(lineno)d - %(levelname)s - %(message)s",
                    datefmt="%H:%M:%S")

from glycresoft_sqlalchemy.search_space_builder import search_space_builder


def test_main():
    digest = search_space_builder.parse_digest("./datafiles/KK-Keratin-type1-prospector.xml")
    print digest
    constant_mods, variable_mods = (["Carbamidomethyl (C)"], ["Deamidated (Q)", "Deamidated (N)"])
    s = search_space_builder.TheoreticalSearchSpaceBuilder(
        "./datafiles/naive_glycopeptide.test_on_20140918_01_isos.glycopeptide_compositions.csv",
        site_list="./datafiles/agp-sitelist-long-names.fa", constant_modifications=constant_mods,
        variable_modifications=variable_mods, n_processes=4, enzyme=["trypsin"])
    s.start()

if __name__ == '__main__':
    test_main()
