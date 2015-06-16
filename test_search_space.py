from glycresoft_sqlalchemy.search_space_builder import search_space_builder

if __name__ == '__main__':
    digest = search_space_builder.parse_digest("./datafiles/KK-Keratin-type1-prospector.xml")
    s = search_space_builder.TheoreticalSearchSpace("./datafiles/ResultOf20140918_01_isos.csv", site_list="./datafiles/sitelist.txt",
                                                    n_processes=6, **digest.__dict__)
    s.run()
