from glycresoft_sqlalchemy.scoring import target_decoy


def test_main():
    target_decoy.TargetDecoyAnalyzer("./datafiles/integrated_omics_simple.db", 1, 2).start()

if __name__ == '__main__':
    test_main()
