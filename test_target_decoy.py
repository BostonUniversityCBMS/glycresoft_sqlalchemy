from glycresoft_sqlalchemy.scoring import target_decoy

import sys

if __name__ == "__main__":
    target_decoy.TargetDecoyAnalyzer(*sys.argv[1:]).run()
