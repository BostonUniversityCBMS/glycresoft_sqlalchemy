import warnings

# Required information for generating peptides that can be cleaved by a given protease
enzymes = {
    "trypsin": {
        "cleavage_start": [""],
        "cleavage_end": ["K", "R"],
        "name": "trypsin"
    }
}


def get_enzyme(name):
    try:
        return enzymes[name.lower()]
    except:
        warnings.warn("Could not identify protease {}".format(name))
        return name
