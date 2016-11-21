from lxml import etree

from glycresoft_sqlalchemy.utils import simple_repr, make_struct
from glycresoft_sqlalchemy.structure import glycan, sequence, residue

uri_template = "http://www.uniprot.org/uniprot/{accession}.xml"


class UniProtFeatureBase(object):
    __repr__ = simple_repr


class SignalPeptide(UniProtFeatureBase):
    feature_type = 'signal peptide'

    def __init__(self, start, end):
        self.start = start
        self.end = end

    @classmethod
    def fromxml(cls, feature):
        return cls(
            int(feature.find(".//{http://uniprot.org/uniprot}begin").attrib['position']),
            int(feature.find(".//{http://uniprot.org/uniprot}end").attrib['position']))


class ModifiedResidue(UniProtFeatureBase):
    feature_type = 'modified residue'

    def __init__(self, position, modification):
        self.position = position
        self.modification = modification

    @classmethod
    def fromxml(cls, feature):
        return cls(
            int(feature.find(".//{http://uniprot.org/uniprot}position").attrib['position']),
            feature.attrib["description"])


def handle_glycosylation_site(feature):
    position = int(feature.find(".//{http://uniprot.org/uniprot}position").attrib['position'])
    glycosylation_type = glycan.GlycosylationType[feature.attrib["description"].split(" ")[0]]
    return glycan.GlycosylationSite(parent=None, position=position, site_types=[glycosylation_type])


UniProtProtein = make_struct("UniProtProtein", ("sequence", "features", "names"))


def get_etree_for(accession):
    tree = etree.parse(uri_template.format(accession=accession))
    return tree


def get_features_for(accession):
    tree = etree.parse(uri_template.format(accession=accession))
    seq = tree.find(".//{http://uniprot.org/uniprot}entry/{http://uniprot.org/uniprot}sequence").text.replace("\n", '')
    names = [el.text for el in tree.findall(
        ".//{http://uniprot.org/uniprot}protein/*/{http://uniprot.org/uniprot}fullName")]
    features = []
    for tag in tree.findall(".//{http://uniprot.org/uniprot}feature"):
        feature_type = tag.attrib['type']
        if feature_type == 'signal peptide':
            features.append(SignalPeptide.fromxml(tag))
        elif feature_type == "modified residue":
            features.append(ModifiedResidue.fromxml(tag))
        elif feature_type == "glycosylation site":
            features.append(handle_glycosylation_site(tag))
    return UniProtProtein(seq, features, names)
