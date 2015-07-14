try:
    import urllib2 as url_parser  # python 2
except ImportError:
    import urllib.parse as url_parser  # python 3
from glycresoft_ms2_classification.utils.xml_utils import ET
import json


def uri_decode(uri):
    return url_parser.unquote(uri)


def uri_decode_list(uri_list=None):
    if uri_list is None:
        return None
    return map(uri_decode, uri_list)


class MSDigestParameters(object):
    '''A port of our quick and dirty Protein Prospector MSDigest parser for just
    modifications'''
    @classmethod
    def parse(cls, xml_path):
        tree = ET.parse(xml_path)
        parser = tree.getroot()

        constant_modifications = []
        variable_modifications = []
        missed_cleavages = []
        enzyme = []

        for parameter_node in parser.iterfind(".//parameters"):
            constant_modifications.extend(uri_decode(mod.text) for mod in parameter_node.iterfind(".//const_mod"))
            variable_modifications.extend(uri_decode(mod.text) for mod in parameter_node.iterfind(".//mod_AA"))
            missed_cleavages.extend(int(node.text) for node in parameter_node.iterfind(".//missed_cleavages"))
            enzyme.extend(node.text for node in parameter_node.iterfind("enzyme"))

        return MSDigestParameters(constant_modifications, variable_modifications, missed_cleavages, enzyme)

    def __init__(self, constant_modifications=None, variable_modifications=None, missed_cleavages=None, enzyme=None):
        self.constant_modifications = constant_modifications
        self.variable_modifications = variable_modifications
        self.missed_cleavages = missed_cleavages
        self.enzyme = enzyme

    def __repr__(self):
        return json.dumps(self.__dict__, indent=2)

if __name__ == '__main__':
    import sys
    print(MSDigestParameters.parse(sys.argv[1]))
