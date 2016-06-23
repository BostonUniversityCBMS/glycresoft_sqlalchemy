from .uniprot_features import get_features_for as _get_uniprot
from .fasta import ProteinFastaFileParser, Protein

from glycresoft_sqlalchemy.utils.collectiontools import decoratordict

from StringIO import StringIO
import requests

database_handler = decoratordict()

efetch_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


@database_handler("DB source NCBI")
def ncbi_get_proteins(accession_list):
    response = requests.get(
        efetch_url, params={
            "db": "protein", "id": ','.join(accession_list), "rettype": 'fasta'})
    response.raise_for_status()
    return list(ProteinFastaFileParser(StringIO(response.text)))


@database_handler("DB source UniProt")
def uniprot_get_proteins(accession_list):
    out = []
    for acc in accession_list:
        data = _get_uniprot(acc)
        out.append(Protein(protein_sequence=data.sequence, name=data.names[0]))
    return out
