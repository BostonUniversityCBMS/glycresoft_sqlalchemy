import re
import logging

from pyteomics import mzid
from sqlalchemy import func
from glycresoft_sqlalchemy.data_model import (
    DatabaseManager, Protein, InformedPeptide, ExactMS1GlycopeptideHypothesis)
from glycresoft_sqlalchemy.structure import sequence, modification, residue
from .enzyme import expasy_rules
from glycresoft_sqlalchemy.utils.database_utils import get_or_create
logger = logging.getLogger("mzid")

Sequence = sequence.Sequence
Residue = residue.Residue
Modification = modification.Modification
ModificationNameResolutionError = modification.ModificationNameResolutionError

MzIdentML = mzid.MzIdentML
_local_name = mzid.xml._local_name
peptide_evidence_ref = re.compile(r"(?P<evidence_id>PEPTIDEEVIDENCE_PEPTIDE_\d+_DBSEQUENCE_)(?P<parent_accession>.+)")

PROTEOMICS_SCORE = ["PEAKS:peptideScore", "mascot:score", "PEAKS:proteinScore"]
WHITELIST_GLYCOSITE_PTMS = ["Deamidation"]


class MultipleProteinMatchesException(Exception):
    def __init__(self, message, evidence_id, db_sequences, key):
        Exception.__init__(self, message)
        self.evidence_id = evidence_id
        self.db_sequences = db_sequences
        self.key = key


def protein_names(mzid_path, pattern=r'.*'):
    pattern = re.compile(pattern)
    parser = Parser(mzid_path, retrieve_refs=True, iterative=False, build_id_cache=True)
    for protein in parser.iterfind(
                "ProteinDetectionHypothesis", retrieve_refs=True, recursive=False, iterative=True):
        name = protein['accession']
        if pattern.match(name):
            yield name


class Parser(MzIdentML):
    def _retrieve_refs(self, info, **kwargs):
        """Retrieves and embeds the data for each attribute in `info` that
        ends in _ref. Removes the id attribute from `info`"""
        multi = None
        for k, v in dict(info).items():
            if k.endswith('_ref'):
                try:
                    info.update(self.get_by_id(v, retrieve_refs=True))
                    del info[k]
                    info.pop('id', None)
                except MultipleProteinMatchesException, e:
                    multi = e
                except:
                    is_multi_db_sequence = peptide_evidence_ref.match(info[k])
                    if is_multi_db_sequence:
                        groups = is_multi_db_sequence.groupdict()
                        evidence_id = groups['evidence_id']
                        db_sequences = groups['parent_accession'].split(':')
                        if len(db_sequences) > 1:
                            multi = MultipleProteinMatchesException(
                                "", evidence_id, db_sequences, k)
                            continue
                    # Fall through
                    logger.debug("%s not found", v)
                    info['skip'] = True
                    info[k] = v
        if multi is not None:
            raise multi

    def _insert_param(self, info, param, **kwargs):
        newinfo = self._handle_param(param, **kwargs)
        if not ('name' in info and 'name' in newinfo):
            info.update(newinfo)
        else:
            if not isinstance(info['name'], list):
                info['name'] = [info['name']]
            info['name'].append(newinfo.pop('name'))

    def _find_immediate_params(self, element, **kwargs):
        return element.xpath('./*[local-name()="{}" or local-name()="{}"]'.format("cvParam", "userParam"))

    def _get_info(self, element, **kwargs):
        """Extract info from element's attributes, possibly recursive.
        <cvParam> and <userParam> elements are treated in a special way."""
        name = _local_name(element)
        schema_info = self.schema_info
        if name in {'cvParam', 'userParam'}:
            return self._handle_param(element)

        info = dict(element.attrib)
        # process subelements
        if kwargs.get('recursive'):
            for child in element.iterchildren():
                cname = _local_name(child)
                if cname in {'cvParam', 'userParam'}:
                    self._insert_param(info, child, **kwargs)
                else:
                    if cname not in schema_info['lists']:
                        info[cname] = self._get_info_smart(child, **kwargs)
                    else:
                        info.setdefault(cname, []).append(
                                self._get_info_smart(child, **kwargs))
        else:
            for param in self._find_immediate_params(element):
                self._insert_param(info, param, **kwargs)

        # process element text
        if element.text and element.text.strip():
            stext = element.text.strip()
            if stext:
                if info:
                    info[name] = stext
                else:
                    return stext

        # convert types
        converters = self._converters
        for k, v in info.items():
            for t, a in converters.items():
                if (_local_name(element), k) in schema_info[t]:
                    info[k] = a(v)
        infos = [info]
        try:
            # resolve refs
            if kwargs.get('retrieve_refs'):
                self._retrieve_refs(info, **kwargs)
        except MultipleProteinMatchesException, e:
            evidence_id = e.evidence_id
            db_sequences = e.db_sequences
            key = e.key
            infos = []
            for name in db_sequences:
                dup = info.copy()
                dup[key] = evidence_id + name
                self._retrieve_refs(dup, **kwargs)
                infos.append(dup)

        # flatten the excessive nesting
        for info in infos:
            for k, v in dict(info).items():
                if k in self._structures_to_flatten:
                    info.update(v)
                    del info[k]

            # another simplification
            for k, v in dict(info).items():
                if isinstance(v, dict) and 'name' in v and len(v) == 1:
                    info[k] = v['name']
        out = []
        for info in infos:
            if len(info) == 2 and 'name' in info and (
                    'value' in info or 'values' in info):
                name = info.pop('name')
                info = {name: info.popitem()[1]}
            out.append(info)
        if len(out) == 1:
            out = out[0]
        return out


def remove_peptide_sequence_alterations(base_sequence, insert_sites, delete_sites):
    """
    Remove all the sequence insertions and deletions in order to reconstruct the
    original peptide sequence.

    Parameters
    ----------
    base_sequence : str
        The peptide sequence string which contains a combination
        of insertion and deletions
    insert_sites : list
        A list of (position, None) pairs indicating the position of
        an amino acid insertion to be removed.
    delete_sites : list
        A list of (position, residue) pairs indicating the position
        and identity of amino acids that were deleted and need to be
        re-inserted.

    Returns
    -------
    str
    """
    sequence_copy = list(base_sequence)

    alteration_sites = insert_sites + delete_sites
    alteration_sites.sort()
    shift = 0
    for position, residue_ in alteration_sites:
        if residue_ is None:
            sequence_copy.pop(position - shift)
            shift += 1
        else:
            sequence_copy.insert(position - shift, residue_)
            shift -= 1
    sequence_copy = ''.join(sequence_copy)
    return sequence_copy


def convert_dict_to_sequence(sequence_dict, session, hypothesis_id, enzyme=None, **kwargs):
    base_sequence = sequence_dict["PeptideSequence"]
    try:
        peptide_sequence = Sequence(sequence_dict["PeptideSequence"])
    except residue.UnknownAminoAcidException:
        return 0

    insert_sites = []
    deleteion_sites = []
    counter = 0

    # Keep a count of the number of variable modifications.
    # It may be desirable to threshold sequences for inclusion
    # based on this value.
    modification_counter = 0
    constant_modifications = kwargs.get("constant_modifications", [])
    modification_translation_table = kwargs.get("modification_translation_table", {})

    glycosite_candidates = (sequence.find_n_glycosylation_sequons(
                peptide_sequence, WHITELIST_GLYCOSITE_PTMS))

    if "SubstitutionModification" in sequence_dict:
        subs = sequence_dict["SubstitutionModification"]
        for sub in subs:
            pos = sub['location'] - 1
            replace = Residue(sub["replacementResidue"])
            peptide_sequence[pos][0] = replace
            modification_counter += 1

    if "Modification" in sequence_dict:
        mods = sequence_dict["Modification"]
        for mod in mods:
            pos = mod["location"] - 1
            try:
                try:
                    _name = mod["name"]
                    modification = Modification(_name)
                except KeyError, e:
                    if "unknown modification" in mod:
                        try:
                            _name = mod['unknown modification']
                            if _name in modification_translation_table:
                                modification = modification_translation_table[_name]()
                            else:
                                modification = Modification(_name)
                        except ModificationNameResolutionError:
                            raise KeyError("Cannot find key %s in %r" % (e, mod))
                    else:
                        raise KeyError("Cannot find key %s in %r" % (e, mod))

                try:
                    rule_text = "%s (%s)" % (_name, mod["residues"][0])
                    if (rule_text not in constant_modifications) and not (
                            pos in glycosite_candidates and modification in WHITELIST_GLYCOSITE_PTMS):
                        modification_counter += 1
                except KeyError:
                    modification_counter += 1

                if pos == -1:
                    peptide_sequence.n_term = modification
                elif pos == len(peptide_sequence):
                    peptide_sequence.c_term = modification
                else:
                    peptide_sequence.add_modification(pos, modification)
            except KeyError:
                if "unknown modification" in mod:
                    mod_description = mod["unknown modification"]
                    insertion = re.search(r"(\S{3})\sinsertion", mod_description)
                    deletion = re.search(r"(\S{3})\sdeletion", mod_description)
                    modification_counter += 1
                    if insertion:
                        insert_sites.append((mod['location'] - 1, None))
                    elif deletion:
                        deleteion_sites.append((mod['location'] - 1, mod['residues'][0]))
                    else:
                        raise
                else:
                    raise

    insert_sites.sort()
    deleteion_sites.sort()

    alteration_sites = insert_sites + deleteion_sites
    alteration_sites.sort()

    evidence_list = sequence_dict["PeptideEvidenceRef"]
    # Flatten the evidence list if it has extra nesting because of alternative
    # mzid parsing
    if isinstance(evidence_list[0], list):
        evidence_list = [x for sub in evidence_list for x in sub]
        # for ev in evidence_list:
        #     print ev['accession']
    score = score_type = None
    for k, v in sequence_dict.items():
        if k in PROTEOMICS_SCORE:
            score_type = k
            score = v
            break

    for evidence in evidence_list:
        if "skip" in evidence:
            continue
        parent_protein = session.query(Protein).filter(
            Protein.name == evidence['accession'],
            Protein.hypothesis_id == hypothesis_id).first()
        if parent_protein is None:
            continue
        start = evidence["start"] - 1
        end = evidence["end"]

        sequence_copy = remove_peptide_sequence_alterations(
            base_sequence, insert_sites, deleteion_sites)

        # sequence_copy = list(base_sequence)
        # for i, position in enumerate(insert_sites):
        #     sequence_copy.pop(position - i)
        # sequence_copy = ''.join(sequence_copy)
        found = parent_protein.protein_sequence.find(sequence_copy)
        if found == -1:
            raise ValueError("Peptide not found in Protein\n%s\n%s\n\n%s" % (
                parent_protein.name, parent_protein.protein_sequence, (
                    sequence_copy, sequence_dict
                    )))
        if found != start:
            start = found
            end = start + len(base_sequence)
        try:
            if enzyme is not None:
                missed_cleavages = len(enzyme.findall(base_sequence))
            else:
                missed_cleavages = None
            match = InformedPeptide(
                calculated_mass=peptide_sequence.mass,
                base_peptide_sequence=base_sequence,
                modified_peptide_sequence=str(peptide_sequence),
                count_glycosylation_sites=None,
                count_missed_cleavages=missed_cleavages,
                count_variable_modifications=modification_counter,
                start_position=start,
                end_position=end,
                peptide_score=score,
                peptide_score_type=score_type,
                sequence_length=end - start,
                protein_id=parent_protein.id,
                glycosylation_sites=None,
                other={k: v for k, v in sequence_dict.items() if k not in
                       exclude_keys_from_sequence_dict})
            match.protein = parent_protein

            glycosites = set(match.n_glycan_sequon_sites) | set(sequence.find_n_glycosylation_sequons(
                peptide_sequence, WHITELIST_GLYCOSITE_PTMS))
            match.count_glycosylation_sites = len(glycosites)
            match.glycosylation_sites = list(glycosites)
            session.add(match)
            counter += 1
        except:
            print(evidence)
            raise
    return counter
exclude_keys_from_sequence_dict = set(("PeptideEvidenceRef",))


class Proteome(object):
    def __init__(self, database_path, mzid_path, hypothesis_id=None, hypothesis_type=ExactMS1GlycopeptideHypothesis):
        self.manager = DatabaseManager(database_path)
        self.manager.initialize()
        self.mzid_path = mzid_path
        self.hypothesis_id = hypothesis_id
        self.parser = Parser(mzid_path, retrieve_refs=True, iterative=False, build_id_cache=True)
        self.hypothesis_type = hypothesis_type
        self.enzymes = []
        self.constant_modifications = []
        self.modification_translation_table = {}

        self._load()

    def _load(self):
        self._load_enzyme()
        self._load_modifications()
        self._load_proteins()
        self._load_spectrum_matches()

    def _load_proteins(self):
        session = self.manager.session()
        hypothesis, created = get_or_create(session, self.hypothesis_type, id=self.hypothesis_id)
        session.add(hypothesis)
        session.commit()
        self.hypothesis_id = hypothesis.id
        for protein in self.parser.iterfind(
                "ProteinDetectionHypothesis", retrieve_refs=True, recursive=False, iterative=True):
            seq = protein.pop('Seq')
            try:
                p = Protein(
                    name=protein.pop('accession'),
                    protein_sequence=seq,
                    glycosylation_sites=sequence.find_n_glycosylation_sequons(seq),
                    other=protein,
                    hypothesis_id=self.hypothesis_id)
                session.add(p)
            except residue.UnknownAminoAcidException:
                continue
        session.commit()

    def _load_spectrum_matches(self):
        session = self.manager.session()
        counter = 0
        last = 0
        i = 0
        try:
            enzyme = re.compile(expasy_rules.get(self.enzymes[0]))
        except Exception, e:
            logger.exception("Enzyme not found.", exc_info=e)
            enzyme = None

        for spectrum_identification in self.parser.iterfind(
                "SpectrumIdentificationItem", retrieve_refs=True, iterative=True):
            counter += convert_dict_to_sequence(
                spectrum_identification, session, self.hypothesis_id, enzyme=enzyme,
                constant_modifications=self.constant_modifications)
            i += 1
            if i % 1000 == 0:
                logger.info("%d spectrum matches processed.", i)
            if (counter - last) > 1000:
                last = counter
                logger.info("%d peptides saved.", counter)
        session.commit()

    def peptides(self):
        session = self.manager.session()
        proteins = session.query(Protein).filter(Protein.hypothesis_id == self.hypothesis_id)
        for protein in proteins:
            for informed in protein.informed_peptides:
                yield informed

    def unique_peptides(self):
        session = self.manager.session()
        query = session.query(InformedPeptide, func.count(InformedPeptide)).filter(
            InformedPeptide.protein_id == Protein.id,
            Protein.hypothesis_id == self.hypothesis_id).group_by(InformedPeptide.modified_peptide_sequence)
        return query

    def _load_enzyme(self):
        self.enzymes = list({e['name'].lower() for e in self.parser.iterfind(
            "EnzymeName", retrieve_refs=True, iterative=True)})

    def _load_modifications(self):
        search_param_modifications = list(self.parser.iterfind(
            "ModificationParams", retrieve_refs=True, iterative=True))
        constant_modifications = []

        for param in search_param_modifications:
            for mod in param['SearchModification']:
                try:
                    name = mod['name']
                except KeyError:
                    name = mod['unknown modification']
                    try:
                        Modification(name)
                    except ModificationNameResolutionError:
                        self.modification_conversion_map[name] = modification.AnonymousModificationRule(
                            name, mod['massDelta'])

                residues = mod['residues']
                if mod.get('fixedMod', False):
                    identifier = "%s (%s)" % (name, ''.join(residues).replace(" ", ""))
                    constant_modifications.append(identifier)
        self.constant_modifications = constant_modifications


def protein_names_taskmain():
    import argparse
    app = argparse.ArgumentParser()
    app.add_argument("mzid_path")
    app.add_argument("-p", "--pattern", default=r'.*', required=False)
    args = app.parse_args()
    values = protein_names(**args.__dict__)
    for value in values:
        print value
