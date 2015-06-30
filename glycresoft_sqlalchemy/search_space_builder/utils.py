import re
import csv
import itertools

from glycresoft_ms2_classification.structure import sequence, constants
from glycresoft_ms2_classification.structure import stub_glycopeptides

Sequence = sequence.Sequence
StubGlycopeptide = stub_glycopeptides.StubGlycopeptide


def fragments(sequence):
    """Generate characteristic 'high energy' CID fragments for a given glycopeptide sequence

    Parameters
    ----------
    sequence : Sequence

    Returns
    -------
    oxonium_ions : list
    b_ions : list
    y_ions : list
    b_ions_hexnac : list
    y_ions_hexnac : list
    stub_ions : list
    """
    fragments = zip(*map(sequence.break_at, range(1, len(sequence))))
    b_type = fragments[0]
    b_ions = []
    b_ions_hexnac = []
    for b in b_type:
        for fm in b:
            key = fm.get_fragment_name()
            if re.search(r'b1\+', key) and constants.EXCLUDE_B1:
                # b1 Ions aren't actually seen in reality, but are an artefact of the generation process
                # so do not include them in the output
                continue
            mass = fm.get_mass()
            golden_pairs = fm.golden_pairs
            if "HexNAc" in key:
                b_ions_hexnac.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                b_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    y_type = fragments[1]
    y_ions = []
    y_ions_hexnac = []
    for y in y_type:
        for fm in y:
            key = fm.get_fragment_name()
            mass = fm.get_mass()
            golden_pairs = fm.golden_pairs
            if "HexNAc" in key:
                y_ions_hexnac.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})
            else:
                y_ions.append({"key": key, "mass": mass, "golden_pairs": golden_pairs})

    pep_stubs = StubGlycopeptide.from_sequence(sequence)

    stub_ions = pep_stubs.get_stubs()
    oxonium_ions = pep_stubs.get_oxonium_ions()
    return (oxonium_ions, b_ions, y_ions,
            b_ions_hexnac, y_ions_hexnac,
            stub_ions)


def flatten(iterable):
    return tuple(itertools.chain.from_iterable(iterable))


def concatenate_csv(csv_list, merge_file):
    if isinstance(csv_list, basestring):
        csv_list = [csv_list]
    headered = False
    with open(merge_file, 'wb') as merge_handle:
        outwriter = csv.writer(merge_handle)
        for path in csv_list:
            with open(path, 'rb') as handle:
                reader = csv.reader(handle)
                g = iter(reader)
                if headered:
                    g.next()
                else:
                    outwriter.writerow(g.next())
                    headered = True
                outwriter.writerows(g)
