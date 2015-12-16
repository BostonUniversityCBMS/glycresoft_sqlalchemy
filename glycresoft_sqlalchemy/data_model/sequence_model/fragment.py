from collections import OrderedDict

from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.orm import relationship, backref, make_transient, Query, aliased
from sqlalchemy.ext.hybrid import hybrid_property

from sqlalchemy import (PickleType, Numeric, Unicode, and_,
                        Column, Integer, ForeignKey, UnicodeText, Boolean)

from ..base import Base2 as Base
from ...structure import fragment


class IonDictFacade(object):

    def __getitem__(self, key):
        if key == 'key':
            return self.name
        elif key == 'mass':
            return self.calculated_mass
        else:
            raise KeyError(key)


# OWNER = "TheoreticalGlycopeptide"
OWNER = "Example"


class TheoreticalPeptideProductIon(IonDictFacade, Base):
    __tablename__ = "TheoreticalPeptideProductIon"

    id = Column(Integer, primary_key=True)
    name = Column(Unicode(64), index=True)
    calculated_mass = Column(Numeric(12, 6, asdecimal=False), index=True)

    n_term = Column(Unicode(10))
    c_term = Column(Unicode(10))
    glycosylated = Column(Boolean, index=True)
    sequence_index = Column(Integer, index=True)
    ion_series = Column(Unicode(10), index=True)

    theoretical_glycopeptide_id = Column(Integer, ForeignKey("Example.id"), index=True)
    theoretical_glycopeptide = relationship(OWNER)

    @classmethod
    def from_fragment(cls, frag):
        inst = cls(
            name=frag.name,
            calculated_mass=frag.mass,
            n_term=frag.flanking_amino_acids[0].name,
            c_term=frag.flanking_amino_acids[1].name,
            glycosylated=('HexNAc' in frag.mod_dict),
            sequence_index=frag.position,
            ion_series=frag.name[0]
            )
        return inst

    def __repr__(self):
        return "<TheoreticalPeptideProductIon %s %s-%s %0.3f>" % (
            self.name, self.n_term, self.c_term, self.calculated_mass)


class TheoreticalGlycopeptideStubIon(IonDictFacade, Base):
    __tablename__ = "TheoreticalGlycopeptideStubIon"

    id = Column(Integer, primary_key=True)
    name = Column(Unicode(64), index=True)
    calculated_mass = Column(Numeric(12, 6, asdecimal=False), index=True)
    ion_series = Column(Unicode(10), index=True)
    includes_peptide = Column(Boolean)

    theoretical_glycopeptide_id = Column(Integer, ForeignKey("Example.id"), index=True)
    theoretical_glycopeptide = relationship(OWNER)

    def __repr__(self):
        return "<TheoreticalGlycopeptideStubIon %s %s %0.3f>" % (
            self.name, self.ion_series, self.calculated_mass)

    @classmethod
    def from_simple_fragment(cls, frag):
        inst = cls(
            name=frag.name,
            calculated_mass=frag.mass,
            ion_series=frag.kind,
            includes_peptide=fragment.IonSeries(frag.kind).includes_peptide
            )
        return inst


class HasTheoreticalGlycopeptideProductIons(object):

    @declared_attr
    def theoretical_glycopeptide_stub_ions(cls):
        return relationship(TheoreticalGlycopeptideStubIon, lazy="dynamic")

    @declared_attr
    def theoretical_peptide_product_ions(cls):
        return relationship(TheoreticalPeptideProductIon, lazy="dynamic")

    # Specialized Series Queries For Stub Glycopeptides / Oxonium Ions
    @hybrid_property
    def oxonium_ions(self):
        return self.theoretical_glycopeptide_stub_ions.filter_by(ion_series="oxonium_ion")

    @oxonium_ions.expression
    def oxonium_ions(cls):
        return and_((cls.id == TheoreticalGlycopeptideStubIon.theoretical_glycopeptide_id),
                    (TheoreticalGlycopeptideStubIon.ion_series == "oxonium_ion"))

    @hybrid_property
    def stub_glycopeptides(self):
        return self.theoretical_glycopeptide_stub_ions.filter_by(ion_series="stub_glycopeptide")

    @stub_glycopeptides.expression
    def stub_glycopeptides(cls):
        return and_((cls.id == TheoreticalGlycopeptideStubIon.theoretical_glycopeptide_id),
                    (TheoreticalGlycopeptideStubIon.ion_series == "stub_glycopeptide"))

    # Specialized Series Queries For b/y ions
    @hybrid_property
    def bare_b_ions(self):
        return self.theoretical_peptide_product_ions.filter_by(ion_series="b", glycosylated=False)

    @bare_b_ions.expression
    def bare_b_ions(cls):
        return and_((cls.id == TheoreticalPeptideProductIon.theoretical_glycopeptide_id),
                    (TheoreticalPeptideProductIon.ion_series == "b"),
                    (~TheoreticalPeptideProductIon.glycosylated))

    @hybrid_property
    def bare_y_ions(self):
        return self.theoretical_peptide_product_ions.filter_by(ion_series="y", glycosylated=False)

    @bare_b_ions.expression
    def bare_b_ions(cls):
        return and_((cls.id == TheoreticalPeptideProductIon.theoretical_glycopeptide_id),
                    (TheoreticalPeptideProductIon.ion_series == "y"),
                    (~TheoreticalPeptideProductIon.glycosylated))

    @hybrid_property
    def glycosylated_b_ions(self):
        return self.theoretical_peptide_product_ions.filter_by(ion_series="b", glycosylated=True)

    @glycosylated_b_ions.expression
    def glycosylated_b_ions(cls):
        return and_((cls.id == TheoreticalPeptideProductIon.theoretical_glycopeptide_id),
                    (TheoreticalPeptideProductIon.ion_series == "b"),
                    (TheoreticalPeptideProductIon.glycosylated))

    @hybrid_property
    def glycosylated_y_ions(self):
        return self.theoretical_peptide_product_ions.filter_by(ion_series="b", glycosylated=True)

    @glycosylated_y_ions.expression
    def glycosylated_y_ions(cls):
        return and_((cls.id == TheoreticalPeptideProductIon.theoretical_glycopeptide_id),
                    (TheoreticalPeptideProductIon.ion_series == "y"),
                    (TheoreticalPeptideProductIon.glycosylated))


class Example(HasTheoreticalGlycopeptideProductIons, Base):
    __tablename__ = "Example"
    id = Column(Integer, primary_key=True)
