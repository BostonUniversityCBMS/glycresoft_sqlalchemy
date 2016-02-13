from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.orm import relationship
from sqlalchemy.ext.hybrid import hybrid_property

from sqlalchemy import (Numeric, Unicode, and_,
                        Column, Integer, ForeignKey, Boolean)

from ..base import Base
from ...structure import fragment


class IonDictFacade(object):

    def __getitem__(self, key):
        if key == 'key':
            return self.name
        elif key == 'mass':
            return self.calculated_mass
        else:
            raise KeyError(key)


class ProductIonBase(IonDictFacade):
    id = Column(Integer, primary_key=True)
    name = Column(Unicode(64), index=True)
    calculated_mass = Column(Numeric(12, 6, asdecimal=False), index=True)
    neutral_loss = Column(Unicode(24))
    ion_series = Column(Unicode(24))

    @declared_attr
    def theoretical_glycopeptide_id(self):
        return Column(Integer, ForeignKey("TheoreticalGlycopeptide.id"), index=True)

    @declared_attr
    def theoretical_glycopeptide(self):
        return relationship("TheoreticalGlycopeptide")


class TheoreticalPeptideProductIon(ProductIonBase, Base):
    __tablename__ = "TheoreticalPeptideProductIon"

    n_term = Column(Unicode(10))
    c_term = Column(Unicode(10))
    glycosylated = Column(Boolean, index=True)
    sequence_index = Column(Integer, index=True)

    @classmethod
    def from_fragment(cls, frag):
        inst = cls(
            name=frag.name,
            calculated_mass=frag.mass,
            n_term=frag.flanking_amino_acids[0].name,
            c_term=frag.flanking_amino_acids[1].name,
            glycosylated=('HexNAc' in frag.modification_dict),
            sequence_index=frag.position,
            ion_series=frag.name[0],
            neutral_loss=str(frag.neutral_loss if frag.neutral_loss is not None else "")
            )
        return inst

    def __repr__(self):
        return "<TheoreticalPeptideProductIon %s %s-%s %0.3f>" % (
            self.name, self.n_term, self.c_term, self.calculated_mass)


class TheoreticalGlycopeptideStubIon(ProductIonBase, Base):
    __tablename__ = "TheoreticalGlycopeptideStubIon"

    includes_peptide = Column(Boolean)

    def __repr__(self):
        return "<TheoreticalGlycopeptideStubIon %s %s %0.3f>" % (
            self.name, self.ion_series, self.calculated_mass)

    @classmethod
    def from_simple_fragment(cls, frag):
        inst = cls(
            name=frag.name,
            calculated_mass=frag.mass,
            ion_series=str(frag.kind),
            includes_peptide=fragment.IonSeries(frag.kind).includes_peptide,
            neutral_loss=str(frag.neutral_loss if frag.neutral_loss is not None else "")
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
