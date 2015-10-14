from glycresoft_sqlalchemy.data_model import Hierarchy as Registry
from glypy.structure.monosaccharide import ReducedEnd
from glypy.composition import glycan_composition, composition_transform

composition_source_type = Registry()
reduction_type = Registry()

reduction_type.register(True)(ReducedEnd)
reduction_type.register("reduced")(ReducedEnd)
reduction_type.register("deuteroreduced")(lambda:ReducedEnd("H[2]H"))
reduction_type.register(False)(lambda: False)
reduction_type.register(None)(lambda: False)


derivatization_type = Registry()
derivatization_type.register("permethylation")("methyl")


def set_reducing_end(glycan_obj, reduction_kind):
    glycan_obj.set_reducing_end(reduction_type[reduction_kind]())

def derivatize(glycan_obj, derivatization_kind):
    if derivatization_kind is not None:
        composition_transform.derivatize(glycan_obj, derivatization_type[derivatization_kind])
