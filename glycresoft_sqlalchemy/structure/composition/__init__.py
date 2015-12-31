# from .composition import Composition, calculate_mass, std_ion_comp

from ...utils.memoize import memoize

from glypy.composition.composition import Composition, calculate_mass, ChemicalCompositionError


@memoize()
def composition_to_mass(formula):
    '''Fast, low allocation mass computation'''
    return Composition(formula).mass


class mass(float):
    def __call__(self):
        return self

    @property
    def fmass(self):
        return self
