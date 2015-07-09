from .composition import Composition, calculate_mass, std_ion_comp

from ...utils.memoize import memoize

@memoize()
def composition_to_mass(formula):
    '''Fast, low allocation mass computation'''
    return Composition(formula).mass
