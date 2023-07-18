from typing import Union

# numeric state type
numeric = Union[int, float]

def in_to_feet(d: numeric) -> numeric:
    """
    Convert length from inch into feet
    """

    FEET_CONVERTION = 12
    
    return d / FEET_CONVERTION

def fahrenheit_to_rankine(T: numeric) -> numeric:
    RANKINE_CONSTANT = 460

    return T + RANKINE_CONSTANT

def gcc_to_lbm_per_ft3(density):
    DENSITY_CONSTANT = 62.428

    return density * DENSITY_CONSTANT

def sg_to_lbm_per_ft3(sg: numeric) -> numeric:
    """
    Convert spesific gravity into lbm/ft^3
    """

    return gcc_to_lbm_per_ft3(sg)