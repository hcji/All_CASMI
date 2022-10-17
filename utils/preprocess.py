# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 10:22:09 2022

@author: DELL
"""

import os
import numpy as np
import pandas as pd

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from matchms import Spectrum
from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf, save_as_msp
import matchms.filtering as msfilters


ADDUCT_MASSES = {
    # POSITIVE
    '[M+H]+': 1.007276,
    '[M-H2O+H]+': 1.007276 - 18.010565, '[M+H-H2O]+': 1.007276 - 18.010565,
    '[M-2H2O+H]+': 1.007276 - 2 * 18.010565, '[M+H-2H2O]+': 1.007276 - 2 * 18.010565,
    '[M+Na]+': 22.98922,
    '[M]+': 0,
    '[M+NH4]+': 18.03383,
    '[M+H-NH3]+': 1.007276 - 17.026549,
    # NEGATIVE
    '[M-H]-': -1.007276,
    '[M+Cl]-': 34.969402,
    '[M+FA-H]-': 44.9982, '[M-H+FA]-': 44.9982,
    '[M]-': 0,
    '[M-2H]-': -1.007276
}


def spectrum_processing(s):
    """This is how one would typically design a desired pre- and post-
    processing pipeline."""
    s = msfilters.default_filters(s)
    s = msfilters.add_parent_mass(s)
    s = msfilters.normalize_intensities(s)
    s = msfilters.select_by_mz(s, mz_from=0, mz_to=2000)
    s = msfilters.add_losses(s, loss_mz_from=10.0, loss_mz_to=200.0)
    return s


def get_true_precursor_mz_from_mass(mass, precursor_type):
    """
    Calculate precursor mz based on (exact or monoisotopic) mass and precursor type
    :param mass: scalar, mass of a compound, e.g. monoisotopic or exact mass.
    :param precursor_type: string, precursor type, e.g. '[M+H]+'
    :return: scalar, ion mass / precursor mz
    """
    try:
        return mass + ADDUCT_MASSES[precursor_type]
    except KeyError:
        raise KeyError("Unsupported precursor-type '%s'." % precursor_type)


def get_adduct_from_mass_precursor(exactmass, precursor_mz, ion_mode, tolerence):
    if ion_mode == 'positive':
        k = [k for k in ADDUCT_MASSES.keys() if k[-1] == '+']
    else:
        k = [k for k in ADDUCT_MASSES.keys() if k[-1] == '-']
    adduct = ''
    diff = precursor_mz - exactmass
    for kk in k:
        if abs(diff - ADDUCT_MASSES[kk]) <= tolerence:
            adduct = kk
            break
    return adduct
        
    

