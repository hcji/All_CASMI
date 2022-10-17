# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 09:13:11 2022

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

from utils.preprocess import spectrum_processing, get_true_precursor_mz_from_mass, get_adduct_from_mass_precursor

files = os.listdir('data/CASMI2016_Training')
challenge = pd.read_csv('solutions/casmi_2016_challenge_training.csv')

challenge_ms = np.repeat(None, len(challenge))
for i, f in enumerate(tqdm(files)):
    f_ = 'data/CASMI2016_Training/{}'.format(f)
    spectrum = [s for s in load_from_mgf(f_)][0]
    
    name = challenge.loc[i, 'challengename']
    smi = challenge.loc[i, 'SMILES']
    inchikey = challenge.loc[i, 'INCHIKEY']
    mol = Chem.MolFromSmiles(smi)
    exactmass = AllChem.CalcExactMolWt(mol)
    formula = AllChem.CalcMolFormula(mol)
    precursor_mz = challenge.loc[i, 'PRECURSOR_MZ']
    ion_mode = challenge.loc[i, 'ION_MODE'].lower().replace(' ', '')
    
    adduct = get_adduct_from_mass_precursor(exactmass, precursor_mz, ion_mode, tolerence = 0.1)
    precursor_mz = get_true_precursor_mz_from_mass(exactmass, adduct)
    
    s_new = Spectrum(mz = spectrum.mz, intensities = spectrum.intensities, metadata = {'precursor_mz': precursor_mz})
    s_new.set('name', name)
    s_new.set('formula', formula)
    s_new.set('smiles', smi)
    s_new.set('inchikey', inchikey)
    s_new.set('adduct', adduct)
    s_new.set('ionmode', ion_mode)
    
    challenge_ms[i] = spectrum_processing(s_new)
    
save_as_mgf(list(challenge_ms), 'save/casmi_2016_challenge_training.mgf')

'''
for i, s in enumerate(challenge_ms):
    save_as_msp([challenge_ms[i]], 'save/casmi_2016_training_msp/casmi_2016_challenge_{}.msp'.format(i))
    with open('save/casmi_2016_training_msp/casmi_2016_challenge_{}.msp'.format(i)) as msp:
        lines = msp.readlines()
        lines = [l.replace('_', '') for l in lines]
        lines = [l.replace('ADDUCT', 'PRECURSORTYPE') for l in lines]
    with open('save/casmi_2016_training_msp/casmi_2016_challenge_{}.msp'.format(i), 'w') as msp:
        msp.writelines(lines)
'''    