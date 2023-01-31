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

from utils.preprocess import spectrum_processing, get_true_precursor_mz_from_mass

files = os.listdir('data/CASMI2022_Priority')
challenge = pd.read_csv('solutions/casmi_2022_challenge_priority.csv')

challenge_ms = np.repeat(None, len(challenge))
for f in tqdm(files):
    f_ = 'data/CASMI2022_Priority/{}'.format(f)
    mode = f.split('_')[2][:3]
    challenge_ = challenge[challenge['File'] == f.split('.')[0]]
    spectrums = [s for s in load_from_mgf(f_)]
    for s in spectrums:
        if 'precursor_intensity' not in list(s.metadata.keys()):
            continue
        intensities = s.intensities / max(s.intensities)
        if np.min(np.abs(s.metadata['precursor_mz'] - challenge_['Precursor m/z (Da)'])) < 0.01:
            w = np.argmin(np.abs(s.metadata['precursor_mz'] - challenge_['Precursor m/z (Da)']))
            if abs(challenge_['RT [min]'].values[w] * 60 - s.metadata['retention_time']) < 5:
                i = challenge_.index[w]
                name = 'casmi_2022_challenge_priority_{}'.format(i)
                smi = challenge_.loc[i, 'SMILES']
                inchikey = challenge_.loc[i, 'InChIKey']
                precursor_type = challenge_.loc[i, 'Adduct'].replace(' ', '')
                mol = Chem.MolFromSmiles(smi)
                exactmass = AllChem.CalcExactMolWt(mol)
                formula = AllChem.CalcMolFormula(mol)
                precursor_mz = get_true_precursor_mz_from_mass(exactmass, precursor_type)
                
                s_new = Spectrum(mz = s.mz, intensities = s.intensities, metadata = {'precursor_mz': precursor_mz})
                s_new.set('name', name)
                s_new.set('formula', formula)
                s_new.set('smiles', smi)
                s_new.set('inchikey', inchikey)
                s_new.set('adduct', precursor_type)
                s_new.set('precursor_intensity', s.metadata['precursor_intensity'])
                if mode == 'pos':
                    s_new.set('ionmode', 'positive')
                else:
                    s_new.set('ionmode', 'negative')
                if challenge_ms[i] is None:
                    challenge_ms[i] = spectrum_processing(s_new)
                else:
                    if challenge_ms[i].metadata['precursor_intensity'] < s.metadata['precursor_intensity']:
                        challenge_ms[i] = spectrum_processing(s_new)

save_as_mgf(list(challenge_ms), 'save/casmi_2022_challenge_priority.mgf')

'''
for i, s in enumerate(challenge_ms):
    save_as_msp([challenge_ms[i]], 'save/casmi_2022_priority_msp/casmi_2022_challenge_priority_{}.msp'.format(i))
    with open('save/casmi_2022_priority_msp/casmi_2022_challenge_priority_{}.msp'.format(i)) as msp:
        lines = msp.readlines()
        lines = [l.replace('_', '') for l in lines]
        lines = [l.replace('ADDUCT', 'PRECURSORTYPE') for l in lines]
    with open('save/casmi_2022_priority_msp/casmi_2022_challenge_priority_{}.msp'.format(i), 'w') as msp:
        msp.writelines(lines)
'''