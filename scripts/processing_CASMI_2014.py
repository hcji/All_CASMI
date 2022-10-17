# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 09:13:11 2022

@author: DELL
"""

import requests
import numpy as np
import pandas as pd

from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from matchms import Spectrum
from matchms.importing import load_from_mgf
from matchms.exporting import save_as_mgf, save_as_msp

from utils.preprocess import spectrum_processing, get_true_precursor_mz_from_mass, get_adduct_from_mass_precursor

challenge = pd.read_csv('solutions/casmi_2014_challenge.csv')

challenge_ms = np.repeat(None, len(challenge))
for i in tqdm(challenge.index):
    link = 'http://casmi-contest.org/2014/Challenge2014/Challenge{}/{}_MSMS.txt'.format(i+1, i+1)
    assay = 'http://casmi-contest.org/2014/Challenge2014/Challenge{}/{}_Experimental_details.txt'.format(i+1, i+1)
    precursor = 'http://casmi-contest.org/2014/Challenge2014/Challenge{}/{}_MS.txt'.format(i+1, i+1)
    
    try:
        ms = requests.get(link).text.split('\n')
        ms = [s.replace('\r', '') for s in ms if s != '']
        ms = [s.replace(' ', '') for s in ms if s != '']
        mz = np.array([float(s.split('\t')[0]) for s in ms if s != ''])
    except:
        continue
    intensities = np.array([float(s.split('\t')[1]) for s in ms if s != ''])
    
    details = requests.get(assay).text
    if 'positive' in details.lower().split(' '):
        ion_mode = 'positive'
    else:
        ion_mode = 'negative'

    for r in requests.get(precursor).text.split('\n'):
        try:
            precursor_mz = float(r.split('\t')[0])
            continue
        except:
            pass
        
    name = challenge.loc[i, 'Challenge']
    inchi = challenge.loc[i, 'InChI']
    mol = Chem.MolFromInchi(inchi)
    smi = Chem.MolToSmiles(mol)
    inchikey = challenge.loc[i, 'InChIkey']
    exactmass = AllChem.CalcExactMolWt(mol)
    formula = AllChem.CalcMolFormula(mol)
    
    adduct = get_adduct_from_mass_precursor(exactmass, precursor_mz, ion_mode, tolerence = 0.1)
    if adduct == '':
        if ion_mode == 'positive':
            adduct = '[M+H]+'
        else:
            adduct = '[M-H]-'
    precursor_mz = get_true_precursor_mz_from_mass(exactmass, adduct)
    
    s_new = Spectrum(mz = mz, intensities = intensities, metadata = {'precursor_mz': precursor_mz})
    s_new.set('name', name)
    s_new.set('formula', formula)
    s_new.set('smiles', smi)
    s_new.set('inchikey', inchikey)
    s_new.set('adduct', adduct)
    s_new.set('ionmode', ion_mode)
    
    challenge_ms[i] = spectrum_processing(s_new)
challenge_ms = [s for s in challenge_ms if s is not None]

save_as_mgf(list(challenge_ms), 'save/casmi_2014_challenge_test.mgf')

'''
for i, s in enumerate(challenge_ms):
    save_as_msp([challenge_ms[i]], 'save/casmi_2014_msp/casmi_2016_challenge_{}.msp'.format(i))
    with open('save/casmi_2014_msp/casmi_2014_challenge_{}.msp'.format(i)) as msp:
        lines = msp.readlines()
        lines = [l.replace('_', '') for l in lines]
        lines = [l.replace('ADDUCT', 'PRECURSORTYPE') for l in lines]
    with open('save/casmi_2014_msp/casmi_2014_challenge_{}.msp'.format(i), 'w') as msp:
        msp.writelines(lines)
'''    