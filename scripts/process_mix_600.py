# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 09:34:04 2022

@author: DELL
"""


import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import inchi

from matchms.exporting import save_as_mgf
from matchms.importing import load_from_mgf
from core.identification import spectrum_processing
from core.msdial import load_MS_DIAL_Peaklist

'''
mix_msdial = load_MS_DIAL_Peaklist('example/600mix_pos.txt')
mix_msdial = [spectrum_processing(s) for s in mix_msdial]
mix_msdial = [s for s in mix_msdial if s is not None]
save_as_mgf(mix_msdial, 'example/600mix_pos.mgf')
spectrums = [s for s in load_from_mgf('example/600mix_pos.mgf')]
'''

get_key = lambda x: inchi.MolToInchiKey(Chem.MolFromSmiles(x))

mix_600 = pd.read_excel('example/SigmaMetList.xlsx')
mix_msdial_matching = pd.read_csv('example/600mix_pos.txt', sep = '\t')

true_smiles = mix_600['SMILES'].values
true_keys = [get_key(x)[:14] for x in true_smiles if Chem.MolFromSmiles(x) is not None]
ms_dial_matching_key = [s[:14] for s in mix_msdial_matching['InChIKey'].values if s is not np.nan]
ms_dial_matching = np.intersect1d(true_keys, ms_dial_matching_key)


import os
deepmass_path = 'example/DeepMASS_Result'
deepmass_top1_key, deepmass_top5_key = [], []
for f in os.listdir(deepmass_path):
    r = pd.read_csv(deepmass_path + '/' + f)
    for i in r.index:
        if i == 0:
            deepmass_top1_key.append(r.loc[i,'InChIKey'][:14])
        if i < 5:
            deepmass_top5_key.append(r.loc[i,'InChIKey'][:14])
deepmass_top1 = np.intersect1d(deepmass_top1_key, ms_dial_matching_key)
deepmass_top5 = np.intersect1d(deepmass_top5_key, ms_dial_matching_key)



