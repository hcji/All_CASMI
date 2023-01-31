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
from rdkit.Chem import AllChem, inchi
from matchms.importing import load_from_mgf

spectrums = [s for s in load_from_mgf('save/casmi_2016_challenge_test.mgf')]

sirius_path = "comparison/casmi2016_test/SIRIUS"
sirius_files = [name for name in os.listdir(sirius_path) if os.path.isdir(os.path.join(sirius_path, name))]
sirius_index = [1 + int(i.split('_')[-2]) for i in sirius_files]

deepmass_path = "D:/DeepMASS2_GUI/result/casmi2016_test"
deepmass_files = [name for name in os.listdir(deepmass_path)]
deepmass_index = [int(i.split('.')[-2][-3:]) for i in deepmass_files]

sirius_candidate = []
ranking_result = []
for s in tqdm(spectrums):
    name = s.metadata['compound_name']
    index = int(name.split('-')[-1])
    true_key = s.metadata['inchikey'][:14]
    
    sirius_file = "/{}/structure_candidates.tsv".format(sirius_files[sirius_index.index(index)])
    sirius_file = sirius_path + sirius_file
    try:
        sirius_result = pd.read_csv(sirius_file, sep='\t')
    except:
        sirius_result = None
    
    # collect candidates
    if sirius_result is not None:
        for j in sirius_result.index:
            nam = sirius_result.loc[j, 'name']
            smi = sirius_result.loc[j, 'smiles']
            mol = Chem.MolFromSmiles(smi)
            formula = AllChem.CalcMolFormula(mol)
            exactmass = AllChem.CalcExactMolWt(mol)
            inchikey = inchi.MolToInchiKey(mol)
            sirius_candidate.append([nam, smi, formula, inchikey, exactmass])
    
    # rank of sirius
    if sirius_result is not None:
        sirius_n = len(sirius_result)
        sirius_key = np.array([k for k in sirius_result['InChIkey2D']])
        sirius_rank = np.where(sirius_key == true_key)[0]
        if len(sirius_rank) == 0:
            sirius_rank = float('inf')
        else:
            sirius_rank = sirius_rank[0] + 1
    else:
        sirius_n = 0
        sirius_rank = float('inf')
    
    # rank of deepmass
    deepmass_file = "/{}".format(deepmass_files[deepmass_index.index(index)])
    deepmass_file = deepmass_path + deepmass_file
    deepmass_result = pd.read_csv(deepmass_file)
    deepmass_key = np.array([k[:14] for k in deepmass_result['InChIKey']])
    deepmass_n = len(deepmass_key)
    deepmass_rank = np.where(deepmass_key == true_key)[0]
    if len(deepmass_rank) == 0:
        deepmass_rank = float('inf')
    else:
        deepmass_rank = deepmass_rank[0] + 1

    ranking_result.append([name, true_key, sirius_rank, sirius_n, deepmass_rank, deepmass_n])

ranking_result = pd.DataFrame(ranking_result, columns = ['Challenge', 'True Inchikey2D', 'SIRIUS Ranking', 'SIRIUS Candidates',
                                                         'DeepMASS Ranking', 'DeepMASS Candidates'])



'''
biodatabase = pd.read_csv('D:/DeepMASS2_GUI/data/MsfinderStructureDB-VS15-plus-GNPS.csv')
sirius_candidate = pd.DataFrame(sirius_candidate, columns=['Title', 'SMILES', 'Formula', 'InChIkey', 'Exact mass'])
sirius_candidate['Short InChIKey'] = [k[:14] for k in sirius_candidate['InChIkey']]
kk = []
for i in tqdm(sirius_candidate.index):
    inchikey = sirius_candidate.loc[i, 'Short InChIKey']
    if inchikey not in biodatabase['Short InChIKey'].values:
        kk.append(i)
        
biodatabase = biodatabase.append(sirius_candidate.loc[kk,:])
biodatabase = biodatabase.sort_values(by = 'Exact mass')

biodatabase.to_csv('D:/DeepMASS2_GUI/data/MsfinderStructureDB-VS15-plus-GNPS.csv', index=False)
'''