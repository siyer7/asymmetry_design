# utils_data.py

import glob
import ast
import numpy as np
import pandas as pd

ROWS_TO_KEEP = {'rhia': 40}
DEFAULT_ROWS = 240

CONDITION_MAP = {
    'A': {'baseline': 'baseline', 'curv_comp': 'pen_comp', 'flat_comp': 'rew_comp'},
    'C': {'baseline': 'baseline', 'curv_comp': 'pen_comp', 'flat_comp': 'rew_comp'},
    'B': {'baseline': 'baseline', 'curv_comp': 'rew_comp', 'flat_comp': 'pen_comp'},
    'D': {'baseline': 'baseline', 'curv_comp': 'rew_comp', 'flat_comp': 'pen_comp'},
}

def load_and_process_subject(subj_id):
    nrows = ROWS_TO_KEEP.get(subj_id, DEFAULT_ROWS)
    file_path = glob.glob(f'../results/*{subj_id}*')[0]
    df = pd.read_csv(
        file_path,
        nrows=nrows,
        converters={'positions': ast.literal_eval}
    )
    df = df.sort_values(by='trial_key').reset_index(drop=True)
    df['flipped'] = np.where(df['shape_order'] == 'curv_flat', 0, 1)
    df['chosen_pos'] = df['positions'].str[-1]
    sess_letter = df.loc[0, 'sess_type'][0]
    sign = +1 if sess_letter in ('A', 'C') else -1
    df['condition'] = df['condition'].map(CONDITION_MAP[sess_letter])
    df['chosen_pos_flipped'] = sign * df['chosen_pos']
    df['div_true'] = sign * df['div_pos']
    df['class_true'] = np.where(df['valence'] == 'rew', 1, 0)
    df['class_pred'] = np.where(df['chosen_pos_flipped'] > df['div_true'], 1, 0)
    df['chosen_pos_aligned'] = df['chosen_pos_flipped']
    df['stim_pos_aligned'] = sign * df['stim_pos']
    df['div_true_aligned'] = df['div_true']
    return df
