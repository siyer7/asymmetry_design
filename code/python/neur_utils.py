import numpy as np
import pandas as pd


def get_pt_metadata(psychopy_df, patient, verbose=True):
    pt_data_dir = f'../../results/2025{int(patient)}/records/processed_data'
    pt_psychopy_df = psychopy_df.loc[psychopy_df['subj'] == patient].reset_index(drop=True)
    pt_neur_df = pd.read_parquet(f'{pt_data_dir}/df_neurs.parquet')
    pt_num_neurs = len(pt_neur_df)
    if verbose:
        print(f'patient={patient}', f'num_trials={len(pt_psychopy_df)}', f'num_neurons={pt_num_neurs}\n')
    return pt_psychopy_df, pt_neur_df, pt_num_neurs


def get_pt_epoch_spike_data(patient, epoch, verbose=True):
    pt_data_dir = f'../../results/2025{int(patient)}/records/processed_data'
    pt_epoch_spikes = np.load(f'{pt_data_dir}/{epoch}_spikes.npy', allow_pickle=True)
    pt_epoch_num_spikes = np.array([[len(pt_epoch_spikes[t, n]) for n in range(pt_epoch_spikes.shape[1])]
                                    for t in range(pt_epoch_spikes.shape[0])])
    pt_epoch_FRs = np.load(f'{pt_data_dir}/{epoch}_FRs.npy', allow_pickle=True)
    pt_epoch_bins = np.load(f'{pt_data_dir}/{epoch}_bin_centers.npy', allow_pickle=True)
    if verbose:
        print(f'patient={patient}, epoch={epoch}')
        print(f'spikes (trials, neurons): {pt_epoch_spikes.shape}')
        print(f'FRs (trials, neurons, bins): {pt_epoch_FRs.shape}')
    return pt_epoch_spikes, pt_epoch_num_spikes, pt_epoch_FRs, pt_epoch_bins


def get_contrast_in_epoch(psychopy_df, contrast, verbose=True):
    cont_trials, cont_labels = [], []
    if contrast == 'boundary_context':
        for cond in ['curv_comp', 'baseline', 'flat_comp']:
            cont_trials.append(psychopy_df[psychopy_df['condition'] == cond].index)
            cont_labels.append(cond)
    elif contrast == 'shape_class':
        cont_trials.append(psychopy_df[psychopy_df['shape'] == 'curv'].index)
        cont_trials.append(psychopy_df[psychopy_df['shape'] == 'flat'].index)
        cont_labels += ['curv', 'flat']
    elif contrast == 'valence':
        cont_trials.append(psychopy_df[psychopy_df['stim_pos_aligned'] > psychopy_df['div_pos_aligned']].index)
        cont_trials.append(psychopy_df[psychopy_df['stim_pos_aligned'] < psychopy_df['div_pos_aligned']].index)
        cont_labels += ['gain', 'loss']
    elif contrast == 'ambiguity':
        cont_trials.append(psychopy_df[~psychopy_df['uncertainty']].index)
        cont_trials.append(psychopy_df[psychopy_df['uncertainty']].index)
        cont_labels += ['certain', 'uncertain']
    elif contrast == 'resp_dir':
        cont_trials.append(psychopy_df[psychopy_df['chosen_pos'] > psychopy_df['div_pos']].index)
        cont_trials.append(psychopy_df[psychopy_df['chosen_pos'] < psychopy_df['div_pos']].index)
        cont_labels += ['right', 'left']
    elif contrast == 'normed_RT':
        median_rt = psychopy_df['normed_RT'].median()
        cont_trials.append(psychopy_df[psychopy_df['normed_RT'] <= median_rt].index)
        cont_trials.append(psychopy_df[psychopy_df['normed_RT'] > median_rt].index)
        cont_labels += ['fast', 'slow']
    elif contrast == 'outcome':
        for val, label in zip([3, 1, -1, -3], ['3 coins', '1 coin', '-1 coin', '-3 coins']):
            cont_trials.append(psychopy_df[psychopy_df['outcome'] == val].index)
            cont_labels.append(label)
    elif contrast == 'binarized_context':
        cont_trials.append(psychopy_df[psychopy_df['condition'] == 'baseline'].index)
        cont_trials.append(psychopy_df[psychopy_df['condition'] != 'baseline'].index)
        cont_labels += ['base', 'comp']
    else:
        raise ValueError(f'Invalid contrast: {contrast}')
    if verbose:
        print(f'contrast: {contrast}')
    return cont_trials, cont_labels
