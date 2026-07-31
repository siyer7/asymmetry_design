import numpy as np
import pandas as pd


def get_pt_metadata(psychopy_df, patient, verbose=True):
    pt_data_dir = f'../../outputs/processed_data/2025{int(patient)}'
    pt_psychopy_df = psychopy_df.loc[psychopy_df['subj'] == patient].reset_index(drop=True)
    pt_neur_df = pd.read_parquet(f'{pt_data_dir}/df_neurs.parquet')
    pt_num_neurs = len(pt_neur_df)
    if verbose:
        print(f'patient={patient}', f'num_trials={len(pt_psychopy_df)}', f'num_neurons={pt_num_neurs}\n')
    return pt_psychopy_df, pt_neur_df, pt_num_neurs


def get_pt_epoch_spike_data(patient, epoch, verbose=True):
    pt_matrices_dir = f'../../outputs/processed_data/2025{int(patient)}/matrices'
    pt_epoch_spikes = np.load(f'{pt_matrices_dir}/{epoch}_spikes.npy', allow_pickle=True)
    pt_epoch_num_spikes = np.array([[len(pt_epoch_spikes[t, n]) for n in range(pt_epoch_spikes.shape[1])]
                                    for t in range(pt_epoch_spikes.shape[0])])
    pt_epoch_FRs = np.load(f'{pt_matrices_dir}/{epoch}_FRs.npy', allow_pickle=True)
    pt_epoch_bins = np.load(f'{pt_matrices_dir}/{epoch}_bin_centers.npy', allow_pickle=True)
    if verbose:
        print(f'patient={patient}, epoch={epoch}')
        print(f'spikes (trials, neurons): {pt_epoch_spikes.shape}')
        print(f'FRs (trials, neurons, bins): {pt_epoch_FRs.shape}')
    return pt_epoch_spikes, pt_epoch_num_spikes, pt_epoch_FRs, pt_epoch_bins


def get_contrast_in_epoch(psychopy_df, contrast, verbose=True):
    if contrast != 'shape_class':
        raise ValueError(f'Invalid contrast: {contrast}')
    cont_trials = [psychopy_df[psychopy_df['shape'] == 'curv'].index,
                   psychopy_df[psychopy_df['shape'] == 'flat'].index]
    cont_labels = ['curv', 'flat']
    if verbose:
        print(f'contrast: {contrast}')
    return cont_trials, cont_labels
