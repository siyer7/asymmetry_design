import numpy as np
import argparse
import pdb
import os
import sys
import pickle
import pandas
import json
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from sklearn.manifold import MDS
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from matplotlib.patches import Ellipse
import pandas as pd
import statsmodels.api as sm


rootpath = os.path.join(os.getcwd(), '/home/nuttidalab/Documents/verbalInstructionTask/dataAnalysis_shared/')
sys.path.append(rootpath)
from helpers.custom_loadmat import custom_loadmat
from helpers.get_rates_from_spikes import get_rates_from_spikes


AREA_LABELS = {
    1: 'MTL', 2: 'MTL', 3: 'MTL', 4: 'MTL',
    5: 'MFC', 6: 'MFC', 7: 'MFC', 8: 'MFC',
    9: 'vmPFC', 10: 'vmPFC',
    11: 'VTC', 12: 'VTC'
}


def flatten_units(allUnits, sessions):
    """
    allUnits: list over sessions, each entry is sessionData.neuralData.unitCell
             which is typically an array-like of unit structs.
    Returns:
        units_flat: list of unit structs in concatenation order
        unit_session_idx: array mapping global unit -> session index
        unit_within_session_idx: array mapping global unit -> within-session unit index
    """
    units_flat = []
    unit_session_idx = []
    unit_within_session_idx = []

    for sI, unitCell in enumerate(allUnits):
        # unitCell may be list-like or np array; iterate robustly
        nU = len(unitCell)
        for uI in range(nU):
            units_flat.append(unitCell[uI])
            unit_session_idx.append(sI)
            unit_within_session_idx.append(uI)

    return units_flat, np.array(unit_session_idx), np.array(unit_within_session_idx)

def parse_unit_info(unitStruct):
    ui = np.asarray(unitStruct.get('unitInfo', np.array([]))).ravel()
    if ui.size < 4:
        return None
    # unitInfo is typically [unit#, wire#, cluster#, areaCode]
    unit_num, wire_num, cluster_num, area_code = ui[:4].astype(int)
    area_label = AREA_LABELS.get(area_code, f"area{area_code}")
    return unit_num, wire_num, cluster_num, area_code, area_label

def get_hist_sliding(timestamps, period, bin_size, step_size):
    """
    Python analogue of getHist(timestamps, Period, binSize, stepSize).

    timestamps: 1D array-like of spike times (seconds)
    period: (t_start, t_end) in seconds
    bin_size: bin size in seconds
    step_size: fraction of bin_size for sliding window shift (e.g., 1/16)
              step = bin_size * step_size

    Returns:
        density: 1D array length = number of kept bins within period
        bin_centers: 1D array of bin centers (seconds), aligned with density
    """
    t0, t1 = float(period[0]), float(period[1])
    step = bin_size * step_size

    # MATLAB: bins = Period(1):binSize*stepSize:Period(2)+binSize;
    edges = np.arange(t0, t1 + bin_size + 1e-12, step, dtype=float)

    # MATLAB: keep = bins>=Period(1) & bins<=Period(2);
    keep = (edges >= t0 - 1e-12) & (edges <= t1 + 1e-12)

    ts = np.asarray(timestamps, dtype=float).ravel()
    ts = ts[np.isfinite(ts)]

    # MATLAB histc(timestamps, bins) returns counts for each edge bin.
    # np.histogram counts into intervals [edges[i], edges[i+1]); so to mimic histc-ish
    # we can histogram on edges with length-1 output, which is what we want anyway.
    # Use histogram over step edges:
    counts, _ = np.histogram(ts, bins=edges)

    # MATLAB does sliding window sum:
    # for i=1:round(1/stepSize)-1
    #   density = density + [temp(i+1:end), NaN(1,i)];
    # end
    # Here, round(1/stepSize) should be int for typical step_size like 1/16.
    temp = counts.astype(float)
    density = temp.copy()

    n_shifts = int(round(1.0 / step_size)) - 1
    for i in range(n_shifts):
        shifted = np.concatenate([temp[(i+1):], np.full((i+1,), np.nan)])
        density = density + shifted

    # Crop to keep window and remove NaNs
    # We need density aligned to edges[:-1] bins; keep mask corresponds to edges,
    # so use keep[:-1] for bins.
    keep_bins = keep[:-1]
    density = density[keep_bins]

    # Bin centers for plotting: edges bin centers for step histogram
    centers = edges[:-1] + step / 2.0
    centers = centers[keep_bins]
    return density, centers

def gaussian_kernel(sig_points: float, n: int) -> np.ndarray:
    """
    Match MATLAB getGaussianKernel(sig, n)
    sig_points: sigma in units of "number of datapoints"
    n: half-width in datapoints; total kernel length = 2n+1
    """
    k = np.arange(-n, n + 1, dtype=float)
    g = (1.0 / np.sqrt(2 * np.pi * sig_points**2)) * np.exp(-(k**2) / (sig_points**2))
    g = g / np.sum(g)
    return g

def get_psth(timestamps_list, period, bin_size=0.2, step_size=1/16, smooth=False,
             smooth_sigma_s=0.4, smooth_halfwidth_bins=2):
    """
    Python analogue of MATLAB getPSTH(timestamps, period, binsize, stepSize, smooth)

    timestamps_list: list of trials, each is 1D array of spike times (seconds)
    period: (t_start, t_end) in seconds
    bin_size: seconds (MATLAB default 200 ms)
    step_size: fraction of bin_size (MATLAB default 1/16)
    smooth: whether to convolve PSTH with Gaussian kernel
    smooth_sigma_s: sigma in seconds for smoothing (MATLAB uses 0.4 sec)
    smooth_halfwidth_bins: n in getGaussianKernel(sig, n) (MATLAB uses 2)

    Returns:
        psth_rate_hz: 1D array (Hz), same length as centers
        density_rate_hz: 2D array (nTrials x nBins) (Hz)
        centers: 1D array of time centers (seconds)
    """
    # Build per-trial densities
    densities = []
    centers_ref = None

    for ts in timestamps_list:
        d, centers = get_hist_sliding(ts, period, bin_size, step_size)
        densities.append(d)
        if centers_ref is None:
            centers_ref = centers
        else:
            # Safety: require same centers
            if len(centers) != len(centers_ref) or np.max(np.abs(centers - centers_ref)) > 1e-9:
                raise RuntimeError("Histogram centers mismatch across trials; check inputs/period.")

    if len(densities) == 0:
        # no trials
        centers_ref = np.arange(period[0], period[1], bin_size * step_size) + (bin_size * step_size)/2
        density = np.zeros((0, len(centers_ref)), dtype=float)
        psth = np.zeros(len(centers_ref), dtype=float)
    else:
        density = np.vstack(densities)  # (nTrials, nBins)
        psth = np.nanmean(density, axis=0)

    # Convert to rate in Hz: MATLAB later divides by binSize (seconds)
    density_rate = density / bin_size
    psth_rate = psth / bin_size

    # Optional smoothing
    if smooth:
        # MATLAB: gaussKernel = getGaussianKernel(0.4/binSize, 2)
        # Here binSize is in seconds already.
        sig_points = smooth_sigma_s / bin_size
        gk = gaussian_kernel(sig_points=sig_points, n=smooth_halfwidth_bins)
        psth_rate = np.convolve(psth_rate, gk, mode='same')

    return psth_rate, density_rate, centers_ref

def _get_alignment_spikes(unitStruct, alignment: str):
    """
    alignment: one of {'trial','instruction','stim1','stim2','decision'}
    Returns: list of 1D arrays (spike times) with length = nTrials
    """
    key_map = {
        'trial': 'trialReferencedSpikes',
        'instruction': 'instructionReferencedSpikes',
        'stim1': 'stim1ReferencedSpikes',
        'stim2': 'stim2ReferencedSpikes',
        'decision': 'decisionReferencedSpikes',
        'stim1_stim2': 'stim1_stim2ReferencedSpikes',
    }
    if alignment not in key_map:
        raise ValueError(f"Unknown alignment '{alignment}'. Valid: {list(key_map.keys())}")

    k = key_map[alignment]
    if k not in unitStruct:
        raise KeyError(f"Unit struct missing field '{k}'. Available keys: {list(unitStruct.keys())[:20]} ...")

    spikes = unitStruct[k]

    # Ensure python list of arrays
    # Depending on custom_loadmat, spikes may be list-of-lists, object array, etc.
    spikes_list = []
    for tr in spikes:
        if tr is None:
            spikes_list.append(np.array([], dtype=float))
        else:
            a = np.asarray(tr).astype(float).ravel()
            # Some loaders represent empty as [nan] or [1]; treat non-finite as empty
            a = a[np.isfinite(a)]
            spikes_list.append(a)
    return spikes_list

def plot_unit_raster_psth(
    unitStruct,
    behaviorData,
    alignment='stim2',
    trialGroups=None,
    group_labels=None,
    timeLimits=None,
    gray_boxes=None,
    bin_size=0.2,
    title=None,
    smooth_pad_s=0.5,   # <-- NEW: padding (seconds) used only for smoothing
):
    """
    unitStruct: a single unit struct (dict-like)
    behaviorData: taskStruct (dict-like) for the session, containing antiTask, respKey, correctResponses, etc.
    trialGroups: list of arrays of trial indices (0-based)
    group_labels: list of strings (same length as trialGroups)
    timeLimits: (tmin, tmax) or None
    gray_boxes: list of (t0, t1) or None
    """
    spikes_list = _get_alignment_spikes(unitStruct, alignment)

    nTrials = len(spikes_list)

    # ---- Default trial groups:
    if trialGroups is None:
        antiTask = np.asarray(behaviorData['antiTask']).ravel()
        respKey = np.asarray(behaviorData['respKey']).ravel()
        trialConditions = np.asarray(behaviorData['trialConditions']).ravel()
        trialConditions_stim1_stim2 = np.asarray(behaviorData['trialConditions_stim1stim2']).ravel()
        correctResponses = np.asarray(behaviorData['correctResponses']).ravel()
        correctTrials = (respKey == correctResponses)        
        targetFeatures_stim1stim2 = behaviorData['targetFeatures_stim1stim2']
        arr = np.array(behaviorData['targetFeatures_stim1stim2'])
        mask_OldNew = np.isin(arr, ['Old', 'New'])
        wasAgeRelevant_stim1stim2 = np.array(behaviorData['wasAgeRelevant']+behaviorData['wasAgeRelevant'])
        wasCountRelevant_stim1stim2 = np.array(behaviorData['wasCountRelevant']+behaviorData['wasCountRelevant'])

        # Early / late instruction
        #g1 = np.where((trialConditions_stim1_stim2 == 1))[0]
        #g2 = np.where((trialConditions_stim1_stim2 == 2))[0]

        #g1 = np.where((antiTask == 0) & correctTrials)[0]
        #g2 = np.where((antiTask == 1) & correctTrials)[0]
        
        # First half / second half
        #g1 = np.arange(192)
        #g2 = np.arange(192, 384)

        # Target features (Old/New)
        g1 = np.where((wasAgeRelevant_stim1stim2 == 1))[0]
        g2 = np.where((wasAgeRelevant_stim1stim2 == 0))[0]

        # Target features (low count/high count)
        #g1 = np.where((wasCountRelevant_stim1stim2 == 1))[0]
        #g2 = np.where((wasCountRelevant_stim1stim2 == 0))[0]

        trialGroups = [g1, g2]
        #group_labels = ['Task (correct)', 'Anti-task (correct)']
        #group_labels = ['Task (correct)', 'Anti-task (correct)']

    if group_labels is None:
        group_labels = [f'group{i+1}' for i in range(len(trialGroups))]

    # Default time limits / gray boxes (match your MATLAB defaults)
    if timeLimits is None:
        if alignment in ['instruction', 'stim1', 'stim2', 'stim1_stim2']:
            timeLimits = (-1.0, 2.0)
        elif alignment == 'decision':
            timeLimits = (-2.0, 1.75)
        else:
            timeLimits = (-0.25, 10.0)

    if gray_boxes is None:
        # Your MATLAB uses [0 0.05] for most alignments; keep it configurable
        gray_boxes = [(0.0, 0.05)]

    # Figure
    fig, (ax_r, ax_p) = plt.subplots(2, 1, figsize=(3.0, 6.0), sharex=True)
    if title is not None:
        fig.suptitle(title, y=0.98)

    # Colors: match your MATLAB spirit (blue/black by default for 2 groups)
    default_colors = ['tab:blue', 'k', 'tab:purple', 'tab:cyan']
    colors = [default_colors[i % len(default_colors)] for i in range(len(trialGroups))]

    # ---- Raster
    cumulativeTrial = 0
    for gI, trials in enumerate(trialGroups):
        for tID in trials:
            if tID < 0 or tID >= nTrials:
                continue
            cumulativeTrial += 1
            ts = spikes_list[tID]
            if ts.size == 0:
                continue
            ys = np.full(ts.shape, cumulativeTrial, dtype=float)
            ax_r.plot(ts, ys, '.', ms=3, color=colors[gI])

    ax_r.set_ylabel('Trial nr. (reordered)')
    ax_r.set_xlabel(f'Time from {alignment} (s)')
    ax_r.set_xlim(timeLimits)

    # gray boxes behind raster
    y0, y1 = ax_r.get_ylim()
    for (t0, t1) in gray_boxes:
        ax_r.axvspan(t0, t1, color='0.87', zorder=-10)
    ax_r.set_ylim(y0, y1)

    # ---- PSTH
    tmin, tmax = timeLimits
    period_display = (tmin, tmax)
    period_padded  = (tmin - smooth_pad_s, tmax + smooth_pad_s)  # <-- padded window

    for gI, trials in enumerate(trialGroups):
        group_spikes = [spikes_list[tID] for tID in trials if 0 <= tID < len(spikes_list)]

        psth_pad, density_pad, centers_pad = get_psth(
            group_spikes,
            period=period_padded,   # <-- compute on padded
            bin_size=bin_size,      # <-- use function arg (was hardcoded 0.2)
            step_size=1/16,
            smooth=True,
            smooth_sigma_s=0.4,
            smooth_halfwidth_bins=2
        )

        # ---- CROP BACK to display limits (so x-axis shows only timeLimits)
        keep = (centers_pad >= tmin) & (centers_pad <= tmax)
        centers = centers_pad[keep]
        psth    = psth_pad[keep]
        density = density_pad[:, keep] if density_pad.ndim == 2 else density_pad[keep]

        # SEM across trials (density is Hz)
        if density.shape[0] > 1:
            sem = np.nanstd(density, axis=0, ddof=1) / np.sqrt(density.shape[0])
        else:
            sem = np.zeros_like(psth)

        ax_p.plot(centers, psth, lw=2, color=colors[gI], label=group_labels[gI])
        ax_p.fill_between(centers, psth - sem, psth + sem, color=colors[gI], alpha=0.2, linewidth=0)

    ax_p.set_ylabel('Firing rate (Hz)')
    ax_p.set_xlabel(f'Time from {alignment} (s)')
    ax_p.set_xlim(timeLimits)          # <-- still only show requested limits
    ax_p.legend(frameon=False)

    # gray boxes behind PSTH
    y0, y1 = ax_p.get_ylim()
    for (t0, t1) in gray_boxes:
        ax_p.axvspan(t0, t1, color='0.87', zorder=-10)
    ax_p.set_ylim(y0, y1)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    return fig



# Visualizing rasters
basefolder = '/home/nuttidalab/Documents/verbalInstructionTask/patientData/'
sessions = ['P81CS_102922','P82CS_011723','P86CS_072723','P90CS_091623','P92CS_021424',
            'P97CS_051324','P98CS_080524','P99CS_081724','P100CS_091524','P102CS_112524','P103CS_111724']

allUnits = []
taskStructs = []

for sI, session in enumerate(sessions):
    data = custom_loadmat(os.path.join(basefolder, session, 'sessionData.mat'))
    taskStructs.append(data['sessionData']['taskStruct'])
    allUnits.append(data['sessionData']['neuralData']['unitCell'])

units_flat, unit_session_idx, unit_within_session_idx = flatten_units(allUnits, sessions)


# Example: pick 5 random "other" units to plot
# Candidate indices for raster plotting

#to_plot = [551,199,352,737,538] # manually selected for nice examples
to_plot = [214,601,45,593,512] # max dec_hid MFC unit


out_dir = os.path.join(basefolder, "rasters")
os.makedirs(out_dir, exist_ok=True)

alignment = 'stim1_stim2'

for gi in to_plot:
    u = units_flat[gi]
    sI = unit_session_idx[gi]
    beh = taskStructs[sI]

    info = parse_unit_info(u)
    if info is None:
        title = f"{sessions[sI]} | globalUnit={gi} | withinSession={unit_within_session_idx[gi]}"
    else:
        unit_num, wire_num, cluster_num, area_code, area_label = info
        title = (f"{sessions[sI]} | globalUnit={gi} | \n unit {unit_num}-{wire_num}-cl{cluster_num} "
                 f"| area={area_label} ({area_code}) \n"
                 f"| withinSession= {unit_within_session_idx[gi]}" )
    
    fig = plot_unit_raster_psth(
        unitStruct=u,
        behaviorData=beh,
        alignment=alignment,
        timeLimits=(-1, 2),
        gray_boxes=[(0.0, 0.05)],
        bin_size=0.2,
        title=title
    )

    out_png = os.path.join(out_dir, f"raster_{sessions[sI]}_g{gi}_{alignment}.png")
    out_pdf = os.path.join(out_dir, f"raster_{sessions[sI]}_g{gi}_{alignment}.pdf")
    fig.savefig(out_png, dpi=150)
    fig.savefig(out_pdf, dpi=150)
    plt.close(fig)
    print(out_png)
