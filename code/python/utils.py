def plot_style(context):
    import seaborn as sns, matplotlib.pyplot as plt, matplotlib as mpl

    sns.set(context=context, style='ticks', palette='deep')
    plt.rcParams['svg.fonttype'] = 'none'     # keep text editable in svg
    mpl.rcParams['xtick.direction'], mpl.rcParams['ytick.direction'] = 'in', 'in'
    mpl.rcParams['axes.spines.top'], mpl.rcParams['axes.spines.right'] = False, False

def norm01(x):
    import numpy as np
    return (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x))

def stim_cmap(base_color, l_range=(.84, .17), s_range=(1, 1)):
    import numpy as np, colorsys, matplotlib.colors as mcolors

    # lets lightness carry stimulus while hue keeps marking context; saturation rides
    # along with lightness so the ramp's two ends separate further than either alone
    h, _, s = colorsys.rgb_to_hls(*mcolors.to_rgb(base_color))
    ls, fs  = np.linspace(*l_range, 256), np.linspace(*s_range, 256)
    return mcolors.LinearSegmentedColormap.from_list(f'{base_color}_stim', [colorsys.hls_to_rgb(h, l, s * f) for l, f in zip(ls, fs)])
