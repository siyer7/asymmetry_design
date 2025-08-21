% Add NPMK toolbox to the path
addpath(genpath('NPMK'));

% Base name of your Blackrock recording
base = '../results/202512/datafile_202512b002_DI';

%% --- Continuous (NSx) ---
% Read .ns6 and scale to µV automatically
nsx = openNSx([base '.ns6'], 'read', 'uV');

data_uV = nsx.Data;                          % [nChan x nSamples], µV, double
Fs      = nsx.MetaTags.SamplingFreq;         % sampling rate in Hz
t0      = double(nsx.MetaTags.Timestamp);    % starting sample
nSamp   = double(nsx.MetaTags.DataPoints);   % total samples
t       = ( (0:nSamp-1) + t0 ) / Fs;         % time vector in seconds

labels  = string({nsx.ElectrodesInfo.Label}).';   % channel labels

% Save as MAT v7.3 (needed for large arrays)
save('../data/neur/pt3_nsx_uV.mat', 'data_uV', 't', 'Fs', 'labels', '-v7.3');

clear nsx

%% --- Spikes/Events (NEV) ---
nev = openNEV([base '.nev'], 'read', 'report');

spkTimes_s = double(nev.Data.Spikes.TimeStamp) / double(nev.MetaTags.TimeRes);
spkChan    = double(nev.Data.Spikes.Electrode);
spkUnit    = double(nev.Data.Spikes.Unit);

save('../data/neur/pt3_nev_spikes.mat', 'spkTimes_s', 'spkChan', 'spkUnit', '-v7.3');
