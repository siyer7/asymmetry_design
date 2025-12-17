addpath(genpath('/home/nuttidalab/Documents/OSort/osort-v4-code'));
subj = '202521';

% where your BL*.mat files live (from convertNSx_toMat)
paths = struct();
% base folder = parent of mat_files/, sort/, figs/, osort_out/
paths.basePath   = sprintf('../../results/%s/osort_mat/', subj);    % must end with /
% input
paths.pathRaw    = [paths.basePath 'nsx2mat/'];
paths.pathOut    = [paths.basePath 'sorted_mats/'];
paths.pathFigs   = [paths.basePath 'figs/'];
paths.timestampspath = paths.basePath;
paths.patientID  = 'P1';   % any label for plots

if ~exist(paths.pathOut,'dir'), mkdir(paths.pathOut); end
if ~exist(paths.pathFigs,'dir'), mkdir(paths.pathFigs); end

% choose which channels (match your BL*.mat files)
matFiles = dir(fullfile(paths.pathRaw,'BL*.mat'));
filesToProcess = sort(arrayfun(@(f) sscanf(f.name,'BL%d.mat'), matFiles));

% % ADDED TO SELECT ONLY CHANNELS THAT NEED RE-SORTING.
% wantedChans    = [210];  % the BL numbers you care about
% filesToProcess = filesToProcess(ismember(filesToProcess, wantedChans));

% filesToProcess  = filesToProcess(1:5); % all would be 1:32
noiseChannels   = [];
groundChannels  = [];
doGroundNormalization = 0;
normalizeOnly   = [];

% alignment (max/min)
filesAlignMax = [];
filesAlignMin = filesToProcess;

% threshold
extractionThreshold = 5;   % adjust as needed
thres = repmat(extractionThreshold,1,length(filesToProcess));

% OSort params
paramsIn = [];
paramsIn.rawFilePrefix       = 'BL';
paramsIn.processedFilePrefix = 'A';
paramsIn.rawFileVersion      = 7;         % Matlab .mat
paramsIn.samplingFreq        = 30000;     % only used if rawFileVersion==3

% display & output
paramsIn.displayFigures      = 0;         % don't pop up windows
paramsIn.outputFormat        = 'png';

% alignment / detection
paramsIn.defaultAlignMethod  = 3;         % 1=max, 2=min, 3=mixed
paramsIn.peakAlignMethod     = 1;         % 1 find Peak, 2 none, etc.
paramsIn.detectionMethod     = 1;         % 1=power, 2=T pos, etc.
dp.kernelSize = 18; 
paramsIn.detectionParams     = dp;

% thresholds & clustering
paramsIn.thresholdMethod     = 1;         % 1=approx, 2=exact
paramsIn.prewhiten           = 0;
paramsIn.minNrSpikes         = 50;        % REQUIRED → fixes your error
paramsIn.blockNrRawFig       = [10 15 20 25];
paramsIn.doGroundNormalization = 0;
paramsIn.tillBlocks          = 999;       % how many 20s blocks to process
paramsIn.doDetection         = 1;
paramsIn.doSorting           = 1;
paramsIn.doFigures           = 1;
paramsIn.noProjectionTest    = 1;
paramsIn.doRawGraphs         = 1;


% run OSort
[normalizationChannels,paramsIn,filesToProcess] = ...
    StandaloneGUI_prepare(noiseChannels,doGroundNormalization,paramsIn, ...
                          filesToProcess,filesAlignMax,filesAlignMin, ...
                          normalizeOnly,groundChannels);

StandaloneGUI(paths, filesToProcess, thres, normalizationChannels, paramsIn);
