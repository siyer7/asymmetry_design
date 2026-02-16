% add blackrock
addpath(genpath('/home/nuttidalab/Documents/OSort/osort-v4-code/code/continuous/blackrock'));

subj = '202518';

rawNSX = sprintf('../../results/%s/raw/Hub1-20251004-160040-001.ns6', subj);
outDir = sprintf('../../results/%s/osort_mat', subj, subj);

NS6 = openNSx(rawNSX,'noread');
channels = NS6.MetaTags.ChannelID;

% Create folder for output
if ~exist(outDir,'dir'); mkdir(outDir); end

% Run converter
convertNSx_toMat(rawNSX, channels, outDir, 1, 'BL');

% NS6hdr  = openNSx(rawNSX,'noread');
% chanIDs = NS6hdr.MetaTags.ChannelID(:)';
% n       = numel(chanIDs);
% 
% convertNSx_toMat(rawNSX, 1:n, outDir, 1, 'BL');   % indices 1..n, not IDs
% 
% % rename BL<index>.mat -> BL<ChannelID>.mat so filenames match your convention
% for ii = 1:n
%     movefile(fullfile(outDir, sprintf('BL%d.mat', ii)), ...
%              fullfile(outDir, sprintf('BL%d.mat', chanIDs(ii))));
% end
