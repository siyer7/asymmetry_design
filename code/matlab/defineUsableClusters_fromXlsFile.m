%
% takes sorting result from excel file and defines usable clusters automatically
%
%
% requires three columns in the excel file: channelNr, threshold, clustersToUse
% channelNr and threshold have to be type numeric
% clustersToUse has to be type text. it is a list of clusters to use, with merges indicated by '+'
% Examples: 128,130,150     (use 3 clusters)
%           128+130,150     (use 2 clusters, merge 128+130)
%
%
% 
%urut/april14

%% ==== define all parameters to process a new session
clear; clc;

%% Excel file containing which clusters to analyze
xls_filename = 'sorting_template.xlsx';

%% Directory (MAKE SURE BASEPATH IS SET CORRECTLY!)
basepath = '/home/nuttidalab/Documents/asymmetry_design/results/202512/osort_mat';

sortDir='sort';
figsDir='figs';
finalDir='final';
Area = 'A'; % A for AIP, B for BA5

% - Patient Parameters
taskStr = fullfile(basepath);

% - Excel Parameters
xlsFile = fullfile(taskStr, xls_filename); %path to XLS file 
sheet = 'Sheet1';  % which worksheet has this session
range=[11:130]; %#ok<NBRAK> % which rows in this worksheet contain all channels to be used (11:90 = Cedars) 
columnChannel=3;
columnThresh=5; 
columnClusterNumbers=6;

% NOTE: code below will reference the  following paths
% basepath\sortDir --> current location of cells and sorted_new mat-files for all clusters
% basepath\figsDir --> current location of png files for all clusters
% basepath\sortDir\final --> where mat files for selected clusters ONLY will be stored
% basepath\figsDir\final --> where png files for selected clusters ONLY will be stored

%% ===
[num,txt,raw] = xlsread(xlsFile, sheet, '','basic');  % range selection does not work in basic mode (unix)

masterTable=[]; %List of all final clusters: Channel Th Cluster

masterMerges=[]; %list all merges: Channel Th MergeTarget MergeSource

for k=range
    disp(k)
    channelNr = (raw{k,columnChannel});
    thresh = raw{k,columnThresh};
    clusters = raw{k,columnClusterNumbers};

    if isnumeric(clusters)
        clusters=num2str(clusters);
    end
    
    if ~isempty(thresh) && ~isnan(thresh) 
        
        %define clusters
        cls = strsplit(clusters, ',');
        clsToUse=[];
        for jj=1:length(cls)
            
            %see if this is a merge operation
            indsPlus = strfind(cls{jj},'+');
            if ~isempty(indsPlus)
                %assume merge was done already, only use first entry (main cl)
                clToUse=str2num( cls{jj}(1:indsPlus-1) );
                
                mergeSources = cls{jj}(indsPlus+1:end);
                
                mergeSources_list = strsplit(mergeSources,'+');
                
                for i=1:length(mergeSources_list)
                    masterMerges = [ masterMerges; [channelNr thresh clToUse str2num(mergeSources_list{i}) ]];
                end
            else
                clToUse = str2num(cls{jj});
            end
            if clToUse>0
                clsToUse = [clsToUse clToUse];
                masterTable = [masterTable; [ channelNr thresh clToUse ]];
            end
        end
        
        disp([num2str(k) ' Using: Ch=' num2str(channelNr) ' Th=' num2str(thresh) ' ClsOrig:' clusters]);
        disp(['ClsParsed:' num2str(clsToUse)]);
    else
        thresh=5;
    end
end

%% execute merges
for j=1:size(masterMerges,1)
    overwriteParams = [ masterMerges(j,:) ];
    % basePath, sortVersion, figsVersion, finalVersion, overwriteParams 
    mergeClusters( taskStr, sortDir, [figsDir num2str(thresh)], finalDir, overwriteParams );
end


%% execute define usable clusters (Error prone due to manual data entry)

% Add section to automatically resume where the error occured if you
% re-read the excel without clearing variables. 
if ~exist('errChan','var') 
    startChan = 1; % Start from beginning
else
    startChan = errChan; % Start from error
end

channelsToProcess = unique( masterTable(:,1) );
for ii = startChan:length(channelsToProcess) % Recommended to not use parfor, for easier troubleshooting
    
    channelNr = channelsToProcess(ii);    
    inds = find(masterTable(:,1)==channelNr);
    
    cls = masterTable(inds,3);
    thresh = masterTable(inds,2);  %all thresholds are the same
    thresh = thresh(1);
    
    overwriteParams = [ channelNr thresh unique(cls') ];
    % basePath, sortVersion, figsVersion, finalVersion, overwriteParams
    if ~defineUsableClusters( taskStr, sortDir, [figsDir num2str(thresh)], finalDir, overwriteParams );
%    if ~defineUsableClusters_v2( [basepath patientID],Area, sortDir, figsDir, finalDir, overwriteParams );
        error('error, fix manually and repeat: Ch%d (%d)',channelNr,j);
        errChan = channelNr; %#ok<UNRCH>
    end
    
end

%% next processing steps
overwriteWarningDisable=1;
rangeToRun=[1:300]; %for Cedars (129:208)
% convertClustersToCells(basepath, patientID, ['/' sortDir '/' finalDir '/'], rangeToRun,overwriteWarningDisable);
convertClustersToCells_v2( Area, basepath, fullfile(sortDir, finalDir), rangeToRun,overwriteWarningDisable)


%%%%%%%%%%%%%%%%%%%%%% STOP HERE - RK %%%%%%%%%%%%%%%%%%%%


