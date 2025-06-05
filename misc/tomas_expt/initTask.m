% This function sets up screen display and task structs which will contain
% important task information

function [taskStruct, dispStruct] = initTask()
%% Some initial global setup
RandStream.setGlobalStream(RandStream('mt19937ar','Seed','shuffle'));
KbName('UnifyKeyNames');

% Subejct ID info
taskStruct.subID = input('Participant number (PxxCS):\n','s');
% % Use eyetracking with Eyelink?
taskStruct.eyeLinkMode = str2double(input('Use Eyelink? 0=no, 1=yes:\n','s'));
% % Use CEDRUS (button box)?
taskStruct.useCEDRUS = str2double(input('Use CEDRUS? 0=no, 1=yes:\n','s'));
% Debug mode? Debug mode means no TTLs are sent, PTB screens are smaller
taskStruct.debug = str2double(input('Debug mode? 0=no, 1=yes:\n','s'));


% taskStruct.subID = '0';
% % Use eyetracking with Eyelink?
% taskStruct.eyeLinkMode = 0;
% % Use CEDRUS (button box)?
% taskStruct.useCEDRUS = 0;
% % Debug mode? Debug mode means no TTLs are sent, PTB screens are smaller
% taskStruct.debug = 1;


% Output folder
taskStruct.outputFolder = fullfile('..', 'patientData','taskLogs');
% Check to see if the output folder exists
if exist(taskStruct.outputFolder, 'dir') == 0
    % Folder does not exist - create it
    mkdir( taskStruct.outputFolder );
end
taskStruct.fileName = [taskStruct.subID '_Sub_' datestr(now, 'mm-dd-yyyy_HH-MM-SS')];

if taskStruct.debug
    Screen('Preference', 'SkipSyncTests', 1); % For debugging
end

%% Setting up task variables
taskStruct.nBlocks = 4;
taskStruct.nTrialsPerBlock = 48; % 48
taskStruct.nTrials = taskStruct.nTrialsPerBlock*taskStruct.nBlocks;
% Trial conditions: 1 = instruction first; 2 = instruction in the middle
% In the current implementation, conditions are blocked; change this vector
% to get different condition ordering
blockOrder_flag = randi(2);
if blockOrder_flag == 1
    taskStruct.trialConditions = [ones(taskStruct.nTrialsPerBlock,1); ...
        2.*ones(taskStruct.nTrialsPerBlock,1); ...
        ones(taskStruct.nTrialsPerBlock,1); ...
        2.*ones(taskStruct.nTrialsPerBlock,1)];
else
    taskStruct.trialConditions = [2.*ones(taskStruct.nTrialsPerBlock,1); ...
        ones(taskStruct.nTrialsPerBlock,1); ...
        2.*ones(taskStruct.nTrialsPerBlock,1); ...
        ones(taskStruct.nTrialsPerBlock,1)];
end

% Relevant axis of each trial (length, color, shape, etc.)
taskStruct.categoryNames = {'Animals','Cars','Faces','Fruits'};
taskStruct.axisNames = {{'Colorful','Count'},{'New','Colorful'},...
    {'New','Identical'},{'Count','Identical'}};
taskStruct.categoryAndAxis = [taskStruct.categoryNames;taskStruct.axisNames];
% Which of the two axes belonging to each category will be used in each
% trial (each category defines a different axis pair)
taskStruct.trialCategories = [1 2 3 4];
taskStruct.trialAxis = [1 2];
taskStruct.antiTask = [0 1];
taskStruct.promptVariant = [0 1];
taskStruct.equivalentVariantID = [0 1];
%% Guaranteeing a balanced distribution of each category x axis combination
[x1,x2,x3,x4,x5] = ndgrid(taskStruct.trialCategories,taskStruct.trialAxis,...
    taskStruct.antiTask,taskStruct.promptVariant,taskStruct.equivalentVariantID);
result = [x1(:) x2(:) x3(:) x4(:) x5(:)];
n = size(result,1);
factor = taskStruct.nTrials/n;
result = repmat(result,factor,1);
% Randomizing trial order
result = result(randperm(size(result,1)),:);
taskStruct.trialCategories = result(:,1);
taskStruct.trialAxis = result(:,2);
taskStruct.antiTask = result(:,3);
taskStruct.promptVariant = result(:,4);
taskStruct.equivalentVariantID = result(:,5);

% Determining which stimuli to use in each trial
taskStruct.stimFolder = fullfile('..', 'stimuli','Task_Stim_Version2');
taskStruct.trialStims = cell(taskStruct.nTrials,2);
taskStruct.trialPairs = zeros(taskStruct.nTrials,1);
taskStruct.stim1Position = nan(taskStruct.nTrials,1);
taskStruct.stim2Position = nan(taskStruct.nTrials,1);
taskStruct.breakTrial = zeros(taskStruct.nTrials,1); % Which trials have a break at the end
% Filling trial cells with pair numbers to be drawn from (without
% replacement) (2 axes x 4 categories x 6 pairs (1 to 6) x 4 repetitions)
taskStruct.pairNumbers = cell(2,length(taskStruct.axisNames));
for i = 1:size(taskStruct.pairNumbers,1)
    for j = 1:size(taskStruct.pairNumbers,2)
        taskStruct.pairNumbers{i,j} = repmat(1:6,1,4);
    end
end
% Filling identical cells with identical trials to be drawn from without
% replacement)
taskStruct.identicalTrials = cell(1,length(taskStruct.axisNames));
for i = 1:length(taskStruct.axisNames)
   if any(~cellfun('isempty',strfind(taskStruct.axisNames{i},'Identical')))
       idx = [zeros(1,taskStruct.nTrialsPerBlock/4) ones(1,taskStruct.nTrialsPerBlock/4)];       
       taskStruct.identicalTrials{i} = idx(randperm(length(idx)));
   end
end
identicalTrialsForReplacement = taskStruct.identicalTrials;
% Response prompts
taskStruct.promptTypes = randi(2,taskStruct.nTrials,1);
taskStruct.leftText = cell(taskStruct.nTrials,1);
taskStruct.rightText = cell(taskStruct.nTrials,1);

%% Instruction setup
% taskStruct.trialInstructions = cell(taskStruct.nTrials,2);
% taskStruct.antiTask = [0 1];
% taskStruct.promptVariant = [0 1];
% taskStruct.equivalentVariantID = [0 1];
% % Guaranteeing a balanced distribution of each task x variant x instruction
% [X,Y] = meshgrid(taskStruct.antiTask,taskStruct.promptVariant);
% result = [X(:) Y(:)];
% n = size(result,1);
% factor = taskStruct.nTrials/n;
% result = repmat(result,factor,1);
% % Randomizing trial order
% result = result(randperm(size(result,1)),:);
% taskStruct.antiTask = result(:,1);
% taskStruct.promptVariant = result(:,2);

%% Trial loop
for tI = 1:taskStruct.nTrials        
    % Selecting one pair number from pair cell and popping value (to
    % enforce no replacement)
    axis = taskStruct.trialAxis(tI);
    category = taskStruct.trialCategories(tI);        
    % Fixing sampling for 1 element (because randsample behaves differently
    % for scalars and vectors...)
    sampledVector = taskStruct.pairNumbers{axis,category};
    if length(sampledVector) == 1
        sampledVector = [sampledVector sampledVector]; %#ok<AGROW>
    end    
    trialPair = randsample(sampledVector,1);
    taskStruct.trialPairs(tI) = trialPair;
    removeIdx = find(taskStruct.pairNumbers{axis,category}==trialPair);
    taskStruct.pairNumbers{axis,category}(removeIdx(1)) = [];
    % Determining if this is an identical stimulus trial in the identical
    % axis
    if strcmp(taskStruct.axisNames{category}{axis},'Identical')
        identicalTrial = identicalTrialsForReplacement{category}(end);
        identicalTrialsForReplacement{category}(end) = [];
    end
    % Loading stimuli
    trialFolder = fullfile(taskStruct.stimFolder,[taskStruct.categoryNames{category}],['Pair',num2str(trialPair)]);
    folderImages = dir(fullfile(trialFolder,'*.jpg'));
    % Sampling 2 random images from folder (without replacement)
    sampledImages = datasample(1:length(folderImages),2,'Replace',false);
    taskStruct.trialStims{tI,1} = fullfile(trialFolder,folderImages(sampledImages(1)).name);
    taskStruct.trialStims{tI,2} = fullfile(trialFolder,folderImages(sampledImages(2)).name);
    % Keeping first image equal to second in identical trials
    if strcmp(taskStruct.axisNames{category}{axis},'Identical') && identicalTrial
        taskStruct.trialStims{tI,2} = fullfile(trialFolder,folderImages(sampledImages(1)).name);
    end
    % Randomly select position of stimuli (1 = left / 2 = right / 3 = center)
    taskStruct.stim1Position(tI) = 3;
    taskStruct.stim2Position(tI) = 3;
    if mod(tI,taskStruct.nTrialsPerBlock) == 0 && tI < taskStruct.nTrials
        taskStruct.breakTrial(tI) = 1;
    end
        
    %% Instructions
    trialAxisName = taskStruct.axisNames{category}{axis};    
    taskStruct.trialInstructions{tI,1} = getInstructionText(category,...
        trialAxisName,taskStruct.antiTask(tI),taskStruct.promptVariant(tI),...
        taskStruct.equivalentVariantID(tI));
    %% Response prompts
    if ~strcmp(trialAxisName,'Identical')
        if taskStruct.promptTypes(tI) == 1
            taskStruct.leftText{tI} = 'First';
            taskStruct.rightText{tI} = 'Second';
        else
            taskStruct.leftText{tI} = 'Second';
            taskStruct.rightText{tI} = 'First';
        end
    else
        if taskStruct.promptTypes(tI) == 1
            taskStruct.leftText{tI} = 'Yes';
            taskStruct.rightText{tI} = 'No';
        else
            taskStruct.leftText{tI} = 'No';
            taskStruct.rightText{tI} = 'Yes';
        end
    end
end
% Given the trial conditions and prompts, get the expected correct response
% for each trial in advance (to compare with actual button presses later)
taskStruct.correctResponses = getCorrectResponses(taskStruct);

% Trial event timing 
% fixation cross (0.9-1s), instruction ([1.5-3]s), stim A (1.5s), stim B (1.5s), response (up to 2s) 
taskStruct.fixationTime = 1;
taskStruct.instructionTimeMin = 2.5; % Instructions have a min/max time window
taskStruct.instructionTimeMax = 4; % Instructions have a min/max time window
taskStruct.stim1Time = 1; % Stimulus presentation time
taskStruct.ISI = 0.8;
taskStruct.stim2Time = 1; % Stimulus presentation time
taskStruct.responseTimeMax = 3; % Trial is considered missed after this time
taskStruct.textHoldoutTime = 0.5; % Time to hold text on screen after button press
taskStruct.ITI = 1; % Intertrial interval
taskStruct.instructionTime = nan(taskStruct.nTrials,1); % Initializing empty timing vectors
taskStruct.responseTime = nan(taskStruct.nTrials,1); % RT in each trial
taskStruct.trialTime = nan(taskStruct.nTrials,1); % Length of each trial for control
% Initializing response vectors to nan (will be replaced if a response is
% produced, left as nans if not)
taskStruct.respKey = nan(taskStruct.nTrials, 1);

%% Setting up Psychtoolbox displays and inputs
taskStruct.pulseKeyCode = [KbName('5'), KbName('5%')];

% Creating CEDRUS button box handle if enabled
if taskStruct.useCEDRUS
    % Making sure that previously opened serial objects are removed
    CedrusResponseBox('CloseAll');
    taskStruct.handle=initCEDRUS;
    taskStruct.leftKey = 4;
    taskStruct.rightKey = 5;
    taskStruct.confirmKey = 3;
else % Regular keyboard is being used
    taskStruct.leftKey = KbName('LeftArrow');
    taskStruct.rightKey = KbName('RightArrow');
    taskStruct.confirmKey = KbName('Space');
end
taskStruct.escapeKey = KbName('q'); % Quit the task
taskStruct.pauseKey = KbName('p'); % Pause the trial
taskStruct.continueKey = KbName('c'); % Continue after pausing

%% Creating display struct
dispStruct = struct();
% Getting screen number; by selecting "max", we ensure stimuli are
% displayed on the intended external screen
screens = Screen('Screens');
screenNumber = max(screens);
dispStruct.screenNumber = screenNumber;
if taskStruct.debug
    screenToUse = [0,0,800,600]; % Smaller screen for debugging
else
    screenToUse = []; % Actual screen for task
end
% Defining the gray background
gray = 80;
% Opening window
[dispStruct.win, dispStruct.winRect] = Screen('OpenWindow', screenNumber, gray, screenToUse);
% Making sure that OpenGL is functional in this system
AssertOpenGL;
% Make more tolerant for noisy systems
[~] = Screen('Preference','SyncTestSettings', 0.005, 50, 0.3, 5);
% Getting screen center and dimensions
dispStruct.centerX = round(dispStruct.winRect(3)/2);
dispStruct.centerY = round(dispStruct.winRect(4)/2);
[dispStruct.width, dispStruct.height] = WindowSize(dispStruct.win);
width = dispStruct.width; height = dispStruct.height;
dx = width/12;
dy = height/5;
% Size of figures on screen, in pixels (4:3)
dispStruct.stimSize = 250;
dispStruct.rewWidth = 466;
dispStruct.rewHeight = 350;
ph = dispStruct.stimSize;
pw = dispStruct.stimSize;
rw = dispStruct.rewWidth;
rh = dispStruct.rewHeight;
% Setting up stimuli presentation rectangles
% center top (1)
verticalRects(1,:) = [width/2-pw/2, height/2-dy-ph/2, width/2+pw/2, height/2-dy+ph/2];        
% center bottom (2)
verticalRects(2,:) = [width/2-pw/2, height/2+dy-ph/2, width/2+pw/2, height/2+dy+ph/2];
% left center (1)
horizontalRects(1,:) = [width/2-dx-pw/2, height/2-ph/2, width/2-dx+pw/2, height/2+ph/2];
% right center (2)
horizontalRects(2,:) = [width/2+dx-pw/2, height/2-ph/2, width/2+dx+pw/2, height/2+ph/2];
% center (rewardRect, enforcing 4:3 proportion)
horizontalRects(3,:) = [width/2-rw/2, height/2-rh/2, width/2+rw/2, height/2+rh/2];
rewardSourceRect = [160, 0, 1120, 720]; % Portion of reward videos displayed
dispStruct.rewardRect = horizontalRects(3,:);
dispStruct.rewardSourceRect = rewardSourceRect;
dispStruct.verticalRects = verticalRects;
dispStruct.horizontalRects = horizontalRects;