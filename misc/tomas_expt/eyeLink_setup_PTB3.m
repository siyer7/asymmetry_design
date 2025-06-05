function [retCodeSuccess,el]=eyeLink_setup_PTB3(window,dummymode, timestampStr)
retCodeSuccess=1; % 1 if succeeded, 0 if failure

%=== adapted from EyelinkImageExample.m demo file

% STEP 2
% Provide Eyelink with details about the graphics environment
% and perform some initializations. The information is returned
% in a structure that also contains useful defaults
% and control codes (e.g. tracker state bit and Eyelink key values).
el=EyelinkInitDefaults(window);

el.targetbeep=1;

if ~EyelinkInit(dummymode)
    fprintf('Eyelink Init aborted.\n');
    cleanup;  % cleanup function
    return;
end

[v vs]=Eyelink('GetTrackerVersion');
fprintf('Running experiment on a ''%s'' tracker.\n', vs );

% open file for recording data
edfFile=[timestampStr '.edf'];

Eyelink('Openfile', edfFile);

disp(['EDF filename: ' edfFile] );

% STEP 4
% Do setup and calibrate the eye tracker

disp('Running EyeLink setup...');
el.backgroundcolour = 255;
EyelinkUpdateDefaults(el);
ret1 = EyelinkDoTrackerSetup(el);


if ret1==0
    % do a final check of calibration using driftcorrection
    % You have to hit esc before return.
    disp('Running EyeLink drift corr...');
    ret2=EyelinkDoDriftCorrection(el);

% STEP 5
    % Start recording eye position
    Eyelink('StartRecording');
    % record a few samples before we actually start displaying
    WaitSecs(0.5);
    
    retCodeSuccess=1;
else
    
    retCodeSuccess=0;
    
end

