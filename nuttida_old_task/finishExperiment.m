% This function wraps up after the session, finishing the experiment

function finishExperiment(taskStruct, dispStruct)
setMarkerIDs;

% Final TTL
if ~taskStruct.debug
    sendTTL(EXPERIMENT_OFF);
end

% Wrapping up EyeLink file
if taskStruct.eyeLinkMode
    disp(['Receving file and store to:' taskStruct.edffilename ... 
        ' to ' taskStruct.edffilename_local]);
    [~] = Eyelink('ReceiveFile', taskStruct.edffilename, taskStruct.edffilename_local);
    writeLog_withEyelink(taskStruct.fidLog, EXPERIMENT_OFF,'',taskStruct.eyeLinkMode);
    Eyelink('Message', 'Regular Stop');
    Eyelink('StopRecording');
    Eyelink('CloseFile');    
    Eyelink('Shutdown'); 
end

% End of session message on screen
blockAccuracy = 100*nansum(taskStruct.correctResponses==taskStruct.respKey)/taskStruct.nTrialsPerBlock;
intermissionScreen(['End of session! \n Your final accuracy was ' num2str(blockAccuracy,3) '%'], ... 
    taskStruct, dispStruct);

return;
end