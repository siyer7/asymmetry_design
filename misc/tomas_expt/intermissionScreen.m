% This function creates a multi-purpose intermission screen to buffer
% different stages of the task

function intermissionScreen(text, taskStruct, dispStruct)
% Displaying same message on MATLAB console
display(text);
% Cleaning screen
Screen('Flip', dispStruct.win);
% Plotting intermission text on center of screen
DrawFormattedText(dispStruct.win, text, 'center', 'center', [255,255,255]);
Screen('Flip', dispStruct.win);
% Waiting for continue key to be pressed (usually 'c')
taskStruct.continueKey;
intermissionOver = false;
while intermissionOver == false
    % Check the keyboard to see if a button has been pressed        
    [~, ~, keyCode] = KbCheck;
    if keyCode(taskStruct.continueKey)
        intermissionOver = true;
    end
end
end