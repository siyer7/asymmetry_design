% This function just tries a bunch of serial ports onto which the 
% button box handle can be created
function handle=initCEDRUS

% handle = CedrusResponseBox('Open', 'COM7',[],0);

try
    handle = CedrusResponseBox('Open', 'COM7',[],0);
catch 
    try
        handle = CedrusResponseBox('Open', 'COM4',[],0);
    catch
    try
        handle = CedrusResponseBox('Open', 'COM3',[],0);
    catch
        try
      handle = CedrusResponseBox('Open', 'COM6',[],0); 
        catch 
            try
              handle = CedrusResponseBox('Open', 'COM5',[],0); 
            catch
                try
                     handle = CedrusResponseBox('Open', 'COM8',[],0); 
                catch
                     handle = CedrusResponseBox('Open', 'COM9',[],0); 
                end
            end
        end
    end
    end
end