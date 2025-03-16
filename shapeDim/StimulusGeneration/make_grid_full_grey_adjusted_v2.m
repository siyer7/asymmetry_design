% This script will save out a larger grid of shapes, to allow for
% parametric manipulation of the difficulty. 
% Saves about 2500 images

% This version allows for background color to be chosen - good for the
% scanner because black shows dust from projector (trying to match these
% parameters closely for behav experiment)

% This set is what we used for final fMRI experiment
% In this version there is a scaling factor of 0.80 applied to dimension 1, 
% this is because in an early pilot version we found that axis 1 was 
% slightly easier than axis 2, so the scaling is meant to equate them. 
% This turns out to be a slight overcorrection.

%%

clear
close all;

root = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/';
root = '/Users/f0064z8/Library/CloudStorage/GoogleDrive-si2442@columbia.edu/My Drive/research/asymmetry_design/input_data/stims'
addpath(genpath(fullfile(pwd,'ODBstims')));

kinds = {'balls','snm','fetus','cater'};

my_back_color = 0.5; % 0.3
my_shape_color = 0.9;

%% Fixed y-coordinate and x-values
fixed_y = 2.00;  % Fixed y-coordinate for all images
x_values = 0:0.1:4;  % x-values from 0.00 to 4.00 in steps of 0.10

%% Loop over stim sets
for kk = [3]

    savedir = fullfile(root);
    if ~isdir(savedir)
        mkdir(savedir);
    end
    kind = kinds{kk};
    mysize = 3000;
    stimtype = '';
    start = 0;
    stop = 5;
    step = 0.1;
    
    % How much to adjust dimension 1 by (less than 1 == shrink)
    scale_by = 0.8;
 
    all_pts_dim2 = start:step:stop;  
    % To scale, going to keep the center of the axis same but multiply by
    % scale factor here
    center = mean(all_pts_dim2);
    all_pts_dim1 = scale_by * (all_pts_dim2 - center) + center;
    [gridx, gridy] = meshgrid(all_pts_dim1, all_pts_dim2);
    all_grid_points = [gridx(:), gridy(:)];

    % These coordinates are just how the images will be named, to make naming
    % consistent with other versions of this set...
    [gridx_savename, gridy_savename] = meshgrid(all_pts_dim2, all_pts_dim2);
    all_grid_points_savename = [gridx_savename(:), gridy_savename(:)];
    
    %% Now create and save the shapes
    for pp = 1:length(x_values)
        % Update x and fixed y for each iteration
        x = x_values(pp);

        % Generate the shape with the fixed y-coordinate and varying x-coordinate
        X = RFCcalc(ODBparams(kind), mysize, stimtype, [x, fixed_y]);

        % Process the image
        image = X{1};
        image(image == 0) = my_back_color;
        image(image == 1) = my_shape_color;

        % Save the image with the correct naming convention
        savepath = fullfile(savedir, sprintf('Shape_%.2f_%.2f.png', x, fixed_y));
        imwrite(image, savepath);

        fprintf('Saved image to %s\n', savepath);
    end

end
