%% Generate full shape spaces for ShapeDim experiment

% This code makes smaller versions of the images
% can be helpful for loading quickly

%%
clear

close all;

root = '/mnt/neurocube/local/serenceslab/maggie/shapeDim/'

addpath(genpath(fullfile(pwd,'ODBstims')));

kinds = {'balls','snm','fetus','cater'};

my_back_color = 0.3;
my_shape_color = 0.9;

%% loop over stim sets
for kk=[3]

    savedir = fullfile(root, 'Stimuli',sprintf('AmpGrid%d_adj_full_grey_small',kk));
    savedir = ('.');
    if ~isdir(savedir)
        mkdir(savedir);
    end
    kind = kinds{kk};
    mysize = 224;
    stimtype = '';
%     start = 0.2;    % min value along each axis
%     stop = 4.8;
    start=0;
    stop=5;
    step = 0.1;
    
    % how much to adjust dimension 1 by (less than 1==shrink)
    scale_by = 0.8;
 
    all_pts_dim2 = start:step:stop;  
    % to scale, going to keep the center of the axis same but multiply by
    % scale factor here
    center=mean(all_pts_dim2);
    all_pts_dim1 = scale_by*(all_pts_dim2-center)+center;
    [gridx,gridy] = meshgrid(all_pts_dim1,all_pts_dim2);
    all_grid_points = [gridx(:),gridy(:)];

    % these coordinates are just how the images will be named, to make naming
    % consistent with other version of this set...
    % but now the physical stimuli these coords correspond to
    % are slightly different. 
    [gridx_savename,gridy_savename] = meshgrid(all_pts_dim2,all_pts_dim2);
    all_grid_points_savename = [gridx_savename(:),gridy_savename(:)];
    %% now actually craete and save the shapes
    for pp=1:size(all_grid_points,1)

        X=RFCcalc(ODBparams(kind),mysize,stimtype,all_grid_points(pp,:));

        image = X{1};
        image(image==0) = my_back_color;
        image(image==1) = my_shape_color;

        
        savepath = fullfile(savedir, sprintf('Shape_%.2f_%.2f.png',all_grid_points_savename(pp,1),all_grid_points_savename(pp,2)));
        imwrite(image,savepath);

        fprintf('saved image to %s\n',savepath);

    end

end