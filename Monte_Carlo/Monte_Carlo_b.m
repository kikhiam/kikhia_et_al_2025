%{
Copyright (c) 2025, Majed Kikhia
All rights reserved.

This source code is licensed under the BSD-style license found in the
LICENSE file in the root directory of this source tree. 

Author: Majed Kikhia
February 2025
%}

%% Calculating the ring volumes for each cell 
% 
% Input: The Matlab table "simulations.mat", which is the ouput of
% Monte_Carlo_a
%
% What does the code do? 
% It calculates the volumes of the concentric rings around each cell.
% The process is computationally heavy because for rings that are partially
% residing in the image, the volume outside the image should be excluded. 
% This is achieved using the Antenna Toolbox which allows to work with
% geometric shapes and caclulate thier intersection areas.
%
% Code structure: 
% The code uses parfor for executing for-loop iterations in parallel 
% on workers, which makes the code 4x faster. The table simulation is a 
% braodcast variable in this code and matlab gives a warning for it but
% code still runs and the warning can be ignored. 
%
% Output: N x 1 cell, where N is the number of images. Each row contains a
% matrix n x 30 matrix where n is the number of cells in that image. These
% matrices save the volume of each ring around each cell in the dataset. 
% 
% Prerequsits: Antenna Toolbox™, Parallel Computing Toolbox
%
% Runtime on demo data = ~ 40 minutes. 

tic

clc
clear

load("Output\simulations.mat")

simulations = simulations(:,["CellLocations", "Num_Microfetti", ...
    "ImageSizeZ", "StepSizeZ", "VoxelSizeXY", "FileName"]);

r0 = 10:10:300; %Radii of the ring 
r1 = r0 - 10; %Radii of the inner circle 
r2 = r0 + 10; %Radii of the outer circle

ring_volumes = cell(height(simulations), 1);


parfor i = 1:height(simulations)
    ring_volumes{i,1} = rings_antenna(simulations.CellLocations{i}, ...
        simulations.Num_Microfetti(i), 1024*simulations.VoxelSizeXY(i)...
        *10^6, 1024*simulations.VoxelSizeXY(i)*10^6, ...
        simulations.StepSizeZ(i)*simulations.ImageSizeZ(i)*10^6, r1, r2);
    fprintf('Image %i was analyzed\n', i);
end


save('Output\ring_volumes', 'ring_volumes')

toc
