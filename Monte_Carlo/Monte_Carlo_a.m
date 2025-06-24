%{
Copyright (c) 2025, Majed Kikhia
All rights reserved.

This source code is licensed under the BSD-style license found in the
LICENSE file in the root directory of this source tree. 

Author: Majed Kikhia
February 2025
%}

%% Monte Carlo Simulaiton 
% 
% Input: spreadsheets that contain the location and the Confetti labels for 
% microglial cells. Each spreadsheet represents one image.
% These spreadsheets also contain some experimental and metadata.
% 
% What does the code do? 
% This code runs Monte Carlo simulations to investigate if microglia
% undergo clonal expansion after MCAo. 
%
% Code structure: 
% Some steps of the code required several hours of compuation. Therefore,
% the code is written in sections. The results for sections with heavy 
% computation are saved in the matlab working folder, and the next section 
% checks for those results before running.The best way to run this code is
% section by section.
%
% Outputs of Monte_Carlo_a:
%   1. A Matlab table containig the sampling vectors for each mouse.
%   2. A Matlab table containig 1000 simulations for each image.
% 
% Note: running time is estimated for running the demo dataset on a loptop
% with the following configurations:
%   Processor: AMD Ryzen 7 5825U with Radeon Graphics   2.00 GHz
%   Installed RAM: 16.0 GB (15.4 GB usable)
% 
% % Running-time for demo data = ~ 70 seconds

tic

clc 
clear

%% Creating sampling vectors

% Read all data in one datastore
ds = datastore(fullfile('..', 'Data', 'Cell_locations', '*.xls'));   % relative path to 'Data/Cell_locations'. Change if needed.
ds.VariableTypes = {'double', 'double', 'double', 'categorical', ...
    'double', 'categorical', 'categorical', 'categorical', ...
    'categorical', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'categorical'}; % datastores do not accept logical variables! So Stroke here is a categorical variable.

% Read the datastore as a table
T = readall(ds);

% Exclude Iba-1 cells from the table
idx = T.Channel == 'Iba-1'; 
T = T(~idx,:);

% Indexing the data by animal and by imgage file name
[Ca, ~, ica] = unique(T.Animal);
[~, iai, ~] = unique(T.FileName);

% Collecting the sample vectors from the table
% A sample vector is a vector containing the color (channel) of all
% labelled clees in one mouse. 

Num_animals = height(Ca); 
Num_images = height(iai);
Num_cells = height(ica);

sampling_vectors = cell(Num_animals,1); 

for i = 1:Num_animals
    sampling_vectors{i,1} = T.Channel(T.Animal == Ca(i));
end 

SVT = table('Size',[Num_animals 2],'VariableTypes',{'categorical', ...
    'categorical'}); 
SVT.Properties.VariableNames = {'Animal', 'sampling_vectors'}; 

SVT.Animal = Ca;
SVT.sampling_vectors = sampling_vectors; % This table contain the name of each mouse and its sampling vector next to it. 

save('Output\SVT', 'SVT')

fprintf('Sampling-vectors were created and saved as SVT.mat\n');
%% Arranging the data required for the simulation

% Creating a table to save the results for each image together with
% metadata

fprintf('Creating the table "simulation" to save the results for each image together with the metadata\n');

simulations =  table('Size',[Num_images 14],'VariableTypes',{'cell', ...
    'cell', 'cell', 'double', 'categorical', 'categorical', ...
    'categorical', 'categorical', 'double', 'double', 'double', ...
    'double', 'double', 'categorical'});
simulations.Properties.VariableNames = {'CellLocations' 'RealColor' ...
    'Simulations' 'Num_Microfetti' 'Experiment' 'Group' 'Animal' ...
    'Stroke' 'ImageSeries' 'ImageOrder' 'ImageSizeZ' 'StepSizeZ' ...
    'VoxelSizeXY' 'FileName'};

% Copying metadata
for i = 1:Num_images   
    simulations.Experiment(i) = T.Experiment(iai(i));
    simulations.Animal(i) = T.Animal(iai(i));
    simulations.Group(i) = T.Group(iai(i));
    simulations.Stroke(i) = T.Stroke(iai(i));
    simulations.ImageSeries(i) = T.ImageSeries(iai(i));
    simulations.ImageOrder(i) = T.ImageOrder(iai(i));
    simulations.ImageSizeZ(i) = T.ImageSizeZ(iai(i));
    simulations.StepSizeZ(i) = T.StepSizeZ(iai(i));
    simulations.VoxelSizeXY(i) = T.VoxelSizeX(iai(i));
    simulations.FileName(i) = T.FileName(iai(i));
end

% Calculating the number of Microfetti positive cells per image and copying
% the locations and the color of cells to the new table
for i = 1:Num_images-1
    simulations.Num_Microfetti(i) = iai(i+1) - iai(i);
    simulations.CellLocations{i} = T{iai(i):(iai(i+1)-1), 1:3};
    simulations.RealColor{i} = T{iai(i):(iai(i+1)-1),"Channel"};
end
simulations.Num_Microfetti(end) = Num_cells - iai(end) +1;
simulations.CellLocations{end} = T{iai(end):Num_cells, 1:3};
simulations.RealColor{end} = T{iai(end):Num_cells,"Channel"};

clear T ds Ca iai ica idx % to free up memory

% From here the code uses the table "simulations" for further analysis 

% Calculating pairwise-distances between cells
for i = 1:Num_images
 simulations.pdist{i} = squareform(pdist(simulations.CellLocations{i}));
end

%% The simulation of color labelling

% The first step is to create simulated data for each image by randomely
% selecting color labels from the sampling vector of the corresponding 
% mouse equal to the total number of Confetti-postive cells in the image.
% The results are saved in the "simulation" table. 

Num_simulations = 1000;

fprintf('Creating simulated color labels from the sampling-vectors.\n This takes some time. Please be patient!\n');

for i = 1:Num_images
    for j = 1:Num_simulations
        simulations.Simulations{i}(:,j) = datasample(SVT.sampling_vectors{SVT.Animal == simulations.Animal(i)},simulations.Num_Microfetti(i));
    end
end

save('Output\simulations', 'simulations')

fprintf('The simulated color labels were created and saved in the table simulations.mat\n Please run the code "Monte_Carlo_6b.m"\n');

toc