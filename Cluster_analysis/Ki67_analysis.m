%{
Copyright (c) 2025, Majed Kikhia
All rights reserved.

This source code is licensed under the BSD-style license found in the
LICENSE file in the root directory of this source tree. 

Author: Majed Kikhia
February 2025
%}

%% Cluster analysis with DBSCAN
% 
% Input: one spreadsheet that contain the locations, the Confetti labels, 
% the Ki-67 label of all analyzed cell.
% 
% What does the code do? 
% It extract metadata from the file name, apply DBSCAN analysis, calculate 
% some gerenal Ki-67 statistics for each animal, calculate a proliferation 
% index for each clone, and plot the data
% 
% Output: Supplementary Fig. 4b, c
%
% Prerequsits: Statistics and Machine Learning Toolbox
%
% Running-time for demo data = ~ 5 seconds

tic

clc
clear

%% Reading the data files and creating an empty output table

data = readtable(fullfile('..', 'Data', 'Ki-67.xls'), VariableNamingRule='preserve'); % relative path to 'Data/Ki-67.xls'. Change if needed.
keep_varaibles = ["X 1: Position X", "Y 1: Position Y", "Ki-67", "ID", "Original Component Name", "Original Image Name"];
data = data(:, keep_varaibles);

data = renamevars(data, ["X 1: Position X", "Y 1: Position Y", "Original Component Name",  "Original Image Name"], ["Position_X", "Position_Y", "Confetti", "img_name"]);

data = convertvars(data,@iscellstr,"string");

data.Confetti = categorical(data.Confetti); 
%% Inserting metadata from file names

% Inserting timepoint

for i = 1:height(data)
    if contains(data.img_name(i), "2d")
        data.Timepoint{i} = '2d';
    elseif contains(data.img_name(i), "1w")
        data.Timepoint{i} = '1w';
    elseif contains(data.img_name(i), "2w")
        data.Timepoint{i} = '2w';
    end
end

data.Timepoint = categorical(data.Timepoint); 

% Inserting Animal ID and converting Ki-67 to logical
for i = 1:height(data)
    data.Animal{i} = strcat('KGE-', data.img_name{i}(14:17));   
end

data.("Ki-67") = strcmp(data.("Ki-67"), "1"); % convert to logical

data.Animal = categorical(data.Animal); 
data.img_name = categorical(data.img_name); 

%rearranging the table  
data = data(:, [1 2 5 3 4 7 8 6]);

%% DBSCAN analysis

% Indexing based on image name and Confetti label  
columnsOfInterest = data(:, {'Confetti', 'img_name'});
[~, ~, uniqueIndices] = unique(columnsOfInterest, 'rows');

data.confetti_img_id = uniqueIndices; 

for i = 1:max(data.confetti_img_id)
    idx = (data.confetti_img_id == i);
    data.dbscan_id(idx) = dbscan(data{idx, (1:2)},50,1);
end

% Indexing based on clone IDs 
columnsOfInterest = data(:, {'Confetti', 'img_name', 'dbscan_id'});
[~, ~, uniqueIndices] = unique(columnsOfInterest, 'rows');

data.clone_id_dataset = uniqueIndices;

img_names = categories(data.img_name);

for i = 1:length(img_names)
    idx = (data.img_name == img_names{i});
    columnsOfInterest = data(idx, {'img_name', 'clone_id_dataset'});
    [~, ~, uniqueIndices] = unique(columnsOfInterest, 'rows');
    data.clone_id_img(idx) = uniqueIndices;
end

%% Scatter plot for an example image with Confetti colors

selected_image = '20241107_KGE-4774_L7_1w_e_Series009.ims';
idx = (data.img_name == selected_image);

colors = [0 1 1; 0 1 0; 1 0 0; 1 1 0];

scatterColors = colors(data.Confetti(idx), :);

img_dim_xy = 369.41; 

figure()
scatter(data.Position_X(idx), data.Position_Y(idx), 50, scatterColors, "filled")

title(strcat(selected_image, '_ConfettiColors'), 'Interpreter', 'none')

xlim([0 img_dim_xy])
ylim([0 img_dim_xy])
xlabel('X in microns')
ylabel('Y in microns')

%% Scatter plot for same example image with distinct color for clones

colors = hsv(max(data.clone_id_img(idx))); % 'lines', 'parula', or 'hsv' can be used for distinct colors

scatterColors = colors(data.clone_id_img(idx), :);

figure()
scatter(data.Position_X(idx), data.Position_Y(idx), 50, scatterColors, "filled")
title(strcat(selected_image, '_CloneIDs'), 'Interpreter', 'none')
xlim([0 img_dim_xy])
ylim([0 img_dim_xy])
xlabel('X in microns')
ylabel('Y in microns')

%% Calculate for each mouse the total number of Ki-67+ cells and their % of all Confetti+ cells

columnsOfInterest = data(:, {'Animal', 'Timepoint'});
[animal_summary, ~, uniqueIndices] = unique(columnsOfInterest, 'rows');

timepointOrder = {'2d', '1w', '2w'};

animal_summary.Timepoint = categorical(animal_summary.Timepoint, timepointOrder, 'Ordinal', true);

for i = 1:max(uniqueIndices)
    idx = (uniqueIndices == i);
    animal_summary.Ki67_no(i) = sum(data.("Ki-67")(idx));
    animal_summary.Ki67_freq(i) = sum(data.("Ki-67")(idx))*100/length(data.("Ki-67")(idx));  
end

figure()
scatter(animal_summary, "Timepoint", "Ki67_no", "filled", "XJitter", "density")

title('Number of Ki-67^+ Confetti^+ cells per animal')
ylabel('Total number of Ki-67^+ Confetti^+ cells per animal')

figure()
scatter(animal_summary, "Timepoint", "Ki67_freq", "filled", "XJitter", "density")

title('Frequency of Ki-67^+ of all Confetti^+ per animal')
ylim([0 100])
ylabel('Frequency of Ki-67^+ of all Confetti^+ per animal (%)')
%% Calculate Ki-67 percentage for each clone (proliferation index) 

clone_id_dataset = (1:max(data.clone_id_dataset))';
clone_summary = table(clone_id_dataset);

for i = 1:max(data.clone_id_dataset)
    idx = (data.clone_id_dataset == i);
    clone_summary.Ki67_number(i) = sum(data.("Ki-67")(idx));
    clone_summary.Confetti_number(i) = length(data.("Ki-67")(idx));
    x = data.Timepoint(idx);
    y = data.Animal(idx);
    z = data.Confetti(idx);
    clone_summary.Timepoint(i) = x(1);
    clone_summary.Animal(i) = y(1);
    clone_summary.Confetti_color(i) = z(1);
end

clone_summary.Prolifetation_index = (clone_summary.Ki67_number.*100)./clone_summary.Confetti_number;

clone_summary = clone_summary(:, [1 2 3 7 6 5 4]);
%% Plot the correlation clone size vs. proliferation index for all clones for every time point

% In these plots the size of the circles reflects the number of events
% that shared the same two parameters (clone size and proliferation
% index. 

columnsOfInterest = clone_summary(:, {'Confetti_number', 'Prolifetation_index', 'Timepoint'});
[C,ia,ic] = unique(columnsOfInterest, 'rows');

timepointOrder = {'2w', '1w', '2d'};

C.Timepoint = categorical(C.Timepoint, timepointOrder, 'Ordinal', true);

C.value_occurance = accumarray(ic,1);

idx = (C.Timepoint == '2d');

figure()
scatter(C.Confetti_number(idx), C.Prolifetation_index(idx), 10*C.value_occurance(idx), "filled", "MarkerFaceAlpha", 0.5)

title('Correlation of clone size and proliferation index at 2d timepoint')
ylim([0 100])
xlim([0 70])
xlabel('Clone size')
ylabel('Proliferation index (%)')

idx = (C.Timepoint == '1w');

figure()
scatter(C.Confetti_number(idx), C.Prolifetation_index(idx), 10*C.value_occurance(idx), "filled", "MarkerFaceAlpha", 0.5)

title('Correlation of clone size and proliferation index at 1w timepoint')
ylim([0 100])
xlim([0 70])
xlabel('Clone size')
ylabel('Proliferation index (%)')

idx = (C.Timepoint == '2w');

figure()
scatter(C.Confetti_number(idx), C.Prolifetation_index(idx), 10*C.value_occurance(idx), "filled", "MarkerFaceAlpha", 0.5)

title('Correlation of clone size and proliferation index at 2w timepoint')
ylim([0 100])
xlim([0 70])
xlabel('Clone size')
ylabel('Proliferation index (%)')

%% The same previous plots in one 3D plot where time is included as a third dimension

colors = [0.4940 0.1840 0.5560; 0.8500 0.3250 0.0980;0 0.4470 0.7410];
scatterColors = colors(C.Timepoint(:), :);

figure()
s = scatter3(C.Confetti_number(:), C.Timepoint(:), C.Prolifetation_index(:), 10*C.value_occurance(:), scatterColors, "filled", "MarkerFaceAlpha", 0.5);

title('Correlation of clone size and proliferation index at 2d, 1w, and 2w timepoint')
xlabel('Clone size (cells)')
ylabel('Timepoint')
zlabel('Prolifetation index (%)')

