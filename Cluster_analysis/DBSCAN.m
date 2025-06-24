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
% Input: spreadsheets that contain the location and the Confetti labels for 
% microglial cells. Each spreadsheet represents one image.
% These spreadsheets also contain some experimental and metadata.
% 
% What does the code do? 
% It applies DBSCAN for each image and calculate the number of clones and
% the clone size. 
%
% Output: a spreadsheet containing the number of Confetti+ cells, the
% number of clones, and the clone size averaged for each mouse and brain
% hemisphere. 
% It also plot histograms of clone sizes per hemisphere for 
% every animal. Here it is a one aninal per timepoint version of 
% (Supplementary Fig. 2 and Supplementary Fig. 3).  
%
% Prerequsits: Statistics and Machine Learning Toolbox
%
% Running-time for demo data = ~ 20 seconds

tic

clc
clear

%% Reading the data files and creating an empty output table
D = dir(fullfile('..', 'Data', 'Cell_locations', '*.xls')); % relative path to 'Data/Cell_locations'. Change if needed.
filenames = {D(:).name}.';

Output = table('Size',[length(filenames) 36],'VariableTypes',{'categorical', 'categorical', 'categorical', 'categorical', 'logical', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double','double', 'cell', 'cell', 'cell', 'cell'});
Output.Properties.VariableNames = {'FileName', 'Experiment', 'Group', 'Animal', 'Stroke', 'ImageSizeZ', 'StepSizeZ', 'VoxelSizeXY', 'ImageVolumeMM3', 'NumCell_Iba1', 'Density_Iba1', 'NumMicrofetti', 'Density_Microfetti', 'NumTotal_Cluster', 'Density_Total_Cluster', 'AvgSizeClust_Total', 'NumCell_CFP', 'NumSingles_CFP', 'NumClust_CFP', 'AvgSizeClust_CFP', 'NumCell_GFP', 'NumSingles_GFP', 'NumClust_GFP', 'AvgSizeClust_GFP', 'NumCell_YFP', 'NumSingles_YFP', 'NumClust_YFP', 'AvgSizeClust_YFP', 'NumCell_RFP', 'NumSingles_RFP', 'NumClust_RFP', 'AvgSizeClust_RFP', 'SizeClust_CFP', 'SizeClust_GFP', 'SizeClust_YFP', 'SizeClust_RFP'};

%% Applying the DBSCAN for each file (image) and save the results to the output table
for k = 1:length(filenames)
    
    fullname = fullfile('..', 'Data', 'Cell_locations', D(k).name); % relative path to 'Data/Cell_locations'. Change if needed.

    data = readtable(fullname);
    data = convertvars(data,{'Channel', 'Experiment', 'Group', 'Animal', 'ImageOrder', 'FileName'}, 'categorical'); 
    
    idx = data.Channel == 'CFP';
    Position_CFP = data{idx, 1:3};
    
    idx = data.Channel == 'GFP';
    Position_GFP = data{idx, 1:3};
    
    idx = data.Channel == 'YFP';
    Position_YFP = data{idx, 1:3};
    
    idx = data.Channel == 'RFP';
    Position_RFP = data{idx, 1:3};
    
    idx = data.Channel == 'Iba-1';
    Position_Iba1 = data{idx, 1:3};
    
    if isempty(Position_CFP) == 1
        cell_in_cluster_CFP = 0;
        cluster_CFP_idx = 0;
    else
        cluster_CFP_idx = dbscan(Position_CFP,50,2);
        if sum(cluster_CFP_idx) == -length(cluster_CFP_idx)
            cell_in_cluster_CFP = length(cluster_CFP_idx) ;
        else
            cell_in_cluster_CFP = zeros(1,max(cluster_CFP_idx)+1);
            for ii = 1:max(cluster_CFP_idx)
                cell_in_cluster_CFP(1, 1) = sum(cluster_CFP_idx == -1);
                cell_in_cluster_CFP(1, ii+1) = sum(cluster_CFP_idx == ii);
            end
        end
    end
    
    if isempty(Position_GFP) == 1
        cell_in_cluster_GFP = 0;
        cluster_GFP_idx = 0;
    else
        cluster_GFP_idx = dbscan(Position_GFP,50,2);
        if sum(cluster_GFP_idx) == -length(cluster_GFP_idx)
            cell_in_cluster_GFP = length(cluster_GFP_idx) ;
        else
            cell_in_cluster_GFP = zeros(1,max(cluster_GFP_idx)+1);
            for ii = 1:max(cluster_GFP_idx)
                cell_in_cluster_GFP(1, 1) = sum(cluster_GFP_idx == -1);
                cell_in_cluster_GFP(1, ii+1) = sum(cluster_GFP_idx == ii);
            end
        end
    end
    
    if isempty(Position_YFP) == 1
        cell_in_cluster_YFP = 0;
        cluster_YFP_idx = 0;
    else
        cluster_YFP_idx = dbscan(Position_YFP,50,2);
        if sum(cluster_YFP_idx) == -length(cluster_YFP_idx)
            cell_in_cluster_YFP = length(cluster_YFP_idx) ;
        else
            cell_in_cluster_YFP = zeros(1,max(cluster_YFP_idx)+1);
            for ii = 1:max(cluster_YFP_idx)
                cell_in_cluster_YFP(1, 1) = sum(cluster_YFP_idx == -1);
                cell_in_cluster_YFP(1, ii+1) = sum(cluster_YFP_idx == ii);
            end
        end
    end
    
    if isempty(Position_RFP) == 1
        cell_in_cluster_RFP = 0;
        cluster_RFP_idx = 0;
    else
        cluster_RFP_idx = dbscan(Position_RFP,50,2);
        if sum(cluster_RFP_idx) == -length(cluster_RFP_idx)
            cell_in_cluster_RFP = length(cluster_RFP_idx) ;
        else
            cell_in_cluster_RFP = zeros(1,max(cluster_RFP_idx)+1);
            for ii = 1:max(cluster_RFP_idx)
                cell_in_cluster_RFP(1, 1) = sum(cluster_RFP_idx == -1);
                cell_in_cluster_RFP(1, ii+1) = sum(cluster_RFP_idx == ii);
            end
        end
    end
    
    Output.NumCell_Iba1(k) = size(Position_Iba1, 1);
    Output.NumCell_CFP(k) = size(Position_CFP, 1);
    Output.NumCell_GFP(k) = size(Position_GFP, 1);
    Output.NumCell_YFP(k) = size(Position_YFP, 1);
    Output.NumCell_RFP(k) = size(Position_RFP, 1);
    Output.NumMicrofetti(k) = Output.NumCell_CFP(k) + Output.NumCell_GFP(k) + Output.NumCell_YFP(k) + Output.NumCell_RFP(k);
    
    Output.NumSingles_CFP(k) = cell_in_cluster_CFP(1,1);
    Output.NumSingles_GFP(k) = cell_in_cluster_GFP(1,1);
    Output.NumSingles_YFP(k) = cell_in_cluster_YFP(1,1);
    Output.NumSingles_RFP(k) = cell_in_cluster_RFP(1,1);
    
    Output.NumClust_CFP(k) = max(cluster_CFP_idx);
    Output.NumClust_GFP(k) = max(cluster_GFP_idx);
    Output.NumClust_YFP(k) = max(cluster_YFP_idx);
    Output.NumClust_RFP(k) = max(cluster_RFP_idx); 
    
    Output.AvgSizeClust_CFP(k) = mean(cell_in_cluster_CFP(2:end));
    Output.AvgSizeClust_GFP(k) = mean(cell_in_cluster_GFP(2:end));
    Output.AvgSizeClust_YFP(k) = mean(cell_in_cluster_YFP(2:end));
    Output.AvgSizeClust_RFP(k) = mean(cell_in_cluster_RFP(2:end));
    
    
    Output.SizeClust_CFP{k} = sort(cell_in_cluster_CFP(2:end));
    Output.SizeClust_GFP{k} = sort(cell_in_cluster_GFP(2:end));
    Output.SizeClust_YFP{k} = sort(cell_in_cluster_YFP(2:end));
    Output.SizeClust_RFP{k} = sort(cell_in_cluster_RFP(2:end));

    Output.AvgSizeClust_Total(k) = mean(cat(2, Output.SizeClust_CFP{13}, Output.SizeClust_GFP{13}, Output.SizeClust_YFP{13}, Output.SizeClust_RFP{k}));
    
    if isnan(Output.AvgSizeClust_Total(k))
        Output.AvgSizeClust_Total(k) = 0;
    end

    Output.ImageSizeZ(k) = data.ImageSizeZ(1);
    Output.StepSizeZ(k) = data.StepSizeZ(1);
    Output.VoxelSizeXY(k) = data.VoxelSizeX(1);
    Output.ImageVolumeMM3(k) = (Output.StepSizeZ(k)*1000)*(Output.VoxelSizeXY(k)*1000)*(Output.VoxelSizeXY(k)*1000)*Output.ImageSizeZ(k)*1024*1024;
    Output.Experiment(k) = data.Experiment(1);
    Output.Group(k) = data.Group(1);
    Output.Animal(k) = data.Animal(1);
    Output.Stroke(k) = data.Stroke(1);
    Output.FileName(k) = data.FileName(1);

    Output.Density_Iba1(k) = Output.NumCell_Iba1(k)/Output.ImageVolumeMM3(k);
    Output.Density_Microfetti(k) = Output.NumMicrofetti(k)/Output.ImageVolumeMM3(k);   
    
end

idx = Output.NumClust_CFP == -1;
Output.NumClust_CFP(idx) = 0;

idx = Output.NumClust_GFP == -1;
Output.NumClust_GFP(idx) = 0;

idx = Output.NumClust_YFP == -1;
Output.NumClust_YFP(idx) = 0;

idx = Output.NumClust_RFP == -1;
Output.NumClust_RFP(idx) = 0;

idx = isnan(Output.AvgSizeClust_CFP); 
Output.AvgSizeClust_CFP(idx) = 0;

idx = isnan(Output.AvgSizeClust_GFP); 
Output.AvgSizeClust_GFP(idx) = 0;

idx = isnan(Output.AvgSizeClust_YFP); 
Output.AvgSizeClust_YFP(idx) = 0;  

idx = isnan(Output.AvgSizeClust_RFP); 
Output.AvgSizeClust_RFP(idx) = 0;

Output.NumTotal_Cluster = Output.NumClust_CFP + Output.NumClust_GFP + Output.NumClust_YFP + Output.NumClust_RFP;
Output.Density_Total_Cluster = Output.NumTotal_Cluster./Output.ImageVolumeMM3;
Output.NumSingles_Total = Output.NumSingles_CFP + Output.NumSingles_GFP + Output.NumSingles_YFP + Output.NumSingles_RFP;

disp(head(data))
disp(head(Output))

%% Averaging values for each mouse and save the results in a new table

[C,ia,ic] = unique(Output.Animal);

avg_Output = table('Size',[length(C)*2 29],'VariableTypes',{'categorical', 'categorical', 'categorical', 'logical', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double','double', 'cell', 'double'});
avg_Output.Properties.VariableNames = {'Experiment', 'Group', 'Animal', 'Stroke', 'NumCell_Iba1', 'Density_Iba1', 'NumMicrofetti', 'Density_Microfetti', 'NumCell_CFP', 'NumTotal_Cluster', 'Density_Total_Cluster', 'AvgSizeClust_Total', 'NumSingles_CFP', 'NumClust_CFP', 'AvgSizeClust_CFP', 'NumCell_GFP', 'NumSingles_GFP', 'NumClust_GFP', 'AvgSizeClust_GFP', 'NumCell_YFP', 'NumSingles_YFP', 'NumClust_YFP', 'AvgSizeClust_YFP', 'NumCell_RFP', 'NumSingles_RFP', 'NumClust_RFP', 'AvgSizeClust_RFP', 'SizeClust_Total', 'Sum_NumSingles_Total'};

for i = 1:length(C)


    avg_Output.Experiment(i) = Output.Experiment(ia(i));
    avg_Output.Group(i) = Output.Group(ia(i));
    avg_Output.Animal(i) = Output.Animal(ia(i));
    avg_Output.Stroke(i) = 1;

    avg_Output.NumCell_Iba1(i) = mean(Output.NumCell_Iba1(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.Density_Iba1(i) = mean(Output.Density_Iba1(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.NumMicrofetti(i) = mean(Output.NumMicrofetti(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.Density_Microfetti(i) = mean(Output.Density_Microfetti(Output.Animal == C(i) & Output.Stroke == 1));

    avg_Output.NumTotal_Cluster(i) = mean(Output.NumTotal_Cluster(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.Density_Total_Cluster(i) = mean(Output.Density_Total_Cluster(Output.Animal == C(i) & Output.Stroke == 1));
    
    avg_Output.NumCell_CFP(i) = mean(Output.NumCell_CFP(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.NumSingles_CFP(i) = mean(Output.NumSingles_CFP(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.NumClust_CFP(i) = mean(Output.NumClust_CFP(Output.Animal == C(i) & Output.Stroke == 1));
    x_c = Output.SizeClust_CFP(Output.Animal == C(i) & Output.Stroke == 1);
    y_c = cat(2, x_c{:});
    avg_Output.AvgSizeClust_CFP(i) = mean(y_c);

    avg_Output.NumCell_GFP(i) = mean(Output.NumCell_GFP(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.NumSingles_GFP(i) = mean(Output.NumSingles_GFP(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.NumClust_GFP(i) = mean(Output.NumClust_GFP(Output.Animal == C(i) & Output.Stroke == 1));
    x_g = Output.SizeClust_GFP(Output.Animal == C(i) & Output.Stroke == 1);
    y_g = cat(2, x_g{:});
    avg_Output.AvgSizeClust_GFP(i) = mean(y_g);

    avg_Output.NumCell_YFP(i) = mean(Output.NumCell_YFP(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.NumSingles_YFP(i) = mean(Output.NumSingles_YFP(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.NumClust_YFP(i) = mean(Output.NumClust_YFP(Output.Animal == C(i) & Output.Stroke == 1));
    x_y = Output.SizeClust_YFP(Output.Animal == C(i) & Output.Stroke == 1);
    y_y = cat(2, x_y{:});
    avg_Output.AvgSizeClust_YFP(i) = mean(y_y);


    avg_Output.NumCell_RFP(i) = mean(Output.NumCell_RFP(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.NumSingles_RFP(i) = mean(Output.NumSingles_RFP(Output.Animal == C(i) & Output.Stroke == 1));
    avg_Output.NumClust_RFP(i) = mean(Output.NumClust_RFP(Output.Animal == C(i) & Output.Stroke == 1));
    x_r = Output.SizeClust_RFP(Output.Animal == C(i) & Output.Stroke == 1);
    y_r = cat(2, x_r{:});
    avg_Output.AvgSizeClust_RFP(i) = mean(y_r);


    avg_Output.AvgSizeClust_Total(i) = mean(cat(2, y_c, y_g, y_y, y_r));
    avg_Output.SizeClust_Total{i} = cat(2, y_c, y_g, y_y, y_r);
    avg_Output.Sum_NumSingles_Total(i) = sum(Output.NumSingles_Total(Output.Animal == C(i) & Output.Stroke == 1));

    avg_Output.Experiment(length(C)+i) = Output.Experiment(ia(i));
    avg_Output.Group(length(C)+i) = Output.Group(ia(i));
    avg_Output.Animal(length(C)+i) = Output.Animal(ia(i));
    avg_Output.Stroke(length(C)+i) = 0;

    avg_Output.NumCell_Iba1(length(C)+i) = mean(Output.NumCell_Iba1(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.Density_Iba1(length(C)+i) = mean(Output.Density_Iba1(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.NumMicrofetti(length(C)+i) = mean(Output.NumMicrofetti(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.Density_Microfetti(length(C)+i) = mean(Output.Density_Microfetti(Output.Animal == C(i) & Output.Stroke == 0));

    avg_Output.NumTotal_Cluster(length(C)+i) = mean(Output.NumTotal_Cluster(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.Density_Total_Cluster(length(C)+i) = mean(Output.Density_Total_Cluster(Output.Animal == C(i) & Output.Stroke == 0));
    
    avg_Output.NumCell_CFP(length(C)+i) = mean(Output.NumCell_CFP(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.NumSingles_CFP(length(C)+i) = mean(Output.NumSingles_CFP(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.NumClust_CFP(length(C)+i) = mean(Output.NumClust_CFP(Output.Animal == C(i) & Output.Stroke == 0));
    
    x_c = Output.SizeClust_CFP(Output.Animal == C(i) & Output.Stroke == 0);
    y_c = cat(2, x_c{:});
    avg_Output.AvgSizeClust_CFP(length(C)+i) = mean(y_c);

    avg_Output.NumCell_GFP(length(C)+i) = mean(Output.NumCell_GFP(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.NumSingles_GFP(length(C)+i) = mean(Output.NumSingles_GFP(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.NumClust_GFP(length(C)+i) = mean(Output.NumClust_GFP(Output.Animal == C(i) & Output.Stroke == 0));
   
    x_g = Output.SizeClust_GFP(Output.Animal == C(i) & Output.Stroke == 0);
    y_g = cat(2, x_g{:});
    avg_Output.AvgSizeClust_GFP(length(C)+i) = mean(y_g);

    avg_Output.NumCell_YFP(length(C)+i) = mean(Output.NumCell_YFP(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.NumSingles_YFP(length(C)+i) = mean(Output.NumSingles_YFP(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.NumClust_YFP(length(C)+i) = mean(Output.NumClust_YFP(Output.Animal == C(i) & Output.Stroke == 0));

    x_y = Output.SizeClust_YFP(Output.Animal == C(i) & Output.Stroke == 0);
    y_y = cat(2, x_y{:});
    avg_Output.AvgSizeClust_YFP(length(C)+i) = mean(y_y);

    avg_Output.NumCell_RFP(length(C)+i) = mean(Output.NumCell_RFP(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.NumSingles_RFP(length(C)+i) = mean(Output.NumSingles_RFP(Output.Animal == C(i) & Output.Stroke == 0));
    avg_Output.NumClust_RFP(length(C)+i) = mean(Output.NumClust_RFP(Output.Animal == C(i) & Output.Stroke == 0));

    x_r = Output.SizeClust_RFP(Output.Animal == C(i) & Output.Stroke == 0);
    y_r = cat(2, x_r{:});
    avg_Output.AvgSizeClust_RFP(length(C)+i) = mean(y_r);

    avg_Output.AvgSizeClust_Total(length(C)+i) = mean(cat(2, y_c, y_g, y_y, y_r));
    avg_Output.SizeClust_Total{length(C)+i} = cat(2, y_c, y_g, y_y, y_r);
    avg_Output.Sum_NumSingles_Total(length(C)+i) = sum(Output.NumSingles_Total(Output.Animal == C(i) & Output.Stroke == 0));

    avg_Output.AvgSizeClust_CFP(isnan(avg_Output.AvgSizeClust_CFP))=0;
    avg_Output.AvgSizeClust_GFP(isnan(avg_Output.AvgSizeClust_GFP))=0;
    avg_Output.AvgSizeClust_YFP(isnan(avg_Output.AvgSizeClust_YFP))=0;
    avg_Output.AvgSizeClust_RFP(isnan(avg_Output.AvgSizeClust_RFP))=0;
    avg_Output.AvgSizeClust_Total(isnan(avg_Output.AvgSizeClust_Total ))=0;

end

avg_Output = sortrows(avg_Output,{'Group', 'Animal','Stroke'});

disp(head(avg_Output))

% writetable(avg_Output, 'C:\Users\...\avg_Output.csv'); % insert output directory
% writetable(avg_Output, 'C:\Users\...\avg_Output.xls'); % insert output directory

% This table was used for further statical analysis in GraphPad

%% Plotting histograms of clone size for every animal

% Contralateral hemisphere
figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.5, 0.1])
for i = 1:length(C)
    idx = (avg_Output.Animal == C(i) & avg_Output.Stroke == 0);
    subplot(1,6,i); 
    counts = histcounts(cat(2, ones(1, sum(avg_Output.Sum_NumSingles_Total(idx))), sort(avg_Output.SizeClust_Total{idx})), 'Normalization', 'probability');
    bar(counts, 0.7, 'FaceColor', "#0072BD");
    ylim([0, 1]);
    xlim([0, 15]);
    set(gca,'Xtick',0:5:15)
    title(strcat(string(C(i)), {','}, string(avg_Output.Group(idx))))
    % ylabel('Clone size frequency (%)') 
end

% Stroke hemisphere
figure('Units', 'normalized', 'Position', [0.2, 0.2, 0.5, 0.1])
for i = 1:length(C)
    idx = (avg_Output.Animal == C(i) & avg_Output.Stroke == 1);
    subplot(1,6,i); 
    counts = histcounts(cat(2, ones(1, sum(avg_Output.Sum_NumSingles_Total(idx))), sort(avg_Output.SizeClust_Total{idx})), 'Normalization', 'probability');
    bar(counts, 0.7, 'FaceColor', "#A2142F");
    ylim([0, 1]);
    xlim([0, 15]);
    set(gca,'Xtick',0:5:15)
    title(strcat(string(C(i)), {','}, string(avg_Output.Group(idx))))
    % ylabel('Clone size frequency (%)') 
end

toc