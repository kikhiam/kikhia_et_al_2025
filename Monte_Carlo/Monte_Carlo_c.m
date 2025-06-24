%{
Copyright (c) 2025, Majed Kikhia
All rights reserved.

This source code is licensed under the BSD-style license found in the
LICENSE file in the root directory of this source tree. 

Author: Majed Kikhia
February 2025
%}

%% Calculating densities and plotting

% Input:
%   1. The Matlab table "simulations.mat", which is the ouput of
%      Monte_Carlo_a
%   2. The Matlab table "ring_volumes.mat", which is the ouput of
%      Monte_Carlo_b
% 
% What does the code do?
% It calculates the density of cells with the same color for the recorded
% data and the simulated data and plot the results. 
% 
% Output: Fig. 2c (Here a one animal per time point version)
%
% Running time on demo data = ~ 88 minutes

tic

clc
clear

fprintf('Loading required variables for calculating cell densities.\n');

load("Output\simulations.mat")
load("Output\ring_volumes.mat")

%% Calculating the indices of cells in concentric rings around each cell
% Here we calculate the number of cells per ring for each cell in all the
% images by comparing the pairwise distances to the diameter of the inner 
% and the outer circles of 30 concentric rings places around each cell.   

fprintf('Calculating the indices of cells in concentric rings around each cell.\n');

Num_images = height(simulations); 
[~, Num_simulations] = size(simulations.Simulations{1,1});
indices_of_cells_in_ring = cell(30,Num_images);

r0 = 10:10:300; %Radii of the ring 
r1 = r0 - 10; %Radii of the inner circle of rings 
r2 = r0 + 10; %Radii of the outer circle of rings


for j = 1:Num_images
    for i = 1:30
    indices_of_cells_in_ring{i,j} = simulations.pdist{j}>r1(1,i) & simulations.pdist{j}<=r2(1,i);
    end
end

%% Calculating the density of cells with same color in the real data
% Then we check the color of cells in rings and compare it to the color of
% the cell in question (the cell in the center). This is done by the
% "compare_sum" function, which gives the total number of cells with the 
% same color in an array of real or simulated data.  

fprintf('Calculating the density of cells with same color in the real data.\n');

density_cells_in_ring_with_same_color = zeros(Num_images,30);

for j = 1:Num_images
    x = zeros(simulations.Num_Microfetti(j),30);
    if simulations.Num_Microfetti(j)<=1 % If there is one or zero Confetti cell in an image, there are no cells in ring and density is zero.   
        density_cells_in_ring_with_same_color(j,:) = zeros(1,30);
    else
        for i = 1:30
            x(:,i) = compare_sum(simulations.RealColor{j}, simulations.RealColor{j}, indices_of_cells_in_ring{i, j});
        end
        x = x./ring_volumes{j,1};
        density_cells_in_ring_with_same_color(j,:) = mean(x,1);
    end
end

density_cells_in_ring_with_same_color = density_cells_in_ring_with_same_color * 10^6; % Transfer densities to cell in 0.001 mm3

%% Calculating the density of cells with same color for MC simulations
% This part is similar to the previous section but it caclulate the 
% denities for the simulated data by doing a similar comparison and sum 

% Because of the long running time the results of this section is saved 
% and loaded in the next section. So when the results are available. You
% can skip this section. 


fprintf('Calculating the density of cells with same color in simulated data.\nThis process takes about an hours. Please be patient!\n');

MC_density_simulations = cell(Num_images,1);
y = zeros(Num_simulations,30);

for j = 1:Num_images
    x = zeros(simulations.Num_Microfetti(j),30);
    if simulations.Num_Microfetti(j)<=1 % If there is one or zero Confetti cell in an image, there are no cells in ring and density is zero.   
                MC_density_simulations{j,1} = zeros(1000,30);
    else
        for s = 1:Num_simulations
            for i = 1:30
                x(:,i) = compare_sum(simulations.RealColor{j}, simulations.Simulations{j,1}(:,s), indices_of_cells_in_ring{i, j});
            end
            x = x./ring_volumes{j,1};
            y(s,:) = mean(x,1);
        end
        MC_density_simulations{j,1} = y;
        fprintf('Image %i was analyzed\n', j);
    end
end

save('Output\MC_density_simulations','MC_density_simulations')
clear indices_of_cells_in_ring % To free up memory

%% Averaging by condition
% This part averages the real densities from all images from one mouse 
% by ring. For simulations, it concatenates the data. 

load("Output\MC_density_simulations.mat")

% Averaging images from the same animal
[C,ia,ic] = unique(simulations.Animal);

mean_density_cells_in_ring_with_same_color = zeros(length(C)*2,30);
cat_MC_density_simulations = cell(length(C)*2,1);

for i = 1:length(C)
    mean_density_cells_in_ring_with_same_color(i,:) = mean(density_cells_in_ring_with_same_color((simulations.Animal == C(i) & simulations.Stroke == 'true'), :), 1);
    mean_density_cells_in_ring_with_same_color(length(C)+i,:) = mean(density_cells_in_ring_with_same_color((simulations.Animal == C(i) & simulations.Stroke == 'false'), :), 1);
    cat_MC_density_simulations{i,1} = cat(1, MC_density_simulations{(simulations.Animal == C(i) & simulations.Stroke == 'true'), 1});
    cat_MC_density_simulations{i,1} = cat_MC_density_simulations{i,1} * 10^6; % Transfer densities to cell in 0.001 mm3
    cat_MC_density_simulations{length(C)+i,1} = cat(1, MC_density_simulations{(simulations.Animal == C(i) & simulations.Stroke == 'false'), 1});
    cat_MC_density_simulations{length(C)+i,1} = cat_MC_density_simulations{length(C)+i,1} * 10^6; % Transfer densities to cell in 0.001 mm3
end

%% Adding timepoint and stroke condition to the list of animals for indexing

for i = 1:length(C)
    x = simulations.Group(simulations.Animal == C(i));
    C(i, 2) = x(1);
end 

C(length(C)+1:length(C)*2, :) = C;

C(1:length(C)/2, 3) = 'true';
C((length(C)/2)+1:end, 3) = 'false';

%% Averaging based on timepoint and stroke condition

timepoints = {'2d' '1w' '2w' '4w' '8w' '12w'};

mean_denstiy_group = zeros(length(timepoints)*2, 30);
cat_MC_group = cell(length(timepoints)*2,1);
avg_MC_group = cell(length(timepoints)*2,1);

% Averaging the real data
for i = 1:length(timepoints)
    idx = (C(:,2)==timepoints{i}); 
    mean_denstiy_group(i, :) = mean(mean_density_cells_in_ring_with_same_color((C(:,2)==timepoints{i} & C(:,3) == 'true'),:), 1);
    mean_denstiy_group(length(timepoints) + i, :) = mean(mean_density_cells_in_ring_with_same_color((C(:,2)==timepoints{i} & C(:,3) == 'false'),:), 1);
    cat_MC_group{i, 1} = cat(1, cat_MC_density_simulations{(C(:,2)==timepoints{i} & C(:,3) == 'true'), 1});
    cat_MC_group{length(timepoints) + i, 1} = cat(1, cat_MC_density_simulations{(C(:,2)==timepoints{i} & C(:,3) == 'false'), 1});
end 

% Averaging the simulated data
for i = 1:length(timepoints)
    x = cat(3, cat_MC_density_simulations{(C(:,2)==timepoints{i} & C(:,3) == 'true'), 1});
    avg_MC_group{i, 1} = mean(x,3);
    y = cat(3, cat_MC_density_simulations{(C(:,2)==timepoints{i} & C(:,3) == 'false'), 1});
    avg_MC_group{length(timepoints) + i, 1} = mean(y,3);
end 

clear x y

%% Calculating the median, 2nd and 98th percentiles on averaged data per group

quant_avg_MC_group = cell(length(timepoints)*2,1);

for i = 1:length(quant_avg_MC_group)
    quant_avg_MC_group{i,1} = quantile(avg_MC_group{i,1}, [0.02 0.5 0.98], 1);
end
%% Plotting by group with simulations (Fig. 2c)

fig1 = figure();

tiledlayout(2,3)

warning('off', 'stats:boxplot:BadObjectType');

for i = 1:length(timepoints)   
    nexttile
    plot(20:10:300, mean_denstiy_group(i,2:30), 'Color', '#991b2b', 'LineWidth',2);
    hold on 
    plot(20:10:300, quant_avg_MC_group{i,1}(1,2:30), 'Color', '#f9bda1', 'LineWidth',1);
    plot(20:10:300, quant_avg_MC_group{i,1}(2,2:30), '--', 'Color', '#991b2b', 'LineWidth',1);  
    plot(20:10:300, quant_avg_MC_group{i,1}(3,2:30), 'Color', '#f9bda1', 'LineWidth',1);

    x2 =[20:10:300, fliplr(20:10:300)];
    y2 = [quant_avg_MC_group{i,1}(1,2:30), fliplr(quant_avg_MC_group{i,1}(3,2:30))];
    fill(x2, y2, 'r', 'FaceAlpha',0.1);

    plot(20:10:300, mean_denstiy_group(length(timepoints)+i,2:30), 'Color', '#2d6ba5', 'LineWidth',2);
    plot(20:10:300, quant_avg_MC_group{length(timepoints)+i,1}(1,2:30), 'Color', '#92c7dd', 'LineWidth',1);
    plot(20:10:300, quant_avg_MC_group{length(timepoints)+i,1}(2,2:30), '--', 'Color', '#2d6ba5', 'LineWidth',1);  
    plot(20:10:300, quant_avg_MC_group{length(timepoints)+i,1}(3,2:30), 'Color', '#92c7dd', 'LineWidth',1);

    y3 = [quant_avg_MC_group{length(timepoints)+i,1}(1,2:30), fliplr(quant_avg_MC_group{length(timepoints)+i,1}(3,2:30))];
    fill(x2, y3, 'b', 'FaceAlpha',0.1);


    axis([0 300 0 15])
    title(sprintf('Monte Carlo simulation for time-point %s', timepoints{i}));
    xlabel('Radius of central ring (\mum)') 
    ylabel('Density of Confetti+ microglia in 0.001 mm3') 
    legend({'Stroke','Contralateral'}, 'Location','northeast')
    hold off
end 

toc