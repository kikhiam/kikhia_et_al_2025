%{
Copyright (c) 2025, Majed Kikhia
All rights reserved.

This source code is licensed under the BSD-style license found in the
LICENSE file in the root directory of this source tree. 

Author: Majed Kikhia
February 2025
%}

% This function is used by the code Monte_Carlo_c.m
% Input: two arrays of color labels and a logical matrix of indices of 
% cells per ring. 
% Ouput: the number of cells with same color in the ring

function cells_with_same_color = compare_sum(original, sim, indices)

[~,Num_cells] = size(indices);

cells_with_same_color = zeros(Num_cells,1);

for i = 1:Num_cells
    cells_in_ring = sim(indices(:,i));
    cells_with_same_color(i,1) = sum(cells_in_ring == original(i));
end 

