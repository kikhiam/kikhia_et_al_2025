%{
Copyright (c) 2025, Majed Kikhia
All rights reserved.

This source code is licensed under the BSD-style license found in the
LICENSE file in the root directory of this source tree. 

Author: Majed Kikhia
February 2025
%}

% This function is used by the code Monte_Carlo_b.m
% Input: the location of cells, the number of cells (n), image
% wdith (x), image height (y), image depth (z), and the diameters of the 
% inner edge (r1) and the out edge (r2) of the concentric rings.
% Ouput: a matrix with n x 30 with ring volume of all cells (n) in the image.

function ring_volumes = rings_antenna(cells, N, image_w, image_h, image_d, r1, r2)
xq = [0, 0, 0; image_w, 0, 0; image_w, image_h, 0; 0, image_h, 0];
poly1 = antenna.Polygon('Vertices',xq);
ring_area = zeros(N,30);
for j = 1:N
    for i = 1:30
        if i == 1
            circle2 = antenna.Circle('Center',[cells(j,1) cells(j,2)],'Radius',r2(1,i));
            intersect = and(poly1,circle2);
            ring_area(j,i) = area(intersect);
        else
            circle1 = antenna.Circle('Center',[cells(j,1) cells(j,2)],'Radius',r1(1,i));
            circle2 = antenna.Circle('Center',[cells(j,1) cells(j,2)],'Radius',r2(1,i));
            ring = minus(circle2,circle1);
            intersect = and(poly1,ring);
            ring_area(j,i) = area(intersect);
        end
    end
end
ring_volumes = ring_area*image_d;
end
