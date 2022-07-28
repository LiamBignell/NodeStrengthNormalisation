function O = TheoreticalStreamlining(T, t, prb)
% This function simulates the seeding of streamlines to create the observed
% matrix O from a theoretical connection matrix T and an optional weighting
% vector w - num streamlines per ROI. 
% See acompying pdf for further explaination.
% The purpose of this is to test ActuallyCorrectMatrixFinder on a
% variety of matrices
%
% T is the input theoretical symmetric connection matrix
% t is an optional vector of roi weights, or number of simulated streamlins
% prb is an optional boolian (default F) 
%
%
% Copyright 2022 Liam Bignell
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

%% Setup
% find matrix size
n = size(T,2);
switch nargin
    case 2
        prb = 0;
    case 1
        prb = 0;
        t = ones(n,1);
end

%% Calculations
% Initialise Output
O = zeros(n);

%Do row-wise: output 
if prb
    for i = 1:n
        % generate the random streamlines by edge wetght
        destn = randsample(1:n, t(i), true, T(i,:));
        % count the streamlines and record
        [GC, GR] = groupcounts(destn');
        O(i,GR) = GC;
    end
else
    for i = 1:n
        % get proportion of node strength for each edge, and adjust by t
        O(i,:) = t(i)*T(i,:)/sum(T(i,:));
    end
end
