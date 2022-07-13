function [T, t] = NSN(O)
%
%This function has been developed for use in structural brain network
%analysis, in particulat, to normalise the streamline count matrix (O)
%produced using probabalistic tractography into the 'true' connectivity
%matrix (T) representative of the true connectivity. 
%
%This is based on 
%modelling probabalistic tractography as each streamline random 
%selecting an edge based on its relative weight for the efferent
%node.
%
%
%O is the (unnormalised) observed connectivity matrix, with O(i,j) being
%the number of streamlines from ROI_i to ROI_j
%
%T is the estimate of the "true" connection matrix
%
%t is a vector of the estimated total node strengths of T
%
%values of t and T are relative.
%
% Copyright 2022 Liam Bignell
%
%This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

%% Get the size of O
n = size(O,2);

%% Normalise by waytotal
for i = 1:n
    O(i,:) = O(i,:)/sum(O(i,:));
end

%% Setup the equations which will be solved to estimate t

% Initialise an empty matrix, which will have one row for each eqn
% equivalently, one row for each pair of upper/lower triangle elements
% equivalently, nchoosek(n,2)
A = zeros(nchoosek(n,2),n); 

% Initialise a counter, so we know where to put each eqn
counter = 0;

%loop over the matrix filling A in At=0
for i = 1:n
    for j = (i+1):n
        counter = counter+1;
        A(counter,i) = O(i,j);
        A(counter,j) = -O(j,i);
    end
end

%% Estimate t
% This is done using the least squared error solution

% Find the matrix in the centre of the t'A'At = MSE
AA = A'*A; 

% Find the eigenvalues and eigenvectors (not in that order) of AA
[V,D] = eig(AA); 

% Find which eigenvalue is smallest
[~, i] = min(abs(diag(D))); 

% The corresponding eigenvector will minimise dMSE/dt
t = V(:,i); 

% I don't believe Matlab will ever make t negative under normal 
% circumstances, as O is strictly non-negative. 
% I have not rigourously tested this.
% all elements should share the same sign, so checking the first 
% element should be sufficient
% I have also not rigourously tested this.
if t(1) < 0 
    t = -t;
end

% re-normalise t so the mean value is 1 (rather than 1/n)
t = t.*n;  

%% Find T
T = zeros(n);
for i = 1:n
    T(i,:) = O(i,:)*t(i);
end
