function [Tdebt] = genSparseTransition(znext,z)
%% genSparseTransition.m
% Generates sparse transition matrix given z and policy function znext
% Note that if znext(z^n) = nan, row of Tdebt corresponding to z^n is set
%  to 0.

M = length(z);

ignore = isnan(znext);
temp = 1:M;
keepInd = temp(~ignore);

% assign znext_e to bins defined by z 
[~,~,znext_bin] = histcounts(znext(~ignore),z);

% turn these bin assignments into transition probabilities, appropriately
%  weighted for each edge of the bin as recommended by Kopecky
Tdebt = sparse([keepInd,keepInd],[znext_bin,[znext_bin+1]],...
    [(z(znext_bin+1) - znext(~ignore))./...
    (z(znext_bin+1)-z(znext_bin)),...
    (znext(~ignore) - z(znext_bin))./...
    (z(znext_bin+1)-z(znext_bin))],M,M);