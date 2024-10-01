function [slopes] = getSlopes(matrix,vals)
%% getSlopes.m
% Given a matrix in which each column is defined at vals, approximates 
%  slopes at each point using average of left- and right- slopes in a
%  linear interpolation of column

rows = length(vals);
cols = size(matrix,2);

betweenslopes = (matrix(2:rows,:)-matrix(1:rows-1,:))./...
    (repmat(vals(2:rows)-vals(1:rows-1),1,cols));
slopes = nan(rows,cols);
slopes(1,:) = betweenslopes(1,:); 
slopes(rows,:) = betweenslopes(rows-1,:);
slopes(2:(rows-1),:) = 0.5*(betweenslopes(1:rows-2,:)+betweenslopes(2:rows-1,:));