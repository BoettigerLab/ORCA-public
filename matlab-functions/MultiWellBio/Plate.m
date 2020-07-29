function [ind, parameters] = Plate(rows, columns, varargin)
%--------------------------------------------------------------------------
% ind = Plate(rows, columns)
% This function converts a vector of columns or a vector of rows into
% indices for a 96 well plate
%--------------------------------------------------------------------------
% Outputs:
% ind: The indices for a 96 well plate
%
%--------------------------------------------------------------------------
% Inputs:
% columns: A list of column values
%
% row: A list of row values (A=1, B=2, ...)
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% September 12, 2012
%
% Version 1.0
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'plateDim', 'array', [12 8]};
defaults(end+1,:) = {'indexOrder', {'rows', 'columns'}, 'rows'};
defaults(end+1,:) = {'invertPlate', 'boolean', false}; % Start with A1 or G12

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%--------------------------------------------------------------------------
% Check Length
%--------------------------------------------------------------------------
ind = zeros(length(columns), length(rows));

for i=1:length(columns)
    for j=1:length(rows)
        if ~parameters.invertPlate
            ind(i,j) = sub2ind(parameters.plateDim, columns(i), rows(j));
        else
            ind(i,j) = sub2ind(parameters.plateDim, parameters.plateDim(1)-columns(i) + 1, ...
                parameters.plateDim(2) - rows(j) + 1);
        end
    end
end

%--------------------------------------------------------------------------
% Flatten Indices
%--------------------------------------------------------------------------
switch parameters.indexOrder
    case 'rows'
        ind = reshape(ind, [1 prod(size(ind))]);
    case 'columns'
        ind = reshape(ind', [1 prod(size(ind))]);
end
