function cellstring = ParseCSL( stringlist,varargin )
%cellstring =  ParseCSL(string_of_comma_separated_list)
%----------------------------------------------------------------------
% Takes a string which contains a list of items separated by commas
% parses each comma separated value into a unique element of a cell array
% of strings.  
% spaces before or after commas are ignored.  
%
% Alistair Boettiger
% boettiger.alistair@gmail.com
% February 12th, 2012

defaults = cell(0,3);
defaults(end+1,:) = {'delimeter', 'string', ','}; % example optional par and its default 

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'a string is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);


items = strfind(stringlist,parameters.delimeter);
if ~isempty(items)
    items = [0,items,length(stringlist)+1];
    Nitems = length(items) -1;
    cellstring = cell(1,Nitems);
    for i=1:Nitems
        cellstring{i} = strtrim(stringlist(items(i)+1:items(i+1)-1));
    end
else
    cellstring = {stringlist};
end
