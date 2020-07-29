function sumMap = AddMaps(map1, map2)
% ------------------------------------------------------------------------
% sumMap = AddMaps(map1, map2)
% This function adds to 
%--------------------------------------------------------------------------
% Necessary Inputs
% map1, map2: Container map objects
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% April 18, 2015
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Check for empty maps
%--------------------------------------------------------------------------
if isempty(map2)
    sumMap = map1;
    return;
end
if isempty(map1)
    sumMap = map2;
    return;
end

%--------------------------------------------------------------------------
% Find overlap and distinct subsets of two maps
%--------------------------------------------------------------------------
% cell2mat is slow!
keys1 = cell2mat(keys(map1));
keys2 = cell2mat(keys(map2));
values1 = cell2mat(values(map1));
values2 = cell2mat(values(map2));

[~, nsInds1, nsInds2] = setxor(keys1, keys2);
[sharedKeys, sInds1, sInds2] = intersect(keys1, keys2);

sharedValues = values1(sInds1) + values2(sInds2);
sumMap = containers.Map([keys1(nsInds1) keys2(nsInds2) sharedKeys], ...
    [values1(nsInds1) values2(nsInds2) sharedValues]);


