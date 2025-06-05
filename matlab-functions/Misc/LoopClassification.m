function [resStruct] = LoopClassification(data, nB)
% LoopClassification classifies loop shapes in the data
% This function pass data variable into LoopClassification.py 
% and restructure the Pythonic outputs into a MATLAB struct object.
% 
% This function use parfor to speed up for-loop using parallelization. 
% It has an upfront amortized cost to find MATLAB workers, but overall 
% this significantly speed up the for loop. 
% 
% Note: This function does not reveal nested loop children.
%
% Parameters
% ----------
% data - an (nx1) cell array [n = the number of single cells in the sample]
% each cell entry contains an (mx2) matrix where each row is the coordinate
% of a thresholded contact.
%
% nB - an integer: the number of readout probes 
%
% Returns
% -------
% resStruct - a (nx1) struct containing 5 fields [n = the number of single
% cell]. Fieldnames are nestedLoopBase, crossLoop, isolatedLoop,
% multiContact, and nestedLoopChildren

% note NestedLoopChildren is buggy

nsc = length(data);

nestedLoopBaseCell = cell(nsc, 1);
crossLoopCell = cell(nsc, 1);
isolatedLoopCell = cell(nsc, 1);
multiContactCell = cell(nsc, 1);
nestedLoopChildrenCell = cell(nsc, 1); 

parfor isc = 1:nsc
    dataSc = data{isc};
    if height(dataSc) > 1
        res = struct(pyrunfile("LoopClassification.py", "loopType", data=dataSc, Mb=nB));
        
        nestedLoopBasePy = py.list(keys(res.nestedLoopDict));
        crossLoopPy = py.list(keys(res.crossLoopDict));
        isolatedLoopPy = py.list(res.isolatedLoopSet);
        multiContactPy = py.list(res.multiContactSet);
        nestedLoopChildrenPy = py.list(values(res.nestedLoopDict));

        nLoopChildren = 0;
        for i = 1:length(nestedLoopChildrenPy)
            currBase = nestedLoopChildrenPy{i};
            for j = 1:length(currBase)
                nLoopChildren = nLoopChildren + 1;
            end
        end

        nestedLoopBaseCoord = zeros(int8(py.len(nestedLoopBasePy)), 2);
        crossLoopCoord = zeros(int8(py.len(crossLoopPy)), 2);
        isolatedLoopCoord = zeros(int8(py.len(isolatedLoopPy)), 2);
        multiContactCoord = zeros(int8(py.len(multiContactPy)), 2);
        nestedLoopChildrenCoord = zeros(nLoopChildren, 2);
    
        for i = 1:length(nestedLoopBasePy)
            nestedLoopBaseCoord(i, 1) = int64(nestedLoopBasePy{i}{1})+1;
            nestedLoopBaseCoord(i, 2) = int64(nestedLoopBasePy{i}{2})+1;
        end 
    
        for i = 1:length(crossLoopPy)
            crossLoopCoord(i, 1) = int64(crossLoopPy{i}{1})+1;
            crossLoopCoord(i, 2) = int64(crossLoopPy{i}{2})+1;
        end
    
        for i = 1:length(isolatedLoopPy)
            isolatedLoopCoord(i, 1) = int64(isolatedLoopPy{i}{1})+1;
            isolatedLoopCoord(i, 2) = int64(isolatedLoopPy{i}{2})+1;
        end
    
        for i = 1:length(multiContactPy)
            multiContactCoord(i, 1) = int64(multiContactPy{i}{1})+1;
            multiContactCoord(i, 2) = int64(multiContactPy{i}{2})+1;
        end

        for i = 1:length(nestedLoopChildrenPy)
            currBase = nestedLoopChildrenPy{i};
            for j = 1:length(currBase)
                matrix = [int64(currBase{j}{1})+1 int64(currBase{j}{2})+1];
                nestedLoopChildrenCoord(j, :) = matrix;
            end
        end 
    
        nestedLoopBaseCell{isc, 1} = nestedLoopBaseCoord;
        crossLoopCell{isc, 1} = crossLoopCoord;
        isolatedLoopCell{isc, 1} = isolatedLoopCoord;
        multiContactCell{isc, 1} = multiContactCoord; 
        nestedLoopChildrenCell{isc, 1} = nestedLoopChildrenCoord;
    end
end 

resStruct = struct("nestedLoopBase", nestedLoopBaseCell, "crossLoop", crossLoopCell, ...
    "isolatedLoop", isolatedLoopCell, "multiContact", multiContactCell, ...
    "nestedLoopChildren", nestedLoopChildrenCell);

end 