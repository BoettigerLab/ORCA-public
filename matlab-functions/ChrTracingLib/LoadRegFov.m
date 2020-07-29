function [xyShifts, xyShiftsOut, hasData] = LoadRegFov(regFovCsv)
% xyShifts - a Nx2 table of shifts for aligning FOVs within a hyb.
% xyShiftsOut - a Nx2 table of shifts with NaNs marking unprocessed data
% regFovCsv - an Nx3 table, recording xShift, yShift, tileOrder
%    In this table unprocessed data is entered as NaN for the shift.

tableData = readtable(regFovCsv);
% reading some tables with NaN's from fprintf leads to cell/double issues
if iscell(tableData{2,1})
   xShifts = cellfun(@str2double,tableData.Var1(2:end));
   yShifts = cellfun(@str2double,tableData.Var2(2:end));
   xyShifts = [xShifts,yShifts];
else
   xShifts = tableData.xShift;
   yShifts = tableData.yShift;
   xyShifts = [xShifts,yShifts];
end   
xyShiftsOut = xyShifts;
nTiles = size(xyShifts,1);
hasData = true(nTiles,1);
noDataRows = any(isnan(xyShifts),2); 
hasData(noDataRows) = false; % default state was false
xyShifts(noDataRows,:) = 0;
