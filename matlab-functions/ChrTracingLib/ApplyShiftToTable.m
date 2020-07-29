function spotTable = ApplyShiftToTable(spotTable,otChn,tform3D,varargin)
% adds / updates columns xcShift, ycShift, zcShift by transforming the data
% in spotTable.x,y,z and spotTable.locusX,Y with the 'tform3D' for the
% specified channel 'otChn'
% 
% called by:
%  ChromaticAlignFromTable
% 
% Alistair Boettiger
% CC BY Jan 2019


isOt = spotTable.chn == otChn;
stageXYZ = [spotTable.x(isOt) + spotTable.locusX(isOt),...
            spotTable.y(isOt) + spotTable.locusY(isOt),...
            spotTable.z(isOt)];
newStageXYZ = tforminv(tform3D,stageXYZ(:,1),stageXYZ(:,2),stageXYZ(:,3));
spotTable.xcShift(isOt) = newStageXYZ(:,1) - stageXYZ(:,1);
spotTable.ycShift(isOt) = newStageXYZ(:,2) - stageXYZ(:,2);
spotTable.zcShift(isOt) = newStageXYZ(:,3) - stageXYZ(:,3);