function [fof_core,fof_spot_quality] = AllFitsToFOF(folder,locusName,stepSize,varargin)
%
% 
% Updates
%  1/18/2023 - added fov to spot_quality table


defaults = cell(0,3);
defaults(end+1,:) = {'fovs', 'nonnegative', 0}; 
defaults(end+1,:) = {'rescale', 'nonnegative', 1}; 
defaults(end+1,:) = {'verbose', 'boolean', true}; 
pars  = ParseVariableArguments(varargin, defaults, mfilename);

fitsFound = FindFiles([folder,'fov*_AllFits.csv']);
if pars.fovs == 0
    fovs = 1:length(fitsFound);
elseif isinf(pars.fovs(end))
    fovs = pars.fovs(1):length(fitsFound);
else
    fovs = pars.fovs;
end
nFOV = length(fovs);
[chrName,chrStart] = ParseLocusName(locusName); 
fof_core = cell(nFOV,1);
fof_spot_quality = cell(nFOV,1);
currTrace = 0;
for f=fovs
    fit_data = [folder,'fov',num2str(f,'%03d'),'_AllFits.csv'];
    fit_table = readtable(fit_data);
    x = fit_table.x*pars.rescale + fit_table.locusX;
    y = fit_table.y*pars.rescale + fit_table.locusY;
    z = fit_table.z*pars.rescale;
    readout = fit_table.readout;
    trace_ID = fit_table.s + currTrace;
    brightness = fit_table.h;
    fit_height = fit_table.a;
    fit_background = fit_table.b;
    fit_width_x = fit_table.wx;
    fit_width_y = fit_table.wy;
    fit_width_z = fit_table.wz;
    drift_x = fit_table.xshift;
    drift_y = fit_table.yshift;
    drift_z = fit_table.zshift;
    fov = fit_table.fov;
    nSpots = length(x);
    chr = repmat(chrName,nSpots,1);
    chr_start = (readout-1)*stepSize + chrStart;
    chr_end = readout*stepSize + chrStart;
    fof_core{f} = table(x,y,z,chr,chr_start,chr_end,trace_ID);
    fof_spot_quality{f} = table(brightness,fit_height,fit_background,fit_width_x,fit_width_y,fit_width_z,drift_x,drift_y,drift_z,fov);
    currTrace = currTrace + length(unique(trace_ID));
end
fof_core = cat(1,fof_core{:});
fof_spot_quality = cat(1,fof_spot_quality{:});
totSpots = height(fof_core);
spot_ID = (1:totSpots)';
spot_ID = table(spot_ID);
fof_core = cat(2,spot_ID,fof_core);
fof_spot_quality = cat(2,spot_ID,fof_spot_quality);

if pars.verbose
    fof_core([1:2,totSpots-1:totSpots],:)
    fof_spot_quality([1:2,totSpots-1:totSpots],:)
end