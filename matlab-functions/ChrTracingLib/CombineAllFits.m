function [polys,maps,spotData] = CombineAllFits(dataFolder,varargin)
% load all data files from folder

defaults = cell(0,3);
defaults(end+1,:) = {'showPlots','integer',0};
defaults(end+1,:) = {'byFOV','boolean',false};
% TableToPolymer parameters
defaults(end+1,:) = {'dims',{'xyz','xy','xz'},'xyz'};
defaults(end+1,:) = {'bins','nonnegative',0}; % 0 for auto detect
defaults(end+1,:) = {'chromCorrect','boolean',true};
defaults(end+1,:) = {'computeDistMap','boolean',true};
defaults(end+1,:) = {'sortMethod',{'byHybe','byRead'},'byRead'};
defaults(end+1,:) = {'selectChn','nonnegative',0}; % 0 for all chn
defaults(end+1,:) = {'dataBasedDriftCorrect','boolean',false};
defaults(end+1,:) = {'removeBlank','boolean',false};
defaults(end+1,:) = {'shiftSign','integer',1};
defaults(end+1,:) = {'parallel','integer',1};
pars = ParseVariableArguments(varargin,defaults,mfilename);

allFitFiles = FindFiles([dataFolder,'*AllFits.csv']);
nFits = length(allFitFiles);
polys = cell(nFits,1);
maps = cell(nFits,1);
spotData = cell(nFits,1);
if pars.parallel == 1
    for f=1:nFits % f=1
        [polys{f},maps{f},spotXY] = TableToPolymer(allFitFiles{f},'parameters',pars);
        spotData{f} = [f*ones(size(spotXY,1),1),spotXY];
        if pars.showPlots > 0
            figure(pars.showPlots); clf; 
            imagesc(nanmedian(maps{f},3)); colorbar; caxis([0,800]); 
            title(['f=',num2str(f),' n=',num2str(size(maps{f},3))]);
            pause(.1);    
        end
    end
else
    parfor f=1:nFits % f=1
        [polys{f},maps{f}] = TableToPolymer(allFitFiles{f},'parameters',pars);
    end
end

if ~pars.byFOV
    try
        polys = cat(3,polys{:});
        maps = cat(3,maps{:});
        spotData = cat(1,spotData{:});
    catch
%         matSize = cellfun(@size,spotData,'UniformOutput',false);
%         matSize = cat(1,matSize{:})
        warning('unable to concatinate FOVs, different numbers of hybs found');
        cprintf([1 .5 0],'try setting "bins" equal to number of readouts+replicates expected');
    end
end