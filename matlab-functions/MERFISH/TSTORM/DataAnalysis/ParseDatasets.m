function datasets2 = ParseDatasets(datasets,varargin)
% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'analysisType', {'Analysis','RP_Analysis'}, 'Analysis'};
defaults(end+1,:) = {'numWords', 'positive', []};
defaults(end+1,:) = {'keepSpecialSets', 'boolean', false};


% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'requires flists, mRNAcents');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments([], defaults, mfilename);


%%
% specify data sets of interest
chosenSets = true(1,length(datasets));
if ~isempty(parameters.numWords)
    chosenSets = chosenSets & [datasets.numWords] == parameters.numWords;
end
% if ~parameters.keepSpecialSets
%     chosenSets = chosenSets & strcmp({datasets.SpecialTx},'none');
% end

% specify pipeline to import data from
analysisType = parameters.analysisType; % 'RP_Analysis'; % 'Analysis';

% specify variables to import (only for analysisType = 'Analysis')
% loading all the data types will be slow and I rarely use them all. 
dataToLoad = {'rnaCnts','perCnts','libGenes','libExpect'}; 


% grab just the indicated sets
datasets2 = datasets(chosenSets);

% loop through datasets find most recent analysis folder
% load indicated variables
for i=1:length(datasets2) % i=5
    analysisFolders = dir([datasets2(i).dataPath,filesep,analysisType,'*']);
    if isempty(analysisFolders)
        warning(['no Analysis folders found in ',datasets2(i).dataPath]);
    else
        [~,sortIdx] = sort([analysisFolders(:).datenum]);
        % [~,mostRecent] = max([analysisFolders(:).datenum]);
        sortedAnalysis = strcat([datasets2(i).dataPath,filesep],{analysisFolders(sortIdx).name},[filesep]);
        
        if strcmp(analysisType,'Analysis')
            kk=0;
            while kk<length(sortedAnalysis)
                try 
                    sortedAnalysis = [datasets2(i).dataPath,filesep,analysisFolders(end-kk).name,filesep];
                    load([sortedAnalysis,'cntData.mat'],dataToLoad{:});
                    disp(['Loaded ', sortedAnalysis]); 
                    datasets2(i).analysisPath = analysisFolders(end-kk).name;
                    kk = inf;
                catch
                    kk = kk+1;
                    % disp(['No cntData in ', sortedAnalysis]); 
                end
            end
            for d=1:length(dataToLoad)
                datasets2(i).(dataToLoad{d}) = eval(dataToLoad{d});
            end    
        elseif strcmp(analysisType,'RP_Analysis')
            fpkmReport = LoadByteStream([sortedAnalysis,'FPKMReport.matb']); 
            fields = fieldnames(fpkmReport);
            for d=1:length(fields)
                datasets2(i).(fields{d}) = fpkmReport.(fields{d});
            end    
        end
    end
end
    
    