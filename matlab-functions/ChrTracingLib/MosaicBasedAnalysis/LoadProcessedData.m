function embDataAbd = LoadProcessedData(loadFolder,varargin)
% embDataAbd = LoadProcessedData(loadFolder

defaults = cell(0,3);
defaults(end+1,:) = {'tag','string',''};
defaults(end+1,:) = {'embs','integer',[]};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if isempty(pars.embs)
    numEmbs = length(cellstr(ls([loadFolder,pars.tag,'*_DataStruct.mat'])));
    embs = 1:numEmbs;
else
    embs = pars.embs;
end

embDataAbd = [];
for e = embs
    try
        load([loadFolder,pars.tag,'Emb',num2str(e,'%02d'),'_DataStruct.mat'],'embDat');
        
        % add missing fields to incomplete files to allow concatination
        if ~isfield(embDat,'embSegmentIDs')
            embDat.embSegmentIDs = [];
        end
        if ~isfield(embDat,'embPolygonData')
            embDat.embPolygonData = [];
        end
        
        if isempty(embDataAbd)
            embDataAbd = embDat;
        else
            embDataAbd(end+1) = embDat; %#ok<AGROW>
        end
    catch er
        disp(er.message);
    end     
end
