function [reports,goMaps] = GOcluster(geneList,clustCond,varargin)
% takes a list of genes, clustered into N groups indexed by clustCond, and
% finds GO terms statistically enriched in each group relative to the whole
% 
% inputs list of genes
% list of clusters for the genes
%
% Optional Inputs
% 
% 

% some global here prevent having to reload stuff 
global GO goGenome

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'sortMethod',{'pvalue','count'},'pvalue'};
defaults(end+1,:) = {'goMaps','freeType',[]}; % 
defaults(end+1,:) = {'aspect','string','FPC'};
defaults(end+1,:) = {'goPath','string','U:\GenomeData\GeneOntology\go-basic.obo'};
defaults(end+1,:) = {'goGenome',{'mm','hg','dm'},'mm'};
defaults(end+1,:) = {'mmGO','string','U:\GenomeData\GeneOntology\mgi.gaf'};
defaults(end+1,:) = {'dmGO','string','U:\GenomeData\GeneOntology\fb.gaf'};
defaults(end+1,:) = {'hgGO','string','U:\GenomeData\GeneOntology\go_human.gaf'};
pars = ParseVariableArguments(varargin,defaults,mfilename);
% geneList = '';

goMaps = pars.goMaps;
if isempty(goMaps)

    % Download from http://geneontology.org/ontology/go-basic.obo
    % GO = geneont('live',true); % this step takes a while
    if isempty(GO)
        try
            tic
            GO =geneont('File',pars.goPath);
            disp('finished loading GO basic');
            toc;       
        catch
           error('error loading GO database. Update "goPath" flag to path to go-basic.obo.  You can download it here: Download from http://geneontology.org/ontology/go-basic.obo');
        end
    end

    % Download from http://geneontology.org/gene-associations/gene_association.goa_human.gz
    if strcmp(pars.goGenome,'mm')
        goFile = pars.mmGO;   
    elseif strcmp(pars.goGenome,'hg')
        goFile = pars.hgGO;
    elseif strcmp(pars.goGenome,'dm')
        goFile = pars.dmGO;
    end

    if isempty(goGenome)
        try
            tic
            goGenome = goannotread(goFile,'Aspect',pars.aspect,'Fields',{'DB_Object_Symbol','GOid'}); 
            disp('finished loading GO genome annoations');
            toc
            % F = molecular function, P = process, C = cellular compartment
        catch er
            warning(er.message);
            error('error loading GO genome. Update "goGenome" flag to "mm","hg" or "dm".  You can download the GO GAF files for your genome here: http://geneontology.org/');
        end
    end



    % look-up table of all go term IDs to names and vice versa
    terms = GO.terms;
    allNames = {};
    allOntologies = {};
    allIDs = [];
    for i=1:length(terms)
        allNames{i} = terms(i).name;
        allOntologies{i} = terms(i).ontology;
        allIDs(i) = terms(i).id;
    end
    mapName2ID = containers.Map(allNames, allIDs);
    mapID2Name = containers.Map(allIDs, allNames);
    mapID2Ontology = containers.Map(allIDs, allOntologies);
    % look up table of the GO IDs for all human genes
    mapGenes2GOids = containers.Map();
    for i = 1:numel(goGenome)
        key =  goGenome(i).DB_Object_Symbol;  %
        if isKey(mapGenes2GOids,key)
            mapGenes2GOids(key) = [mapGenes2GOids(key) goGenome(i).GOid];
        else
            mapGenes2GOids(key) = goGenome(i).GOid;
        end
    end

    if pars.verbose
        fprintf('Number of annotated genes related to molecular function is %d.\n',mapGenes2GOids.Count)
        fprintf('Number of unique GO terms associated to annotated genes is %d.\n',numel(unique([goGenome.GOid])))
        fprintf('Number of gene-GO term associations is %d.\n',numel(goGenome))
    end
    %% create look up maps
    % also note which terms are not present in the list

    numGenes = length(geneList);
    mapGOids2Genes = containers.Map();
    keyedGenes = cell(numGenes,1); % create a new list which only has the genes that are in the database.
    missingGenes = cell(numGenes,1); % create a new list which only has the genes that are in the database.
    missingID = false(numGenes,1);
    for i = 1:numGenes % i=43
        try
            goid = getrelatives(GO,mapGenes2GOids(geneList{i}));
            keys =  arrayfun(@(x) mapID2Name(x),goid,'UniformOutput',false,'ErrorHandler',@ErrorFunEmptyString);
            keys = unique(keys); 
            for k=1:length(keys)
                if isKey(mapGOids2Genes,keys{k})
                    currTerms = mapGOids2Genes(keys{k});
                    mapGOids2Genes(keys{k}) = {currTerms{:}, geneList{i}};
                else
                    mapGOids2Genes(keys{k}) = geneList(i);
                end
            end
            keyedGenes{i} = geneList{i};
        catch er
            missingGenes{i}= geneList{i};
            missingID(i) = true;
            disp(['Gene ',geneList{i},' num ',num2str(i),' ',er.message])
        end
    end

    keyedGenes( cellfun(@isempty,keyedGenes)) = [];
    missingGenes( cellfun(@isempty,missingGenes)) = [];
    familyIdx = clustCond;
    familyIdx(missingID) = [];
    libRealGenes = keyedGenes;
else
    libRealGenes = geneList; % maybe error
    familyIdx = clustCond;
    mapName2ID = goMaps.mapName2ID;
    mapID2Name = goMaps.mapID2Name;
    mapID2Ontology = goMaps.mapID2Ontology;
    mapGenes2GOids = goMaps.mapGenes2GOids; 
    mapGOids2Genes = goMaps.mapGOids2Genes;
end

%% get N most signficant GO terms for each group
numTerms =30;

numGroups = max(familyIdx);
numRealGenes = length(libRealGenes);
reports = cell(1,numGroups);

m = GO.Terms(end).id;           % gets the last term id
geneschipcount = zeros(m,1);    % a vector recording for each GO-term the number of times it was used in the combined geneList
genesclustercount = zeros(m,numGroups); %  a vector recording for each GO-term the number of times it was used in the cluster
for i = 1:numel(libRealGenes)
    if isKey(mapGenes2GOids,libRealGenes{i})
        goid = mapGenes2GOids(libRealGenes{i}); % getrelatives(GO,mapGenes2GOids(libRealGenes{i}));
        geneschipcount(goid) = geneschipcount(goid) + 1;
        genesclustercount(goid,familyIdx(i)) = genesclustercount(goid,familyIdx(i)) + 1;
    end
end
%% show most common GO terms among the entire gene list
% clc;
% most = find(geneschipcount > 40)
% vs = cell(length(most),1);
% for i=1:length(most)
%     vs{i} = mapID2Name(most(i));
% end
% disp(vs)

%%

pvalues = ones(m,numGroups);
grp = cell(numGroups,1); 
for g=1:numGroups
    grp{g} = libRealGenes(familyIdx==g);
    pvalues(:,g) = hygepdf(genesclustercount(:,g),numRealGenes,...
                           length(grp{g}),geneschipcount);
end

clc;
topPvalues = zeros(numTerms,numGroups);
topEvalues = zeros(numTerms,numGroups);
topGoIds = zeros(numTerms,numGroups);
clusterCnt = zeros(numTerms,numGroups);
globalCnt = zeros(numTerms,numGroups);
topGoNames = cell(numGroups,1);
for g=1:numGroups % g = 8; 
    enrichment = genesclustercount(:,g)/length(grp{g})./(geneschipcount/numRealGenes); 
    enrichment(isnan(enrichment)) = 0; 
    if strcmp(pars.sortMethod,'count')
        relAbund =  genesclustercount(:,g);
        relAbund(relAbund < 5) = 0;
        absAbund = geneschipcount;
        abundScore = relAbund./absAbund;
        abundScore( pvalues(:,g) < .05 | isnan(abundScore)) = 0;
        [as,sortIdxP] = sort(abundScore,'descend'); 
    elseif strcmp(pars.sortMethod,'pvalue')
        idx = enrichment > 1  & pvalues(:,g) < .05; % 
        genePvalues =  pvalues(:,g);
        genePvalues(~idx) = 1;
        maxTerms = min(numTerms,sum(idx));
        [spvalues,sortIdxP] = sort(genePvalues);
    end
    topPvalues(:,g) = spvalues(1:numTerms);
    topGoIds(:,g) = sortIdxP(1:numTerms);
    clusterCnt(1:maxTerms,g) = genesclustercount( sortIdxP(1:maxTerms),g );
    globalCnt(1:maxTerms,g) = geneschipcount(  sortIdxP(1:maxTerms) );
    topEvalues(:,g) = (clusterCnt(:,g)/length(grp{g})) ./ (globalCnt(:,g)/numRealGenes);
    topGoNames{g} = arrayfun(@(x) mapID2Name(x),sortIdxP(1:maxTerms),'UniformOutput',false,'ErrorHandler',@ErrorFunEmptyString);
    topGoNames{g};
    
    disp(['Group ',num2str(g)]);
    disp(sort(grp{g}));
    disp(['Group ',num2str(g)]);
    report = sprintf('GO Term      p-val e-val  counts  definition\n');
    for j = 1:maxTerms % '%s%s\t%-1.4f\t%-d / %-d\t%s\n'
        report = sprintf('%s%s\t%-1.4f\t%-1.4f\t%-d / %-d\t%s\n', report, ...
                    char(num2goid( topGoIds(j,g))),...
                    topPvalues(j,g),...    
                    topEvalues(j,g),...
                    clusterCnt(j,g),....
                    globalCnt(j,g),...
                    topGoNames{g}{j});
    end
    disp(report); 
    reports{g} = report;
end

goMaps.mapName2ID = mapName2ID;
goMaps.mapID2Name = mapID2Name;
goMaps.mapID2Ontology = mapID2Ontology;
goMaps.mapGenes2GOids = mapGenes2GOids;
goMaps.mapGOids2Genes = mapGOids2Genes;
