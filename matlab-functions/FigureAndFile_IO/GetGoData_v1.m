function [GO,humanGO,mapName2ID,mapID2Name,mapGenes2GOids,mapGOids2Genes, mapID2Ontology] = GetGoData_v1(libGenes, aspect)

if nargin < 2
    aspect = 'FPC';
end

% Download from http://geneontology.org/ontology/go-basic.obo
% GO = geneont('live',true); % this step takes a while
GO =geneont('File','\\MORGAN\TSTORMdata2\GenomeData\HumanGenome\go-basic.obo');

% Download from http://geneontology.org/gene-associations/gene_association.goa_human.gz
humanGOfile = '\\MORGAN\TSTORMdata2\GenomeData\HumanGenome\gene_association.goa_human';
humanGO = goannotread(humanGOfile,'Aspect',aspect,'Fields',{'DB_Object_Symbol','GOid'});

% F = molecular function, P = process, C = cellular compartment

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
for i = 1:numel(humanGO)
    key =  humanGO(i).DB_Object_Symbol;  %
    if isKey(mapGenes2GOids,key)
        mapGenes2GOids(key) = [mapGenes2GOids(key) humanGO(i).GOid];
    else
        mapGenes2GOids(key) = humanGO(i).GOid;
    end
end

fprintf('Number of annotated genes related to molecular function is %d.\n',mapGenes2GOids.Count)
fprintf('Number of unique GO terms associated to annotated genes is %d.\n',numel(unique([humanGO.GOid])))
fprintf('Number of gene-GO term associations is %d.\n',numel(humanGO))

% this is lib genes specific
numGenes = length(libGenes);
mapGOids2Genes = containers.Map();
for i = 1:numGenes % i=43
    try
        goid = getrelatives(GO,mapGenes2GOids(libGenes{i}));
        keys =  arrayfun(@(x) mapID2Name(x),goid,'UniformOutput',false,'ErrorHandler',@ErrorFunEmptyString);
        keys = unique(keys); 
        for k=1:length(keys)
            if isKey(mapGOids2Genes,keys{k})
                currTerms = mapGOids2Genes(keys{k});
                mapGOids2Genes(keys{k}) = {currTerms{:}, libGenes{i}};
            else
                mapGOids2Genes(keys{k}) = libGenes(i);
            end
        end
    catch er
        disp(['Gene ',libGenes{i},' num ',num2str(i),' ',er.message])
    end
end
