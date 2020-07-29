function genes = GetFlyGenes(locustxt,varargin)
% figure(5); clf;
% plot([coding.exonStarts,coding.exonStops],[0,0],'k','lineWidth',5); hold on;
% plot([coding.starts,coding.stops],[0,0],'color',[.5 .5 .5],'LineWidth',2); 
% text(coding.namePos,-.2 + zeros(1,length(coding.symbols)), coding.symbols);

%
defaults = cell(0,3);
defaults(end+1,:) = {'showplot','boolean',false};
defaults(end+1,:) = {'fontSize','nonnegative',10};
defaults(end+1,:) = {'rotation','float',90}; 

pars = ParseVariableArguments(varargin,defaults,mfilename); 

%%
[chr,locusStart,locusEnd] = ParseLocusName(locustxt);

%% version 
global dm3_RNA
if isempty(dm3_RNA)
    dm3_RNA = readtable('C:\Data\Fly\dm3\dm3_NascentRNA.txt');
end
dm3_RNA(1:10,:);

onChr = strcmp(dm3_RNA.chrName,chr);
genes = dm3_RNA(onChr,:);
inRange = genes.startSeq > locusStart & genes.endSeq < locusEnd;
genes = genes(inRange,:);

%%
global dm3_gtf flyNameTable

if isempty(dm3_gtf)
    dm3_gtf = readtable('C:\Data\Fly\dm3\dm3_FlyBaseGenes.txt','ReadVariableNames',false);
    dm3_gtf = dm3_gtf(:,[1:5,7,10,12]);
    dm3_gtf.Properties.VariableNames = {'chr','assembly','featureType','featureStart','featureStop','strand','geneID','transcriptID'};
    % dm3_gtf(randperm(height(dm3_gtf),10),:)
end

if isempty(flyNameTable)
    flyNameTable = readtable('C:\Data\Fly\FlyBase_IDs.txt');
end

onChr = strcmp(dm3_gtf.chr,chr);
exons = dm3_gtf(onChr,:);
inRange = exons.featureStart > locusStart & exons.featureStop < locusEnd;
exons = exons(inRange,:);

isExons = StringFind(exons.featureType,'exon');
exons = exons(isExons,:);


if nargout == 0 || pars.showplot
    plot([exons.featureStart,exons.featureStop],[0,0],'k','lineWidth',5); hold on;
    plot([genes.startSeq,genes.endSeq],[0,0],'k','lineWidth',1); hold on;
    plusGenes = strcmp(genes.strand,'+');
    if ~isempty(genes.startSeq(plusGenes))
        plot(genes.endSeq(plusGenes),zeros(length(genes.endSeq(plusGenes)),1),'k>');
        text(genes.startSeq(plusGenes),-.2+zeros(1,length(genes.Symbol(plusGenes))), genes.Symbol(plusGenes),'HorizontalAlignment','right','Rotation',pars.rotation,'FontSize',pars.fontSize);
    end
    if ~isempty(genes.startSeq(~plusGenes))
        plot(genes.startSeq(~plusGenes),zeros(length(genes.endSeq(~plusGenes)),1),'k<');
        text(genes.endSeq(~plusGenes),-.2+zeros(1,length(genes.Symbol(~plusGenes))), genes.Symbol(~plusGenes),'HorizontalAlignment','right','Rotation',pars.rotation,'FontSize',pars.fontSize);
    end
end

% figure(5); clf;
% plot([coding.exonStarts,coding.exonStops],[0,0],'k','lineWidth',5); hold on;
% plot([coding.starts,coding.stops],[0,0],'color',[.5 .5 .5],'LineWidth',2); 
% text(coding.namePos,-.2 + zeros(1,length(coding.symbols)), coding.symbols);


