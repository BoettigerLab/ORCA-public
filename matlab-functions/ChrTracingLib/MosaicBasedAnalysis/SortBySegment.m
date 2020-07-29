function [segMap,segName] = SortBySegment(eMapsAbd,eTablesAbd,varargin)
% separate into Head, PS7 - PS14
% plot if requested

defaults = cell(0,3);
defaults(end+1,:) = {'badReads','integer',[]};
defaults(end+1,:) = {'rs','integer',[]};
defaults(end+1,:) = {'showPlot','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename); 


if isempty(pars.rs)
    rs = 1:size(eMapsAbd,1);
else
    rs = pars.rs;
end
rs(pars.badReads) = [];  % [16,26,45,21,29,31,39];


k = eTablesAbd.emb < 20;
segMap = cell(12,1);
segName = cell(12,1);
for i=1:11
    if i==1
        regName = 'Head';
    else
        regName = ['PS',num2str(5+i-2)];
    end
    reg = strcmp(eTablesAbd.spotCellType,regName) & k;
    segMap{i} =  eMapsAbd(rs,rs,reg);
    segName{i} = regName;
    if pars.showPlot
        map = nanmedian(segMap{i} , 3);
        subplot(2,6,i); imagesc( map );
        title([regName,' n=',num2str(sum(reg))]);
        caxis([100,500]);
    end
end
if pars.showPlot
    colormap(flipud(parula));
    set(gcf,'color','w');
end
