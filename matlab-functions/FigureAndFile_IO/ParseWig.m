function [trk,xs2] = ParseWig(vals,locusTxt,chrNames)

[chr,st,en] = ParseLocusName(locusTxt);
xs = linspace(st,en,10);
i = StringFind(chrNames,chr,'exactly',true); 
in = vals{i}(:,1) > st & vals{i}(:,1) < en;
trk = vals{i}(in,2);
xs2 = linspace(st,en,length(trk));