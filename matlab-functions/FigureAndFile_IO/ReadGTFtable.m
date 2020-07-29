function gtfTable = ReadGTFtable(gtfFile,varargin)

% -------- Load GTF Data -------------------------
 disp('Reading genedata GTF file...'); 
fid = fopen(gtfFile);
fmt = repmat('%s ',1,11);
rawData = textscan(fid,fmt,'delimiter','\t','HeaderLines',5);
disp('geneData loading complete');
%-----------------------------------------------------
%  GTF columns
% 'chr'     1
% 'source'      2
% 'feature'  3
% 'start'         4
% 'stop' 5
% 'dot1'          6
% 'strand'           7
% 'dot2'          8
% 'geneID'        9

chr = rawData{1};
feature = rawData{3};
start = rawData{4};
stop = rawData{5};
strand = rawData{7};
name = rawData{9};
gtfTable = table(chr,feature,start,stop,strand,name);
