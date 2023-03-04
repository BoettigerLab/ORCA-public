function ProbeFastaToBed(libName,varargin)
% take a library fasta file, build a bed file

defaults = cell(0,3);
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'verbose','boolean',false};
defaults(end+1,:) = {'blastDatabase','string','U:\GenomeData\GenomeAssemblies\hg19\hg19.fasta'};
defaults(end+1,:) = {'blastWordSize','integer',35}; 
pars = ParseVariableArguments(varargin,defaults,mfilename); 

%% Load data
% saveFolder = 'C:\Data\Oligos\Library03\';
% libName = [saveFolder,'RNAprobes\','chr21RnaProbes.fasta'];

%% Load probes
[sourceFolder,fname] = fileparts(libName);
if isempty(pars.saveFolder)
    saveFolder = [sourceFolder,filesep];
else
    saveFolder = pars.saveFolder;
end
libProbes = fastaread(libName);


%% Run BLAST ( a bit slow on hg)
[~,genome] = fileparts(pars.blastDatabase);
[~,blastData] = BLAST(libName,pars.blastDatabase,'numThreads',40,'wordSize',pars.blastWordSize,'verbose',pars.verbose);

%% write table

% Get Colormap (based on number of reads)
% nameparts = strsplit(blastData(1).Query,'__');
% read = nameparts{3};
% firstRead = str2double(regexprep(read,{'Read','read','Stv_'},{'','',''}));
% 
% nameparts = strsplit(blastData(end).Query,'__');
% read = nameparts{3};
% lastRead = str2double(regexprep(read,{'Read','read','Stv_'},{'','',''}));
% numReads = lastRead-firstRead+1;


nameparts = cellfun(@(x) strsplit(x,'__'),{blastData.Query},'UniformOutput',false);

namebits = cat(1,nameparts{:});
if length(nameparts{1})==5 % with commonRT
    readNames = namebits(:,3);
else
    readNames = namebits(:,2);
end
% remove 'was188' from readout053was188
probeSwap = cellfun(@(x) strfind(x,'was'), readNames,'UniformOutput',false);
for i=1:length(probeSwap)
    if ~isempty(probeSwap{i})
        readNames{i} = readNames{i}(1:probeSwap{i}-1);
    end
end
 % 
readNums =  cellfun(@(x) str2double(regexprep(x,{'Fiducial','Read','readout','barcode','read','Stv_','B188','B189','B'},{'00','','','','','','','',''})),readNames);  
firstRead = min(readNums);
lastRead = max(readNums); 
numReads = lastRead-firstRead+1;
readInds = readNums-firstRead+1;

try
    readColor = round(255*hsv(  round(numReads*1.1))  ); 
    readColor = readColor(1:numReads,:);
catch
    error('problem creating colormap');
end
% readColor = readColor( randperm(numReads),:);

% other useful stuff
numProbes = length(libProbes);
% read = cell(numProbes,1);

% bed table parts
chromStart = zeros(numProbes,1);
chromEnd = zeros(numProbes,1);
name = cell(numProbes,1);
strand  = cell(numProbes,1);
itemRgb = cell(numProbes,1); 
chrom = cell(numProbes,1); 
for  i=1:length(libProbes)  
    if ~isempty(blastData(i).Hits)
        try
        hitName = strsplit(blastData(i).Hits(1).Name,'|');
        chr = hitName{2};
        chrom{i} = chr;
        coords = blastData(i).Hits(1).HSPs.SubjectIndices;
        % nameparts = strsplit(blastData(i).Query,'__');
        % read{i} = nameparts{3};
        % r = str2double(regexprep(read{i},{'Read','read','Stv_','B188','B189'},{'','','','',''})) - firstRead + 1; 
        r = readInds(i);
        if coords(2) > coords(1)
            chromStart(i) = coords(1);
            chromEnd(i) = coords(2);    
            strand{i} = '+';
        else
            chromStart(i) = coords(2);
            chromEnd(i) = coords(1);
            strand{i} = '-';
        end
        probeFullName = regexprep(blastData(i).Query, {'\.',' '},{'',''});
        nameparts = strsplit(probeFullName,'__');
        if length(nameparts) == 5
            targetNameParts = strsplit(nameparts{4},'_');
            nameShort = [nameparts{3},'_',targetNameParts{end}];
        else
            nameShort = probeFullName;
        end
        name{i} = nameShort;
        itemRgb{i} = [num2str(readColor(r,1)),',',num2str(readColor(r,2)),',',num2str(readColor(r,3))];
        catch er
            warning(er.getReport);
            disp('place debug here'); 
            error('problem here');
        end
    end
end
remove = chromStart==0;

chrom(remove) = [];
chromStart(remove) = [];
chromEnd(remove) = [];
name(remove) = [];
strand(remove) = [];
itemRgb(remove) = [];

numRecords = length(chromStart);
score = ones(numRecords,1); 
thickStart = chromStart;
thickEnd = chromEnd;
    
libBed = table(chrom,chromStart,chromEnd,name,score,strand,thickStart,thickEnd,itemRgb);

%% 
header = ['track name="',fname,'" description="',genome,' ',fname,'" visibility=1 itemRgb="On"'];
bedfile =[saveFolder,fname,'_bed.txt']; 
headerfile = [saveFolder,fname,'_header.txt'];
datafile = [saveFolder,fname,'_data.txt'];
if exist(bedfile,'file')~=0
    delete(bedfile);
end

fid = fopen(headerfile,'w+');
fprintf(fid,'%s\r\n',header);
fclose(fid); 
writetable(libBed,datafile,'delimiter','\t','WriteVariableNames',false);

[~,~] = system(['copy /b ',headerfile,'+',datafile,' ',bedfile]);
delete(datafile);
delete(headerfile);

