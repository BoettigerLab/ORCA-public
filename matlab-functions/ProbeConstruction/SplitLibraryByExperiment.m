function probeLibs = SplitLibraryByExperiment(probeLib,varargin)
% probe library can be either a probe library fasta-file struct

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'delim','string','__'};
defaults(end+1,:) = {'itemProbeName','integer',4};
defaults(end+1,:) = {'writeFasta','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);


if ischar(probeLib)
    probeLibName = probeLib;
    probeLib = fastaread(probeLib);
elseif pars.writeFasta
    warning('writeFasta requires probeLib be a filepath to a fastafile');
    pars.writeFasta = false; 
end

nameparts = cellfun(@(x) strsplit(x,pars.delim),{probeLib.Header},'UniformOutput',false);
libNameParts = cat(1,nameparts{:});
[revPrimers,uIdx] = unique(libNameParts(:,end),'stable');

if pars.verbose
    disp('Library contains revPrimers:');
    disp(revPrimers);
    disp('Library experiments include probeNames:');
    disp(libNameParts(uIdx,pars.itemProbeName));
end

expNames = cellfun(@(x) strsplit(x,'_'), libNameParts(uIdx,pars.itemProbeName),'UniformOutput',false);
expNames = cellfun(@(x) x{1},expNames,'UniformOutput',false);

totExps = length(revPrimers);
totProbes = size(libNameParts,1);
uIdx = [uIdx; totProbes];
probeLibs = cell(totExps,1);
for e=1:totExps
    probeLibs{e} = probeLib(uIdx(e):uIdx(e+1)-1);
end
 
if pars.writeFasta
    for e=1:totExps
        try
        expName = regexprep(probeLibName,'.fasta',['_E',num2str(e),'_',expNames{e},'.fasta']);
        WriteFasta(expName,probeLibs{e},[],'Append',false);
        probeLibs{e} = expName; 
        catch er
            warning(er.message);
        end
    end
end