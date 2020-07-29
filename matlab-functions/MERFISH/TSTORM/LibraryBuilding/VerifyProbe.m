function flag = VerifyProbe(probes,codebook,secondaries,varargin)
%--------------------------------------------------------------------------
% Generate Probes based on SECDEDcodewords
%
flag = 0;
%% Parse variable input
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'geneNamePosition', 'positive',5 };
defaults(end+1,:) = {'flagEndCommonName', 'string','__'};
defaults(end+1,:) = {'troubleshoot', 'boolean',false};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'A MList is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments([], defaults);

g = parameters.geneNamePosition-1; 
flagEndCommonName = parameters.flagEndCommonName;

%-------------------------------------------------------------------------
%% Main Function
%-------------------------------------------------------------------------
OligoSeq = {probes.Sequence};
OligoName = {probes.Header};


%Find out which mRNAs are detected in each round of hybridization
mRNAdetected = cell(length(secondaries),1);
for n = 1:length(secondaries)
    oligosDetected = ~cellfun(@isempty,strfind(OligoSeq,seqrcomplement(secondaries(n).Sequence)));
    sp = cellfun(@(x) strfind(x,' '),OligoName,'UniformOutput',false);
    un = cellfun(@(x) strfind(x,flagEndCommonName),OligoName,'UniformOutput',false);
    mRNAdetected{n} = unique(cellfun(@(sp,un,gname)  gname(sp(g)+1:un(1)-1),...
        sp(oligosDetected),un(oligosDetected),OligoName(oligosDetected),'UniformOutput',false));

    if parameters.troubleshoot % code good for troubleshooting errors  
        idx1 = sp(oligosDetected);
        idx2 = un(oligosDetected);
        mRNAnames = OligoName(oligosDetected);
        j = 1;
        mRNAnames{j}(idx1{j}(g)+1:idx2{j}(1)-1)
    end
end   

%Check result with the codebook
codes = false(length(codebook),length(str2num(codebook(1).Header)));
codeGenes = cell(length(codebook),1); 
for n = 1:length(codebook)
    codes(n,:) = logical(str2num(codebook(n).Header));
    sp = strfind(codebook(n).Sequence,' ');
    codeGenes{n} = codebook(n).Sequence(1:sp-1);
end    
for n = 1:length(str2num(codebook(1).Header)) 
    diffGenes = setdiff(codeGenes(codes(:,n)),mRNAdetected{n}');
    blanks = StringFind(diffGenes,'blank');
    diffGenes(blanks) = [];
    if isempty(diffGenes) ~= 1 
         display(['Hybridization round ' num2str(n) ' disagrees with codebook.']);
         display('codebook expected genes: '), disp(codeGenes(codes(:,n)));
         disp('secondaries detected genes: '), disp(mRNAdetected{n}');
         flag = 1;
    end
end  
end