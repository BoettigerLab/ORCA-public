function [libExpect,libExpectPerHybe,libGenes,libCodes] = ExpectedLocsPerHybe(codebook,varargin)

%% globals
global fpkmData

% Defaults
fpkmPath = '\\tuck\tstorm\GenomeData\A549\';
codebookPath = '';

%--------------------------------------------------------------------------
%% Parse variable input
%--------------------------------------------------------------------------
if nargin > 1
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'fpkmPath'
                fpkmPath = CheckParameter(parameterValue,'string','fpkmPath');
            case 'codebookPath'
                codebookPath = CheckParameter(parameterValue,'string','codebookPath');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%%

if  exist([codebookPath,'libExpect.mat'],'file')
    load([codebookPath,'libExpect.mat']);
else

    % Extract gene names from codebook
    if length(codebook) < 15
        libGenes = {codebook.Sequence};
        numHybes = 8;
    else
        lib4names = {codebook.Sequence}';
        gaps = cellfun(@(x) strfind(x,' '),lib4names,'UniformOutput',false);
        libGenes = cellfun(@(x,y) x(1:y(1)-1),lib4names,gaps,'UniformOutput',false);
        numHybes = 16;
    end


    % -------- Load A549 Data -------------------------
    if isempty(fpkmData)
    disp('Reading Cufflinks Gene Data...'); 
    dfile = 'isoforms.fpkm_tracking';
    fid = fopen([fpkmPath,dfile]);
        textscan(fid,'%s',13,'delimiter','\t'); % skip header;
        fmt = repmat('%s ',1,13);
        fpkmRaw = textscan(fid,fmt,'CollectOutput',true,'TreatAsEmpty','-');
    fclose(fid);
    fpkmRaw = fpkmRaw{1}; 
    fpkmData.fpkm = str2double(fpkmRaw(:,10));
    fpkmData.names = fpkmRaw(:,5);
    end
    %-----------------------------------------------------
    % 'tracking_id'     1
    % 'class_code'      2
    % 'nearest_ref_id'  3
    % 'gene_id'         4
    % 'gene_short_name' 5
    % 'tss_id'          6
    % 'locus'           7
    % 'length'          8
    % 'coverage'        9
    % 'FPKM'            10
    % 'FPKM_conf_lo'    11
    % 'FPKM_conf_hi'    12
    % 'FPKM_status'     13
    %-----------------------------------------------------


    % Compute expected localizations per hybe
    libExpect = zeros(length(libGenes),1);
    for i=1:length(libGenes)
        idx =  find(~cellfun(@isempty,strfind(fpkmData.names, libGenes{i})));
        mtch = cellfun(@length,fpkmData.names(idx)) == length(libGenes{i});
        idx = idx(mtch); 
        libExpect(i) = sum(fpkmData.fpkm(idx));
    end

    libCodes = logical(cell2mat(cellfun(@str2num, {codebook.Header},'UniformOutput',false)'));
    libExpectPerHybe = zeros(numHybes,1);
    for i=1:numHybes % 
        libExpectPerHybe(i) = sum(libExpect(libCodes(:,i)));
    end
    libExpectPerHybe = libExpectPerHybe/sum(libExpectPerHybe);


    if ~isempty(codebookPath)
       save([codebookPath,'libExpect.mat'],'libExpect','libExpectPerHybe','libGenes'); 
    end

end