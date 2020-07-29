function WriteFasta(filename,Header,seq,varargin)
%  WriteFasta(filename,Header,seq,varargin)
%  WriteFasta(filename,GeneFasta,[],'Append',false,'Warnings',true)
% GeneFasta may be a Nx1 structure array or a single structure containing
% two cell arrays.  The fields must be '.Header' and '.Sequence'.

%% Default Parameters
appendFlag = false;
WarningsOn = true; 

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if nargin > 3
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'Append'
                appendFlag = CheckParameter(parameterValue,'boolean','Append');
            case 'Warnings'
                WarningsOn = CheckParameter(parameterValue,'boolean','Warnings');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end



%% Main Function

% bioinfochecknargin(nargin,2,mfilename)  <-- This is unnecessary


if ~ischar(filename)
    error(message('bioinfo:fastawrite:FilenameMustBeString'));
end

% The default behavior is to overwrite the file.
% If the file exists, deliver a warning that data is being appended to
% something that exists or overwriting something that exists.  


if appendFlag
    if exist(filename,'file') && WarningsOn
       warning(['data will be appended to existing file ',filename]);  
    end
    fid = fopen(filename,'at');
else
    if exist(filename,'file') && WarningsOn
       warning(['existing file ',filename,' will now be overwritten']);  
    end
    fid = fopen(filename,'w+');
end

if fid == (-1)
    [theDir, theFile, theExtension] = fileparts(filename);
    if ~isempty(theDir)
        error(message('bioinfo:fastawrite:CouldNotOpenFileinDir', [ theFile, theExtension ], theDir));
    else
        error(message('bioinfo:fastawrite:CouldNotOpenFileinPwd', filename));
    end
end


try

    % Parse do we have a structure or cell array
    if isstruct(Header)
        try 
            if length(Header)> 1
                header = {Header.Header};
                seq = {Header.Sequence};
            else
                header = Header.Header;
                seq = Header.Sequence;
            end
        catch  
            disp(['if passed a structure WriteFasta expects two fields, ',...
                '.Header and .Sequence ']); 
        end
    elseif ischar(Header) % single character array
        header = {Header};
        seq = {seq}; 
    elseif iscell(Header)
        header = Header;
    else 
        error('Header is not a character, cell, or field of structure'); 
    end

    numSequences = length(header); 

    for i=1:numSequences
        currseq = seq{i};
        currheader = header{i};
        len = size(currseq,2);
        if len == 0
            error(message('bioinfo:fastawrite:EmptySeq'));
        end
        % Add the > token if needed
        if isempty(currheader) || currheader(1) ~= '>'
            fprintf(fid,'>%s\n',currheader);
        else
            fprintf(fid,'%s\n',currheader);
        end
        maxcols = 70; % NCBI uses 70 columns for their FASTA files
        for line = 1:ceil(len/maxcols)
            start = ((line - 1) * maxcols) + 1;
            stop = min((line * maxcols),len);
            fprintf(fid,'%s\n',currseq(start:stop));
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
catch le
    % close
	fclose(fid);
	
	% delete only if it was a newly created file
	if ~appendFlag
		delete(filename);
	end

    %rethrow the error
    rethrow(le);
end
