function allhits = MaxBLASThits(BLASTdata,varargin)
%--------------------------------------------------------------------------
% allhits = MaxBLASThits(BLASTdata)
% allhits = MaxBLASThits('Path to Blast Data')
% Returns vector with the maximum number of matching bases for each member
% in the library. 
% 
%--------------------------------------------------------------------------
%% Required Inputs
% BLASTdata - string to BLAST output file
%             or matlab BLAST datastructure from 'blastreadlocal'
%
%--------------------------------------------------------------------------
%% Optional Inputs
% 'maxhits' / double / 14 
%                       -- max hits allowed to keep
% 'no3pMatch' / logical / false
%                        -- sets hit count to maxhits + 1 if the 3' region
%                       of the query matches the target.  Used to make 
%                       selective primers .
% 'verbose' / logical / false
%                       -- prints hits, highlights hits over maxhits using
%                       the warning function. 
% 'exludeFirstHit' / logical /false
%                       -- set true to skip first hit (e.g. if first hit is
%                       expected to be self).  
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Default Parameters
%--------------------------------------------------------------------------
exludeFirstHit = false;
no3pMatch = false;
verbose = false;
maxhits = 14;
primers = '';

%--------------------------------------------------------------------------
% Parse mustHave variables
%--------------------------------------------------------------------------
if nargin < 1
   error([mfilename,' expects 1 input, Blast Data structure or Data filename']);
end


%--------------------------------------------------------------------------
% Parse variable input
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
            case 'maxhits'
                maxhits = CheckParameter(parameterValue,'nonnegative','maxhits');
            case 'no3pMatch'
                no3pMatch = CheckParameter(parameterValue,'boolean','no3pMatch');
            case 'primers'
                primers = parameterValue; % a structure or an empty string
            case 'exludeFirstHit'
                exludeFirstHit = CheckParameter(parameterValue,'boolean','exludeFirstHit');
            case 'verbose'
                verbose = CheckParameter(parameterValue,'boolean','verbose');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

if ischar(BLASTdata) 
    Data = blastreadlocal(BLASTdata, 0);
else
    Data = BLASTdata;
end

if no3pMatch
    if isempty(primers)
        error('parameter no3pMatch requires "primers" option'); 
    end
    primerlength = length(primers(1).Sequence);
end

if exludeFirstHit
    hitstart = 2;
else
    hitstart = 1;
end

%--------------------------------------------------------------------------
%% Main Function
%--------------------------------------------------------------------------

N = length(Data);
allhits = zeros(N,1); 
for n=1:N
    H = length(Data(n).Hits);
    lenhit = zeros(H,1);
    for h=hitstart:H % loop over hits
        HSPs = length( Data(n).Hits(h).HSPs );
        tempdat = zeros(HSPs,1);
        for j = 1:HSPs % loop over HSPs
            tempdat(j) = sum(Data(n).Hits(h).HSPs(1).Alignment(2,:)=='|');
            
            if no3pMatch && ~isempty(primers);
                % Reject sequences that have a perfect match at the 3' base.  
                try
                [~,e] = regexp(primers(n).Sequence,Data(n).Hits(h).HSPs(1).Alignment(1,:) );
                if isempty(e)
                    e=0;
                end
                catch
                    e=0;
                end
                if tempdat(j) > maxhits-1 && e == primerlength
                    tempdat(j) = maxhits + 1; 
                end
            end
            
        end
        lenhit(h) = max(tempdat);  
    end
    lenhit = max(lenhit); 
    if ~isempty(lenhit)
        allhits(n) = lenhit;
    end
    if n<5 && verbose
        disp([Data(n).Query,' has largest hits of ',num2str(lenhit)]);      
    end
    if verbose
        if lenhit > maxhits
                disp([Data(n).Query,' has hits of length ',num2str(lenhit)]);
        end
        if lenhit <= maxhits
            disp([Data(n).Query,' has largest hits of ',num2str(lenhit)]);
        end  
    end
end