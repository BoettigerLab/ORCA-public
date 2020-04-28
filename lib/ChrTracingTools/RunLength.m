function [runsON,runsOFF] = RunLength(map1,varargin)
% takes a matrix of runs 1s and 0s, each row is a different trial
% returns a list of the length of all the runs of 1s 
% Optionally, RunLength(map1,'type','first') returns only the first switch
% from 1 to 0 - i.e. the length of the initial run for each column.
% 

% default pars
defaults = cell(0,3);
defaults(end+1,:) = {'type',{'all','first'},'all'};
defaults(end+1,:) = {'troubleshoot','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% for troubleshooting
if pars.troubleshoot
    figure(1); clf; imagesc(map1);
end

% main function
map2 = diff(map1,1,2);
[nRows,nCols] = size(map1);
if pars.troubleshoot
    figure(2); clf; imagesc(map2);
end

try
switch pars.type
    case 'all'       
        runsON = cell(nRows,1);
        runsOFF = cell(nRows,1);
        for n=1:nRows % n=2
            run10 = find(map2(n,:)==-1)+1; % switched 1->0
            run01 = find(map2(n,:)==1)+1; % switched 0->1
            if isempty(run10) && isempty(run01) % no switches
               if map1(n,1) == 1
                   runsON{n} = nCols;
                   runsOFF{n} = 0;
               else
                   runsON{n} = 0;
                   runsOFF{n} = nCols;
               end
            else
                if isempty(run10)
                    run10 = nCols+1;
                end
                if isempty(run01)
                    run01 = nCols + 1;
                end
                % even number of switches, first switch is ON
                if length(run10) == length(run01)
                    if run10(1) < run01(1)
                        runsON{n} = [run10,nCols+1] - [1,run01];
                        runsOFF{n} = run01 - run10;
                    else
                        runsON{n} = run10 - run01; 
                        runsOFF{n} = [run01,nCols+1] - [1,run10];
                    end
                end
                % 
                if run10(1) < run01(1) && length(run10) > length(run01) % switched 1->0 first
                    runsON{n} = run10 - [1,run01];
                    runsOFF{n} = [run01,nCols+1] - run10;
                elseif  run10(1) > run01(1)  && length(run01) > length(run10) % switched 0->1 first
                    runsON{n} = [run10,nCols+1] - run01;
                    runsOFF{n} = run01 - [1,run10] ;
                end
            end
        end
        runsON = cat(2,runsON{:});
        runsOFF = cat(2,runsOFF{:});
        runsON(runsON==0) = [];
        runsOFF(runsOFF==0) = [];
    case 'first'
        % find first occurrence in each row
        runsON = zeros(nRows,1);
        runsOFF = zeros(nRows,1);
        for n=1:nRows
            runsON(n) = find(map2(n,:)==-1,1);
            runsOFF(n) = find(map2(n,:)==1,1);
        end
end

catch er
    disp(er.getReport);
    disp('er here');
end