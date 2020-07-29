function zshift = ComputeHybeZoffset(rawDataNames,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'cropFrames','integer',14:150}; % this is a bit arbitrary
defaults(end+1,:) = {'troubleshoot','boolean',false}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);

% =============  check zoffset

   showplots = pars.troubleshoot;
   numHybes = size(rawDataNames,1);
   warning('off','MATLAB:table:ModifiedAndSavedVarnames');
   offName = regexprep(rawDataNames{1},'.dax','.off');
   off1 = ReadTableFile(offName,'preview',false);
   off1 = off1(pars.cropFrames,:);
   zshift = zeros(numHybes,2);
   for h=1:numHybes
       offName = regexprep(rawDataNames{h},'.dax','.off');
       offH = ReadTableFile(offName,'preview',false);
       offH = offH(pars.cropFrames,:);
       if showplots
           figure(1); clf; 
           plot(off1.offset,'color','b'); hold on;
           plot(offH.offset,'color','r');
       end
        shiftUp = true; %  median(offH.offset) < median(off1.offset);
        nT = length(off1.offset)-10;
        testZs = zeros(nT,1);
        for z=1:nT
            if shiftUp
                testZs(z) = sum(off1.offset(z:end) - offH.offset(1:end-z+1));
            else
                testZs(z) = sum(off1.offset(1:end-z+1) - offH.offset(z:end));
            end
        end
        if showplots
            figure(2); clf; plot(abs(testZs));
        end
        [v,i] = min(abs(testZs));
        if shiftUp
           zshift(h,1) = ceil(i/2)-1;
           zshift(h,2) = -v;
           xx = 1+i:length(off1.offset)+i;
           yy = offH.offset;
        else
            zshift(h,1) = ceil(i/2)-1;
            zshift(h,2) = v;
            xx = 1:length(offH.offset(i:end));
            yy = offH.offset(i:end);
        end
        if showplots
            figure(1); 
            plot(xx,yy,'m--');
            legend('1','h','corrected');
            pause(1);
        end
   end
