function best = SortPassage(passageAligned,traceLength,varargin)
% 
% see AlignPassage
% 08/09/2024

defaults = cell(0,3); % 
% shared parameters for moving average
defaults(end+1,:) = {'d_start','positive',1}; % Distance used for start separation;  
defaults(end+1,:) = {'minPoints','positive',30}; % min observations in a passage event
defaults(end+1,:) = {'maxBackStepFrac','positive',0.2}; % max points in trace to record
defaults(end+1,:) = {'maxMiss','positive',0.2}; % max points in trace to record
pars = ParseVariableArguments(varargin,defaults,mfilename);
% pars = ParseVariableArguments([],defaults,mfilename);


[pa,ismiss] = fillmissing(passageAligned,"linear",2,"EndValues","none");
nMiss = sum(ismiss,2);
figure(1); clf; imagesc(passageAligned);

a = diff(pa,[],2);    
an = a; an(a>0) = nan; % back steps
ap = a; ap(a<0) = nan; % closer
b = abs(nansum(an,2));

[~,idx] = sort(b); % sort on least amount of backstep
pa2 = passageAligned(idx,:);

% sort on combination of low back step, sufficient length, low miss fraction   
best =find( nansum(ap,2) > pars.d_start & abs(nansum(an,2)) < .2*pars.d_start & nMiss./traceLength < pars.maxMiss & traceLength > pars.minPoints );

% figure(3); clf; plot(passageAligned(1:6,:)','.-'); 
% box off; 
% %         title(['processive motion, dis=',num2str(dis_kb1(e))]); ylabel('3D dist. (nm)');
% xlabel('frame (0.5s)')
% proc_trace =500-1+traceNum(best);
% proc_t_start= traceStart(best);
% proc_t_end = traceStart(best)+traceLength(best);
