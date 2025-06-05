function [ksd,p_val] = Reversibility(distTraceMat,varargin)
% inputs 
% distTraceMat  nT, T
% see also

defaults = cell(0,3);
% data loading
defaults(end+1,:) = {'tau','integer',3}; 
defaults(end+1,:) = {'showPlots','integer',1}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);

tau = pars.tau;


[nT,T] = size(distTraceMat);
combos = rem(floor([0:2^tau-1]'*pow2(-(tau-1):0)),2);
combos(combos==0)=-1;
n_fwd = zeros(1,2^tau); 
n_rev = zeros(1,2^tau);
for n=1:nT
    distTrace = distTraceMat(n,:);
    vdir = nan(T,tau);
    distRev = fliplr(distTrace);
    vrev = nan(T,tau); 
    for t=1:T-tau
        v = distTrace(t:t+tau);
        vdir(t,diff(v)>0) = 1; % to handle missing data
        vdir(t,diff(v)<0) = -1;
        v = distRev(t:t+tau);
        vrev(t,diff(v)>0) = 1;
        vrev(t,diff(v)<0) = -1;
    end
    drop = isnan(vdir(:,1));
    vdir(drop,:) = [];
    drop = isnan(vrev(:,1));
    vrev(drop,:) = [];
    % count occurrences
     %  data might not exhibit all possible combinations.  
    [nC,tau] = size(combos);
    nfwd = zeros(1,nC);
    nrev = zeros(1,nC);
    for c=1:nC
        ia = ismember(vdir,combos(c,:),'rows');
        nfwd(c) = sum(ia);
        ia = ismember(vrev,combos(c,:),'rows');
        nrev(c) = sum(ia);
    end
    n_fwd = n_fwd + nfwd;
    n_rev = n_rev + nrev;
end
p_fwd = n_fwd/(sum(n_fwd));
p_rev = n_rev/(sum(n_rev));
ksd = sum( p_fwd.*log(p_fwd./p_rev));

p_val = zeros(1,length(n_fwd));
for i=1:length(n_fwd)
    p_val(i) = myBinomTest(n_fwd(i),n_fwd(i)+n_rev(i),0.5,'two'); % obs, total, p, two-sided
end

if pars.showPlots
    figure(pars.showPlots); clf; 
    subplot(2,1,1); bar(-log10(p_val)); ylabel('-log10 p-value'); xlim([0,length(p_val)+1]); title('multi-step direction persistence')
    ym = max([max(-log10(p_val)),2]); ylim([0,ym])
    subplot(2,1,2); bar([n_fwd;n_rev]'); ylabel('n-Obs'); xlim([0,length(p_val)+1])
end