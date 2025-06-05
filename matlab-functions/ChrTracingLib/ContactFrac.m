function [contFreq, nObs] = ContactFrac(distMap,varargin)
% [contFreq,nObs] = ContactFrac(distMap);
%    returns the contact frequency map, contFreq, and the matrix of the
%    number of observations for each pair, nObs. 
% [contFreq,nObs] = ContactFrac(distMap,'threshold',250);
%    points closer than the threshold 250 nm are declared in contact. 
% [contFreq,nObs] = ContactFrac(distMap,'relativeDist',.5); 
%    points closer than half the median distance among all points are 
%    declared in contact 

%% Defaults
defaults = cell(0,3);
defaults(end+1,:) = {'relativeDist','nonnegative',.7}; % relative to median distance between all regions in spot.  0 for off.
defaults(end+1,:) = {'threshold','nonnegative',150}; % absolute distance in nm.  
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% Main Function
[nReads,~,nSpots] = size(distMap);
if pars.threshold > 0
    contCount = sum(distMap < pars.threshold,3);
else
    medDist = nanmedian(reshape(distMap,nReads^2,nSpots),1);
    contCount = zeros(nReads,nReads,nSpots);
    for s=1:nSpots
        contCount(:,:,s) = distMap(:,:,s) < pars.relativeDist*medDist(s);
    end
end
nObs = nSpots - sum( isnan(distMap),3);
contFreq = sum(contCount,3)./nObs;


% convert nobs to p-value
%  xi squared or something, there is a contact, no-contact, and N
% 
