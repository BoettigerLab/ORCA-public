function ratio = InsulationScore(a,varargin)
% Compute Insulation Score, compares intra-domain contacts (upstream or
% downstream) to inter-domain contacts (between the regions) on maps that
% have been distance normalized.  
% the output score is reported as a ratio, which is > 1 at boundaries.
% 
%  score = InsulationScore(map);
%  score = InsulationScore(map,'x',30,'w',5), returns the insulation score
%           at position boundary position 30 with a box of size 5 pixels
%          Note x is the last pixel in the upstream window and excluded
%          from the downstream window. 
% score = InsulationScore(map,'norm',false) 
%       uses the existing data to normalize distance effects first. 
% 
% See Also: InsulationWindow -- computes the difference between upstream
% and downstream average contacts.  
% 

defaults = cell(0,3);
defaults(end+1,:) = {'norm','boolean',true};
defaults(end+1,:) = {'method',{'median','mean'},'median'};
defaults(end+1,:) = {'dataType',{'contact','distance'},'contact'};
defaults(end+1,:) = {'w','integer',10}; % window
defaults(end+1,:) = {'showProgress','boolean',false};
defaults(end+1,:) = {'progressStep','integer',1e3};
defaults(end+1,:) = {'x','integer',0}; % points at which to compute.   
pars = ParseVariableArguments(varargin,defaults,mfilename);

% window size
w=pars.w;

% remove polymer distance effects 
if pars.norm
    norm = NormMap(a,'max',4*w+2);
    a = a./norm;
end

% position at which to compute score
if pars.x ==0
    N = size(a,1)-2*w;
    xs = 1:N;
else
    xs = pars.x;
    N = length(xs);
end
ratio = zeros(N,1);

k=0;
for x=xs
    k=k+1;
    b1 = x:x+w;
    b2 = x+w+1:x+2*w;
    L  = a(b1,b1);
    L = triu(L,1); 
    L(L==0) = nan;
    R = a(b2,b2);
    R = triu(R,1);
    R(R==0) = nan;
    X = a(b1,b2);
    X = tril(X,-1);
    X(X==0) = nan;
    if rem(k,pars.progressStep)==0 && pars.showProgress
        figure(10); clf; 
        subplot(2,3,1); imagesc(R); 
        subplot(2,3,2); imagesc(L);
        subplot(2,3,3); imagesc(X);
        subplot(2,3,4); imagesc(a([b1,b2],[b1,b2]));
    end
    if strcmp(pars.method,'median')
        inter = nanmedian(X(:));
        intra = [nanmedian(L(:)),nanmedian(R(:))];
        intra_max = max(nanmedian(L(:)),nanmedian(R(:)));
    elseif strcmp(pars.method,'mean')
        inter = nanmean(X(:));
        intra = [nanmean(L(:)),nanmean(R(:))];
        intra_max = max(nanmean(L(:)),nanmean(R(:)));
    end
%     if strcmp(pars.dataType,'contact')
%         intra_max = max(intra);
%     elseif strcmp(pars.dataType,'distance')
%         intra_max = min(intra);
%     end
    ratio(k) = intra_max./inter;
end