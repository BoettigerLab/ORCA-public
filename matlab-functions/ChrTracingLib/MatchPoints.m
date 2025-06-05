function [matched,cost,costDif] = MatchPoints(xyz1,xyz2,varargin)
% Pair N points
% input xyz1,xyz2 or xy1 xy2 (2D and 3D inputs allowed, pdist sorts it
% out). 
% The points will be sorted such that the second 
% 
% output matched
% matched, a 2 x N matrix listing all the matched pairs between the lists
% costDif, the difference in cost between the best solution 
%  
% 
% example
% matched = MatchPoints(p1,p2); 
% plot([p1(matched(:,1),1),p2(matched(:,2),1)]',...
%     [p1(matched(:,1),2),p2(matched(:,2),2)]','k-');
% hold on;
% plot(p1(matched(:,1),1),p1(matched(:,1),2),'ro'); hold on;
% plot(p2(matched(:,2),1),p2(matched(:,2),2),'b+');


defaults = cell(0,3);
defaults(end+1,:) = {'rank','integer',1};
defaults(end+1,:) = {'maxDist','nonnegative',0};
pars = ParseVariableArguments(varargin,defaults,mfilename);



n1 = size(xyz1,1);
n2 = size(xyz2,1);

% this method uses matlabs function 'matchpairs'
%    this doesn't require equal number of points
% fast computation of all pairwise distances
mp = squareform(pdist([xyz1;xyz2]));
if pars.maxDist ~=0
    mp(mp>pars.maxDist) = pars.maxDist;
end
mp = mp(1:n1,n1+1:n1+n2);
costUnmatched = sum(mp(:)); % can't be infinite, or all unmatched soln are equal. 
matched = matchpairs(mp,costUnmatched); % minimize cost (distance) matrix

cost = diag(mp(matched(:,1),matched(:,2)));
costDif = [];

if nargout > 2 && n1 == n2
    n=n1;
    % % compute total distance of all assignments
    % %   A valid assignment takes 1 row-position from each column of mp. 
    per = perms(1:n);
    dis = nan(size(per,1),1);
    for p=1:size(per,1)
        dis(p) = sum(diag(mp(1:n,per(p,:))));
    end
    [vs,is] = sort(dis,'ascend'); % in place of 'min', so we can compare next nearest 
    i = is(pars.rank);
    costDif = vs(pars.rank+1)-vs(pars.rank);
    cost = diag(mp(1:n,per(i,:))); %sum this get: vs(1);
    matched = [1:n; per(i,:)]';
end

