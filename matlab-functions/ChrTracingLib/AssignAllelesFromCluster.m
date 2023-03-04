function sortPts = AssignAllelesFromCluster(inputPts,varargin)
%  clusters 3D arrays of points into 2 groups by minimizing the intra-group
%  distance among all points which were initially detected in both groups.
%  Singlets are assigned uniquely to the group they are closest to.
% 
% 
%% Input
%   rawPts, a Nx4 matrix: x,y,z,positionID,
%% Output
%   sortPts, a Nx5 matrix: x,y,z,positionID,clusterLabel
%% Notes 
% 12/15/2020 - Alistair Boettiger
%  for future, may want to add a flag that allows singlets to be assigned
%  to both groups (treated as overlapping points).  This would be
%  espeically useful with the pairing data. 

%% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'maxMatch','integer',16};
defaults(end+1,:) = {'figScatter','integer',20};
defaults(end+1,:) = {'figMap','integer',21};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'remove',{'random','closest'},'random'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% Main function
try

% find points which occur at least twice
% (should do something about repeat data that occurs 4 times, or 3)
rawPts = [inputPts, (1:size(inputPts,1))']; % add a sorting index in col 5
rd = rawPts(:,4);
[v,n]= occurrences(rd);  % count occurrences
os = ismember(rd, v(n==2)');
jnt = rawPts(os,:);      % the points on both alleles
sngl = rawPts(~os,:);    % the points found on only 1 allele
tPts = size(jnt,1);
jnt = [jnt, (1:tPts)']; % add and index for sorting
% joint is now Nx6: x,y,z,name,origOrderID,jntListID
% sngl is now Nx5: x,y,z,name,origOrderID

% Downsample to avoid overload or freeze
%   This is a brute force algorithm over all possible clusterings
later = [];
totPairs = ceil(tPts/2);
if totPairs > pars.maxMatch
    if pars.verbose
       warning(['Requested assignment with ',num2str(totPairs),' matched points' ...
           ' will require ',num2str(2^totPairs) ' brute force combinations.  '...
           'Downsampling to ',num2str(pars.maxMatch),' matched points.'])
    end
   % need to remove points in pairs
   nRemove = ceil(totPairs - pars.maxMatch);
   [val,idx] = unique(jnt(:,4),'sorted');
   % randomly remove
   if strcmp(pars.remove,'random')
       idxRem = randperm(length(val));
       val = val(idxRem);
       selRemove = ismember(jnt(:,4),val(1:nRemove));
   % remove least separated
   elseif strcmp(pars.remove,'closest')
       allDist = squareform(pdist(jnt(idx,1:3)));
       [~,idxRem] = sort(sum(allDist,1));
       val = val(idxRem);
       selRemove =  ismember(jnt(:,4),val(1:nRemove));
   end
   % save the rest for later. we'll still process them as pairs
   later = jnt(selRemove,:);
   jnt = jnt(~selRemove,:);
   tPts = size(jnt,1);
   jnt(:,6)  =(1:tPts)'; % add and index for sorting
end

% % plot current distance map (mostly for troubleshooting) 
% mp = squareform(pdist(jnt(:,1:3)));
% figure(2); clf; imagesc(mp);

% create two matching lists, sorted by ID-number
%   each list has a full set of points (i.e. a potential contig chromosome)
[val,id1] = unique(jnt(:,4),'sorted');
is1 = false(tPts,1);
is1(id1) = true;
list2 = jnt(~is1,6);
[~,i] = sort( jnt(~is1,4) );
id2 = list2(i);

% compute all shuffles between these two lists
numDataBits = length(val);
codebook = boolean( rem(floor([0:2^numDataBits-1]'*pow2(-(numDataBits-1):0)),2) ); %#ok<NBRAK>
nCs = size(codebook,1); % total number of combos
b1 = repmat(id1',nCs,1);
b2 = repmat(id2',nCs,1);
combos = nan(size(codebook));
combos(codebook) = b1(codebook);
combos(~codebook) = b2(~codebook);
aveDist = zeros(nCs,1);
for n=1:nCs   % this is slow since it is an enormous loop. try Cost functions like in my live trace code? 
    k1 = combos(n,:);
    xyz1 = jnt(k1,1:3);
    k2 = 1:tPts;
    k2(k1) = [];
    xyz2 = jnt(k2,1:3);
    dlist = cat(2,pdist(xyz1),pdist(xyz2));
    aveDist(n) = nanmean(dlist);
end
[~,i] = min( aveDist );
k1 = combos(i,:);
k2 = 1:size(jnt,1);
k2(k1) = [];
cent1 = nanmean(jnt(k1,1:3),1);
cent2 = nanmean(jnt(k2,1:3),1);

% Assign single points to the nearest centroid
nS = size(sngl,1);
grp = zeros(nS,2);
for t=1:nS
    if grp(t,1)==0
        [grp(t,1),grp(t,2)] = knnsearch([cent1; cent2],sngl(t,1:3));
    end
end

% address postponed points
if ~isempty(later)
    later(:,6) = 0;
    val = unique(later(:,4),'sorted');
    for v=1:length(val)
        i = find(ismember(later(:,4),val(v)));
        test = later(i,1:3);
        aa = sum((test(1,:) - cent1).^2);
        ab = sum((test(1,:) - cent2).^2);
        ba = sum((test(2,:) - cent1).^2);
        bb = sum((test(2,:) - cent2).^2);
        [~,o] = min([aa+bb,ab+ba]);
        if o==1
            later(i(1),end) = 1;
            later(i(2),end) = 2;
        else
            later(i(1),end) = 2;
            later(i(2),end) = 1;
        end
        later(i(3:end),end) = 0;
    end
    l1 = later(later(:,end)==1,:);
    l2 = later(later(:,end)==2,:);
else
    l1 = []; l2 = [];
end
        

% combine data into a sorted list, append cluster assignment
c1 = [jnt(k1,1:5), 1*ones(length(k1),1)];
c2 = [jnt(k2,1:5), 2*ones(length(k2),1)];
s1 = [sngl(grp(:,1)==1,:), 1*ones(sum(grp(:,1)==1),1)];
s2 = [sngl(grp(:,1)==2,:), 2*ones(sum(grp(:,1)==2),1)];
sortPts = [c1; l1; s1; c2; l2; s2]; 

if pars.figMap
       figure(pars.figMap); 
       subplot(1,2,2);
       imagesc( squareform( pdist( sortPts(:,1:3) ) ) );  colorbar;
       title('sorted');
end

[~,i] = sort(sortPts(:,5)); % sort back to original order, keep group assingments
sortPts = sortPts(i,:);

if pars.figScatter
       figure(pars.figScatter); 
       r = sortPts(:,end)==1;
       g = sortPts(:,end)==2;
       b = sortPts(:,end)==0;
       plot3(sortPts(r,1),sortPts(r,2),sortPts(r,3),'r.'); hold on;
       plot3(sortPts(g,1),sortPts(g,2),sortPts(g,3),'g.'); hold on;
       text(sortPts(r,1),sortPts(r,2),sortPts(r,3),cellstr(num2str(sortPts(r,4))),'color','r');
       text(sortPts(g,1),sortPts(g,2),sortPts(g,3),cellstr(num2str(sortPts(g,4))),'color','g');
       plot3(c1(:,1),c1(:,2),c1(:,3),'ro'); hold on;
       plot3(c2(:,1),c2(:,2),c2(:,3),'go'); 
       plot3(sortPts(b,1),sortPts(b,2),sortPts(b,3),'b.'); hold on;
       hold off;
       pause(.1);
end

if pars.figMap
       figure(pars.figMap); 
       subplot(1,2,1);
       imagesc( squareform( pdist( sortPts(:,1:3) ) ) ); colorbar;
       title('original');
       pause(.1); 
end

catch er
    warning(er.getReport);
    disp('place debug here');
    error(er.getReport);
end
