function allPts = UniformStepPoly(apol)
% adds squiggles between polymer steps so that an equal amount of path
% length is incorporated between each step

% isMiss = inan(apol(:,1));
apol = fillmissing(apol,'linear');

stps = sqrt(sum(diff(apol,[],1).^2,2));
targetLength = max(stps);
  stepSize = min(stps);

M = size(apol,1);
allPts = cell(M,1);
for m = 1:M-1
  a = apol(m,:);
  b = apol(m+1,:);
  ab = [a;b];
  gapLength = pdist(ab);

  ptsToAdd = ceil((targetLength-gapLength)/stepSize);
  d = b-a;
  o = null(d(:).'); % an orthogonal vector to b-a.  There are 2, we need only 1
  addDist = stepSize;
  gain = sqrt( ((addDist+gapLength)/2)^2-(gapLength/2)^2)/2;
  ps = linspace(0,1,2*ptsToAdd+1);
  mid = zeros(2*ptsToAdd,3);
  for p=1:length(ps)-1
      if rem(p,3)==0
        mid(p,:) = a + ps(p+1)*(b-a) + gain*o(:,1)';
      elseif rem(p,3)==1
          mid(p,:) = a + ps(p+1)*(b-a) + gain*o(:,2)';
      else
        mid(p,:) = a + ps(p+1)*(b-a);
      end
  end
  newPts = [a;mid;b];
%   newDist = cumsum(sqrt(sum(diff(newPts,[],1).^2,2)));
%   newDist-gap;
%   
  allPts{m} = newPts(1:end-1,:);

%   figure(11); clf; plot3(newPts(:,1),newPts(:,2),newPts(:,3));
%   figure(10); clf; plot(a(:,1),a(:,2),'.'); hold on;
%   plot(b(:,1),b(:,2),'.'); plot(mid(:,1),mid(:,2),'.');

end
allPts = cat(1,allPts{:});
% 
% figure(11); clf;
% subplot(2,1,1); hist(sqrt(sum(diff(allPts,[],1).^2,2)),0:20:1000);
% subplot(2,1,2); hist(sqrt(sum(diff(apol,[],1).^2,2)),0:20:1000);
% 
% figure(11); clf; plot3(allPts(:,1),allPts(:,2),allPts(:,3),'.-')
