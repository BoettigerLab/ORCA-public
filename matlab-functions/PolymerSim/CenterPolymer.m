function B = CenterPolymer(B,varargin)
% B = CenterPolymer(B);  center the polymer at zero
% B = CenterPolymer(B,x); center at a fixed corner x.   Example: x = (3,4,5)
% B = CenterPolymer(B,'min'); center in the positive octant. 
%
% updated to allow polymer stack for zero center


[nB,nD,nC] = size(B);
cm = nanmean(B);
% center at zero
if isempty(varargin)
   for c=1:nC
       cm = nanmean(B(:,:,c));
       B(:,:,c) = B(:,:,c) - repmat(cm,nB,1);
   end
% center in positive octant
elseif strcmp(varargin{1}, 'min')
   for d = 1:nD
        B(:,d) = B(:,d) - min(B(:,d));
   end
else % center at fixed corner
    corner = varargin{1};
    for d = 1:nD
        B(:,d) = B(:,d) - min(B(:,d)) + corner(d);
    end
end 

    


