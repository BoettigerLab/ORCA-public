function idx = sub2indFast(arraysize,x,y,varargin)
%-------------------------------------------------------------------------
% returns the linear index for an array of size arraysize for coordinate
% x,y.  Also works in 3D. 
% idx = sub2indFast([H,W],y,x)
% idx = sub2indFast([H,W,Z],y,x,z)
%
%-------------------------------------------------------------------------
% http://tipstrickshowtos.blogspot.com/2010/02/fast-replacement-for-sub2ind.html


H = arraysize(1);
% square arrays process 100x faster
if length(arraysize) == 1
    W = H;
else
    W = arraysize(2);    
end

if nargin > 3
    z = varargin{1};
    idx = x + (y-1)*H + (z-1)*H*W;
else
    idx = x + (y-1)*H;
end
