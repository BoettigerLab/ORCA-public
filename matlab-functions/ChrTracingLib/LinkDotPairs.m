function [xdat,ydat,zdat] = LinkDotPairs(dots1,dots2,varargin)

[numPts,dims] = size(dots1);

try
    xdat = reshape([dots1(:,1) dots2(:,1) NaN*ones(numPts,1)]',numPts*3,1);
    ydat = reshape([dots1(:,2) dots2(:,2) NaN*ones(numPts,1)]',numPts*3,1);
    if dims == 3
        zdat = reshape([dots1(:,3) dots2(:,3) NaN*ones(numPts,1)]',numPts*3,1);
    else
        zdat = NaN;
    end
catch er
    disp(er.getReport)
    warning('input dots1 and dots2 must be Nx2 or Nx3 matrices.  Unable to link lists'); 
    xdat = NaN;
    ydat = NaN;
    zdat = NaN;
end