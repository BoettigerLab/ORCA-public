function warp = WarpProps(tform,varargin)
defaults = cell(0,3);
defaults(end+1,:) = {'imsize','freeType',[]};
pars = ParseVariableArguments(varargin,defaults,'WarpProps');

if isempty(pars.imsize)
    % see how it rotates shifts and scales the unit vector at the origin
    [x,y] = transformPointsInverse(tform,[0 1],[0 0]);
else
    % see how it rotates shifts and scales the unit vector at the image center
    [x,y] = transformPointsInverse(tform,[0 1]+pars.imsize/2,[0 0]+pars.imsize/2);
end

    
dx = x(2) - x(1); 
dy = y(2) - y(1); 
warp.rotationAngle = (180/pi) * atan2(dy, dx);
warp.scaleImage = 1 / sqrt(dx^2 + dy^2);
warp.xshift = x(1);
warp.yshift = y(1);
