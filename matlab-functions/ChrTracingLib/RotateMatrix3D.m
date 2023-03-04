function [mat3D_t] = RotateMatrix3D(mat3D,tx,tz,varargin)
% rotate around z-axis, then x-axis,

% mat3D = mat3D(3:end,:,:);
defaults = cell(0,3);
defaults(end+1,:) = {'method',{'bilinear','nearest'},'bilinear'};
defaults(end+1,:) = {'cropType',{'loose','crop'},'crop'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

[y1,x1,z1] = size(mat3D);
mat3D_t = mat3D;
mat3D_t = imrotate(mat3D_t,tz,pars.method,pars.cropType);
mat3D_t = imrotate(permute(mat3D_t,[3,1,2]),tx,pars.method,pars.cropType);
mat3D_t = permute(mat3D_t,[2,3,1]); % back straight, just for me to keep
% mat3D_t = imrotate(permute(mat3D_t,[3,2,1]),tx,pars.method,pars.cropType);
% mat3D_t = permute(mat3D_t,[3,2,1]); % back straight, just for me to keep
% mat3D_t = imrotate(permute(mat3D_t,[3,1,2]),ty,method);

% resize
[y2,x2,z2]  = size(mat3D_t);
y3 = floor((y2-y1)/2);
x3 = floor((x2-x1)/2);
z3 = floor((z2-z1)/2);
if y3 < 0
    y3 = 0;
end
if x3 < 0
    x3 = 0;
end
if z3 < 0
    z3 = 0;
end
try
    mat3D_t = mat3D_t(1+y3:end-y3,1+x3:end-x3,1+z3:end-z3);
catch
   disp('debug here'); 
end


    