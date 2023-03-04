function [xyz_t, axisVecs] = RotateVector3D(xyz,tx,ty,tz,varargin)
% xyz_t = RotateVector3D(xyz,tx,ty,tz,varargin)
% rotates the vector xyz in 3D space. 
% tx = degrees around the x-axis
% ty = degrees around the y-axis
% tz = degrees around the z-axis
% optional: spotCenter - rotate around this point 

if length(varargin) > 0
    spotCenter = varargin{1};
else
    spotCenter = zeros(3,1);
end
if numel(xyz)==3
    xyz = Column(xyz);
end

% convert degrees to radians
tx = tx*2*pi/360;
ty = ty*2*pi/360;
tz = tz*2*pi/360;

% rotation matrices
rx =  @(tx) ...
      [1       0       0;
      0     cos(tx) -sin(tx);
      0     sin(tx) cos(tx)];
ry = @(ty) ...
    [cos(ty)   0   sin(ty);
        0       1       0;
      -sin(ty)  0   cos(ty)];
rz = @(tz) ...
    [cos(tz) -sin(tz)  0;
      sin(tz) cos(tz)   0;
        0       0       1];
    
% shift center    
xyz = xyz - spotCenter;

% rotate
r = rz(tz)*ry(ty)*rx(tx);
xyz_t = r*xyz;

% shift center back
xyz_t = xyz_t + spotCenter;

% Axis
axisVecs = cat(2, r*[1;0;0],  r*[0;1;0],  r*[0;0;1]);

    