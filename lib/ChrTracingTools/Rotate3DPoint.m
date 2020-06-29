function [x_rot,y_rot,z_rot]= Rotate3DPoint(xyz,theta,cm)
% rotate 3D data (an Nx3 matrix) a total angle theta (tx,ty,tz)
% around the point cm.  
% angles must be in radians
% 

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);
% tx = theta(:,1);
% ty = theta(:,2);
% tz = theta(:,3); 
tx = theta(1);
ty = theta(2);
tz = theta(3); 

cmx = cm(1); %  mean(x);
cmy = cm(2); % mean(y);
cmz = cm(3); % mean(z);

x = x - cmx;
y = y - cmy;
z = z - cmz;

% enforce row vector
sx=size(x);
if sx(1)>sx(2)
   x=x';
end
sy=size(y);
if sy(1)>sy(2)
   y=y';
end
sz=size(z);
if sz(1)>sz(2)
   z=z';
end
xyz=[x;y;z];
% tx=0.3;
% ty=0.2;
% tz=0.4;
Rx=[1,0,0;...
    0,cos(tx),sin(tx);...
    0,-sin(tx),cos(tx)];
Ry=[cos(ty),0,sin(ty);...
    0,1,0;...
    -sin(ty),0,cos(ty)];
Rz=[cos(tz),sin(tz),0;...
    -sin(tz),cos(tz),0;...
    0,0,1];
xyz_rot=Rz*Ry*Rx*xyz;
x_rot=xyz_rot(1,:);
y_rot=xyz_rot(2,:);
z_rot=xyz_rot(3,:);
if sx(1)>sx(2)
   x_rot=x_rot';
end
if sy(1)>sy(2)
   y_rot=y_rot';
end
if sz(1)>sz(2)
   z_rot=z_rot';
end


x_rot = x_rot  + cmx;
y_rot = y_rot + cmy;
z_rot = z_rot + cmz;

if nargout == 1
    x_rot = [x_rot,y_rot,z_rot];
end
end
