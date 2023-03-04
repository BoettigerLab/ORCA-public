function Rg = RadiusOfGyration(B)
% compute the radius of gyration for the 3D-polymer specified by B
% Rg = RadiusOfGyration(B)
% B is a N x 3 matix


B(isnan(B(:,1)),:) = [];

N = size(B,1);
dim = size(B,2);


 rm = mean(B,1 );
 if dim == 3
    Rd = (B(:,1)-rm(1)).^2 + (B(:,2)-rm(2)).^2 + (B(:,3)-rm(3)).^2;
    Rg = sqrt( sum( Rd)/N);
 elseif dim == 2
    Rd = (B(:,1)-rm(1)).^2 + (B(:,2)-rm(2)).^2;
    Rg = sqrt( sum( Rd)/N);
 else
     disp('warning: empty polymer, Rg=NaN'); 
     Rg = NaN; 
 end