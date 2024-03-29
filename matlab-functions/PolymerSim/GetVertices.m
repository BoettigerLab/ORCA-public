function occupiedVertices = GetVertices(B,varargin)

N = size(B,1);
if nargin == 1
    M = 2*N+2;
else
    M = varargin{1};
end

      
X =   [2*B(:,1)-1,...
      2*B(:,1)-1,...
      2*B(:,1)+1,...
      2*B(:,1)+1,...
      2*B(:,1)-1,...
      2*B(:,1)-1,...
      2*B(:,1)+1,...
      2*B(:,1)+1];
Y =  [ 2*B(:,2)+1,...
      2*B(:,2)+1,...
      2*B(:,2)+1,...
      2*B(:,2)+1,...
      2*B(:,2)-1,...
      2*B(:,2)-1,...
      2*B(:,2)-1,...
      2*B(:,2)-1];
Z =   [2*B(:,3)+1,...
      2*B(:,3)-1,...
      2*B(:,3)+1,...
      2*B(:,3)-1,...
      2*B(:,3)+1,...
      2*B(:,3)-1,...
      2*B(:,3)+1,...
      2*B(:,3)-1];
              
occupiedVertices = sub2indFast(4*M,X,Y,Z);


      
     
