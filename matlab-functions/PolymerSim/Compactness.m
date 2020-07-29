function c = Compactness(xyz)
warning('function under development');

% xyz =  double([chrData(n).flists{d}.xc*npp,chrData(n).flists{d}.yc*npp,chrData(n).flists{d}.zc]);
% figure(3); clf; plot3(xyz1(:,1),xyz1(:,2),xyz1(:,3),'.','color','r'); pause(.1);

 
 while length(xyz) > 1E4;
     xyz = xyz(1:2:end,:);
 end
 [nPts,dim] = size(xyz);

pY = pdist(xyz);
pD = squareform(pY);
rg = RadiusOfGyration(xyz);
if dim == 3
    c = median(sum(pD < rg))/nPts;% / rg.^3 ;
elseif dim == 2
    c = median(sum(pD < rg))/nPts;% / rg.^2 ;
else
    error('invalid dimension');
end
