function [correctClusters,spotMatrix] = PlotGeneLocations(mRNAcents,Mcents,uniqueMsg,hammingDis,codebook,xf,yf,libGenes,maxDtoCentroid,showbkd)


numGenes = length(codebook);
numHybes = size(Mcents,2);
clrmap = hsv(numHybes); 


clear spotMatrix;
cntsPerSpot = cell(numGenes,1);
spotMatrix = cell(numGenes,1); % figure(40); clf;
for h = 1:numHybes
    plot(0,0,'.','color',clrmap(h,:),'MarkerSize',10); hold on;
end
colormap(clrmap);
legend(num2str( (1:numHybes)' )); hold on; 
 

correctClusters = cell(numGenes,1);
for i=1:numGenes
    correctClusters{i} = ismember(Mcents,...
            uniqueMsg((hammingDis(:,i)==1|hammingDis(:,i)==0),:),'rows');
end


for h = 1:numHybes
    if showbkd
        plot(xf{h},yf{h},'.','color',clrmap(h,:),'MarkerSize',1); hold on;
    end
    for n=1:numGenes
        [idx,dist] = knnsearch( mRNAcents(correctClusters{n},:) , [xf{h},yf{h}] );
        inRange = dist<maxDtoCentroid;  
        cntsPerSpot{n}{h} = hist(idx(inRange),1:size(mRNAcents(correctClusters{n},:),1) );
        spotMatrix{n} = cat(1,spotMatrix{n},cntsPerSpot{n}{h});
        
        plot(xf{h}(inRange),yf{h}(inRange),'.','color',clrmap(h,:)); hold on;
        text(mRNAcents(correctClusters{n},1),mRNAcents(correctClusters{n},2),libGenes{n});
    end
end