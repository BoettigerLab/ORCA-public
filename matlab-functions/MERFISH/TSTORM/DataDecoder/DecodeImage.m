function [uniqueMsg,hammingDis,freqUniqueMsg] = DecodeImage(wordsDetected,codebook)

%%
showPlots = false;

[uniqueMsg,~,ind] = unique(wordsDetected,'rows'); 
freqUniqueMsg = histc(ind,1:numel(uniqueMsg(:,1)));
[sortedUniqueMsg,sIdx] = sort(freqUniqueMsg);
ascUniqueMsg = flipud(sortedUniqueMsg);
ascIdx = flipud(sIdx);

Ngenes = length(codebook); 
codes = zeros(Ngenes,length(str2num(codebook(1).Header))); %#ok<*ST2NM>
for n = 1:Ngenes
    codes(n,:) = str2num(codebook(n).Header);
end 
hammingDis = zeros(length(uniqueMsg),Ngenes);
for i=1:length(uniqueMsg)
    for n=1:Ngenes
        hammingDis(i,n) = sum(codes(n,:) ~= uniqueMsg(i,:));
    end    
end

if showPlots
    fHand = figure(7);clf;
    aHand = axes('parent', fHand);
    hold(aHand, 'on')
    colors = lines(4);
    for i = 1:numel(ascUniqueMsg)
        if any(hammingDis(ascIdx(i),:) == 0)
            bar(i, ascUniqueMsg(i), 'parent', aHand, 'facecolor', colors(3,:));
        elseif any(hammingDis(ascIdx(i),:) == 1)     
            bar(i, ascUniqueMsg(i), 'parent', aHand, 'facecolor', colors(2,:));
        elseif any(hammingDis(ascIdx(i),:) == 2)     
            bar(i, ascUniqueMsg(i), 'parent', aHand, 'facecolor', colors(1,:));
        elseif any(hammingDis(ascIdx(i),:) == 3)     
            bar(i, ascUniqueMsg(i), 'parent', aHand, 'facecolor', colors(1,:));        
        else
            bar(i, ascUniqueMsg(i), 'parent', aHand, 'facecolor', colors(4,:));
        end    
    end
    x = 1:length(uniqueMsg);
    for i=1:length(uniqueMsg)
        uLabel = text(x(i),freqUniqueMsg(ascIdx(i))+10,num2str(uniqueMsg(ascIdx(i),:)));
        set(uLabel,'rotation',90)
    end    
    set(gca,'xticklabel','')
end