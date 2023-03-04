function pmap = PValMap(map1,map2)   

[nB,~,nC] = size(map1);
pmap_T = zeros(nB,nB);
for a=1:nB
    for b=a+1:nB
        try
            pmap_T(a,b) = ranksum( squeeze(map1(a,b,:)),squeeze(map2(a,b,:)) );
        catch
        end
    end
end
pmap= pmap_T + pmap_T';
