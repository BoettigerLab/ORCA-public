function s1 = StackLocalDistances(currMaps,a,b)
% a short function to combine all the pairs of points around a specific
% pair of interest, to increase the amount of data available (using local
% distances as a proxy). 

t = 1; % tile size - could become optional variable in future

s1 = {}; % output. will convert to vector
k = 0;
for p=-t:t
    for q=-t:t
        k=k+1;
        s1{k} = squeeze(currMaps(a+p,b+q,:));
    end
end
s1 = cat(1,s1{:});