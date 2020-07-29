function polymer = RandomPolymer(N)
% returns a non-self-avoiding random polymer
polymer = [];
while isempty(polymer)
    polymer = rand(3*N,3)*2 - 1; % -1;
    dists = sqrt( sum(polymer.^2,2) );
    polymer = polymer(dists <= 1,:);
    dists = sqrt( sum(polymer.^2,2) );
    polymer = polymer./repmat(dists,1,3);
    polymer = cumsum(polymer);
    try
        polymer = polymer(1:N,:);
    catch
        polymer = [];
    end
end