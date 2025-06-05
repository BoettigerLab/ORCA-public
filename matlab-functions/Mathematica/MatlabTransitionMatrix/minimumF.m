function minexpr=minimumF(q,exps,indxcodes)
% Note this requires the reduced function for the reduced transition matrix WR
M=WR(q);
MR=M(indxcodes,:);
while ~det(MR)
    MR=MR+rand(size(MR))/10^8;
end
A=inv(MR)*exps(indxcodes)';
zerovec=M*A-exps';
minexpr=zerovec'*zerovec;
end