function [newlink,validMove] = BFMmove(B,c)
% [newlink,validMove] = BFMmove(B,c)
% Description:
% compute a bond fluctuation move for random walk polymer speicified by B
% at link specified by c. 
% Outputs: 
% newlink, 3-vector, returns new position for link
% validMove, logical, returns whether move broke the polymer backbone.            

N = size(B,1); 

stepDir = randi([1,3],1);
if stepDir == 1
    newstep = [2*round(2*(rand-.5)),0,0];
elseif stepDir == 2
    newstep = [0 ,2*round(2*(rand-.5)),0];
elseif stepDir == 3
    newstep = [0,0,2*round(2*(rand-.5))];
end  
newlink = B(c,:) + newstep;

% Compute bond lengths. 
if c>1
    bondlength1 = ( newlink(1,1) - B(c-1,1) )^2 + ...
                  ( newlink(1,2) - B(c-1,2) )^2 + ...
                  ( newlink(1,3) - B(c-1,3) )^2;
else
    bondlength1 = 0;
end
if c<N
    bondlength2 = ( newlink(1,1) - B(c+1,1) )^2 + ...
                  ( newlink(1,2) - B(c+1,2) )^2 + ...
                  ( newlink(1,3) - B(c+1,3) )^2;
else
    bondlength2 = 0;
end
validMove = bondlength1 <= 12 && bondlength2 <= 12;