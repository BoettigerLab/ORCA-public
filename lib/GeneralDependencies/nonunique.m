function [B,ia,values,n] = nonunique(A)
%--------------------------------------------------------------------------
% B = nonunique(A) 
%     Returns the elements of A that are non-unique.  A is NxM which has N 
%     potentionally non-unique entries of M values each.  
% [B,ia] = nonunique(A)
%     Also returns the indices ia such  that A(ia)=B
%
% With code from: Joe.
% (http://www.mathworks.com/matlabcentral/newsreader/view_thread/30051)
%--------------------------------------------------------------------------

% Sample Data
%-------------------------------------
% A = rand(5,2);
% A = [A; A(1:3,:); A(1:2,:); A(2,:)];

[N,M] = size(A);
if N<1
    A = A';
    [N,M] = size(A);
end

% sort
A1=sortrows([A,[1:N]']);
% then find duplicates position:
d = abs(diff(A1(:,1:M)));
J=find(sum(d,2)==0);
% Pick the numbers:
ia1=unique([J(:)',J(:)'+1]); % the unique is needed when there are more than two duplicates.
ia = A1(ia1,3);
 B = A(ia,:); 
% B=A1(ia1,1:M); % troubleshooting check
    


    
