function a = Asphericity(pts)


% https://boulderschool.yale.edu/sites/default/files/files/rudnick_notes.pdf

% https://en.wikipedia.org/wiki/Gyration_tensor

% % S[m,n] = 
% pts = 10*randn(100,3) + [100*rand(100,1),zeros(100,2)];
% figure(1); clf; plot(pts(:,1),pts(:,2),'k.');

[N,dims] = size(pts);
pts = CenterPolymer(pts);
% S = zeros(dims,dims);
% for i=1:dims
%     for j=1:dims
%         S(i,j) = pts(:,i)'*pts(:,j)/N;
%     end
% end
toRemove = isnan(pts(:,1));
pts(toRemove,:) = [];
N = size(pts,3);
S=pts'*pts/N; % build the gyration tensor

[~,v,~] = svd(S);
v= v./nansum(v(:));
a = v(1,1) - .5*(v(3,3)+v(2,2));