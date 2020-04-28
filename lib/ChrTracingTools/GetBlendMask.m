function blendMask = GetBlendMask(existingData,varargin)

[nRows,nCols] = size(existingData);

defaults = cell(0,3);
defaults(end+1,:) = {'blendRadius','integer',round(.1*min([nRows,nCols]))};
defaults(end+1,:) = {'troubleshoot','boolean',false};
defaults(end+1,:) = {'symmetric','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% % 
% im = zeros(100,100);
% im(40:end,1:3) = 1;
% im(20:60,1:30) = 1;
% im(12:20,80:end) = 1;
% existingData = im;
% figure(8); clf; imagesc(existingData);
%%
try 

%%
imE = existingData > 0;
% figure(7); clf; imagesc(imE); colorbar;
imEdge = edge(imE,'sobel');
imEdge = imEdge>0;
% figure(8); clf; imagesc(imEdge); colorbar;


[nRows,nCols] = size(imE);
x1 = find(sum(imEdge,1) > 5,1,'first');
x2 = find(sum(imEdge,1) > 5,1,'last');
y1 = find(sum(imEdge,2) > 5,1,'first');
y2 = find(sum(imEdge,2) > 5,1,'last');

if pars.symmetric
    rr = min([pars.blendRadius,x1,nCols-x2,y1,nRows-y2]);
    H = [0,linspace(0,1-1/rr,rr),1,fliplr(linspace(0,1-1/rr,rr)),0];
    Hx = repmat(H,rr*2+3,1);
    Hy = repmat(H',1,rr*2+3);
    H =  Hx.*Hy;
else
    rx = min([pars.blendRadius,x1,nCols-x2]);
    ry = min([pars.blendRadius,y1,nRows-y2]);
    Hx = [0,linspace(0,1-1/rx,rx),1,fliplr(linspace(0,1-1/rx,rx)),0];
    Hx = repmat(Hx,ry*2+3,1);
    Hy = [0,linspace(0,1-1/ry,ry),1,fliplr(linspace(0,1-1/ry,ry)),0]';
    Hy = repmat(Hy,1,rx*2+3);
    H =  Hx.*Hy;
end
% figure(6); clf; imagesc(H); colorbar;

Y = filter2(H,imEdge);
norm = mode(Y(imEdge));
if isnan(norm) || norm==0
    blendMask = zeros(nRows,nCols);
else
    Y = Y./norm; Y(Y>1) = 1;
    % figure(8); clf; imagesc(Y); colorbar;

    blendMask = imE.*(1-Y);
end
% 
% figure(8); clf; imagesc(1-Y); colorbar;
% figure(9); clf; imagesc(blendMask); 

catch er
    disp(er.getReport);
    disp('debug now');
    error('problem here');
end

if pars.troubleshoot 
    figure(9); clf; subplot(1,2,1); 
    imagesc(blendMask); colorbar;
    subplot(1,2,2); imagesc(1-blendMask); colorbar;
end
