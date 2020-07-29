function sld = CounterSlider(im1,varargin)
disp('view output.UserData to access count info');
defaults = cell(0,3);
defaults(end+1,:) = {'color','colormap','y'};
pars  = ParseVariableArguments(varargin,defaults,mfilename);
           
% figure(1); clf; imagesc(imS);

fs = uifigure('Position',[500,500,400,50]); 
pnl = uipanel(fs,'Position',[10,10,400,40]);
sld = uislider(pnl,'Position',[60 30 300 3]);
txt = uilabel(pnl,'Position',[10,10,30,30]);
if ~isempty(pars.color)
   txt.Text = string(pars.color); 
end
sld.Value = 50;
hold on; 
sldData.image = im1;
sldData.pars = pars;
sld.UserData.sH = [];
sld.ValueChangedFcn = @(sld,event) UpdateSpotCount(sld,sldData);
% sld.DeleteFcn = @(sld,event) ExportValues;  % disp('close when done'); 

function [sldData] = UpdateSpotCount(sld,sldData)
% ----- apply simple absolute threshold

rnaIm = sldData.image;
imMax = max(rnaIm(:));
theta = sld.Value/100*imMax;
bw = imregionalmax(rnaIm); 
bw2 = bw;
bw2(rnaIm<theta)  = 0;
[h,w] = size(rnaIm(:,:,1));
[y,x] = ind2sub([h,w],find(bw2));
brightness = rnaIm(bw2); 
if length(x) <= 1
    disp('no spots found!');
else
    sld.UserData.spotCounts = table(x,y,brightness);
    % delete(findobj(gca,'color',sldData.pars.color)); 
    if ~isempty(sld.UserData.sH)
       delete(sld.UserData.sH); 
    end
    hold on;
    sld.UserData.sH = plot(x,y,'o','color',sldData.pars.color);
    hold off;
    % assignin('base','sldData',sldData); % overwrites
end