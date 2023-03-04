function PlotMQ(data,varargin)
% data is O-replicates by N-data points. PlotMQ takes the median and
% quantile range over O and renders these as errorbars. 

defaults = cell(0,3);
defaults(end+1,:) = {'quantiles','array',[.025,.975]}; % default is 95% confidence interval
defaults(end+1,:) = {'color','colormap','k'}; % 
defaults(end+1,:) = {'plotOptions','array',{}};
pars = ParseVariableArguments(varargin,defaults,mfilename);

[o,n] = size(data);
if o>1
    med  = Row(quantile(data,.5,2));
    lo =Row(quantile(data,pars.quantiles(1),2));
    hi = Row(quantile(data,pars.quantiles(2),2));
else
    med = data;
end

plot(med,'.','color',pars.color,'MarkerSize',10); 
if o>1
    hold on; x=1:length(lo); plot([x;x],[lo;hi],'-','color',pars.color);
end
