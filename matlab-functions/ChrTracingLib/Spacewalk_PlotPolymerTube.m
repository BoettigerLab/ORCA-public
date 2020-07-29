function Spacewalk_PlotPolymerTube(poly,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'figHandle','freeType',105};
defaults(end+1,:) = {'fillmissing','boolean',true};
defaults(end+1,:) = {'tubeRadius','positive',20};
defaults(end+1,:) = {'sphereRadius','positive',30};
defaults(end+1,:) = {'showSpheres','boolean',true};
defaults(end+1,:) = {'showTube','boolean',true};
defaults(end+1,:) = {'w','positive',1200}; % box size
defaults(end+1,:) = {'maxJump','positive',inf}; % box size
defaults(end+1,:) = {'numColors','positive',[]}; % box size
defaults(end+1,:) = {'cent','array',[]}; % box size
defaults(end+1,:) = {'colormap','colormap','hsvCut'}; 
defaults(end+1,:) = {'showLabels','boolean',false}; 
defaults(end+1,:) = {'labels','cell',{}}; 
defaults(end+1,:) = {'alpha','fraction',.5}; 
defaults(end+1,:) = {'lightOn','boolean',true}; 
defaults(end+1,:) = {'fit','boolean',true};
defaults(end+1,:) = {'fontSize','positive',20};
defaults(end+1,:) = {'autoAxisLabels','boolean',false};
defaults(end+1,:) = {'method',{'spline','skip'},'skip'};
defaults(end+1,:) = {'center','boolean',true};
defaults(end+1,:) = {'number','boolean',false}; % obsolete
pars = ParseVariableArguments(varargin,defaults,mfilename);

figure(pars.figHandle); clf;
PlotPolymerTube(poly,'parameters',pars);