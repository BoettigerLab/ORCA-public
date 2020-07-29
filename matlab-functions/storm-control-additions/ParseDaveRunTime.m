function totalTime = ParseDaveRunTime(daveRunXml,kiloryXml,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};

pars = ParseVariableArguments(varargin,defaults,mfilename);

% kiloryXML = 'E:\Leslie\2017-05-18_WT_14-16hrs_Bxc_10kb\Settings\Kilroy_Configuration_15minHybe.xml';
% xml = 'E:\Leslie\2017-05-18_WT_14-16hrs_Bxc_10kb\Settings\Run_full.xml';

% allow either a xml path or an imported xml-structure to be read. 
if ~isstruct(daveRunXml)
    r = xml2struct(daveRunXml);
else
    r= daveRunXml;
end

% flatten all movies
allMoviePars = {};
allMovieFrames = [];
k = 0;
for i=1:length(r.sequence.branch)
    for j=1:length(r.sequence.branch{i}.branch)
        k = k+1;
        allMoviePars{k} = r.sequence.branch{i}.branch{j}.DASetParameters.parameters.Text; %#ok<*SAGROW>
        allMovieFrames(k) = str2double(r.sequence.branch{i}.branch{j}.DATakeMovie.length.Text); %#ok<*AGROW>
    end
end

% load parameter files
parsFiles = unique(allMoviePars);
folder = fileparts(daveRunXml);

movieTimes = zeros(length(parsFiles),1);
for i=1:length(parsFiles)
    try
        p = xml2struct([folder,filesep,parsFiles{i}]);
    catch
        error(['failed to find parameter file ',folder,filesep,parsFiles{i}]);
    end
    fps = 1/str2double(p.settings.camera1.exposure_time.Text);
    isPar_i = StringFind(allMoviePars,parsFiles{i},'exactly',true);
    movieTimes(i) = sum(allMovieFrames(isPar_i))/fps;
end

if pars.verbose
    for i=1:length(parsFiles)
        disp([parsFiles{i} ' (hh:mm:ss) :']);
        disp( duration(seconds(movieTimes(i)),'Format','hh:mm:ss') );
    end
end

movieTime = sum(movieTimes(:));

allValveCommands = cellfun(@(x) x.Text, r.sequence.DAValveProtocol,'UniformOutput',false)';
fluidTime = ParseKilroyTime(kiloryXml,allValveCommands);

totalTime = movieTime+fluidTime;

disp('total run time: (hh:mm:ss) ');
disp(duration(seconds(totalTime),'Format','hh:mm:ss'));
disp(['Expected completion: ',datestr(now+seconds(totalTime),'mmmm dd, yyyy HH:MM AM')]);




