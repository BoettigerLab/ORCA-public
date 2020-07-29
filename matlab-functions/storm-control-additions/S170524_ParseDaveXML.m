% function ParseXML

xml = 'E:\Leslie\2017-05-18_WT_14-16hrs_Bxc_10kb\Settings\dave_Recipe_Bxc_Bleach_488.xml';
xml = 'E:\Leslie\2017-05-18_WT_14-16hrs_Bxc_10kb\Settings\dave_Recipe_Bxc.xml';


d = xml2struct(xml)

% this is the position loop
d.recipe.loop_variable.Attributes
d.recipe.loop_variable.file_path

d.recipe.command_sequence.valve_protocol{1}
allValveCommands = cellfun(@(x) x.Text, d.recipe.command_sequence.valve_protocol,'UniformOutput',false)';

%% Parse the run file
xml = 'E:\Leslie\2017-05-18_WT_14-16hrs_Bxc_10kb\Settings\Run_full.xml'

% 
r = xml2struct(xml)
allValveCommands = cellfun(@(x) x.Text, r.sequence.DAValveProtocol,'UniformOutput',false)';

% flatten all movies
allMoviePars = {};
allMovieFrames = [];
k = 0;
for i=1:length(r.sequence.branch)
    for j=1:length(r.sequence.branch{i}.branch)
        k = k+1;
        allMoviePars{k} = r.sequence.branch{i}.branch{j}.DASetParameters.parameters.Text;
        allMovieFrames(k) = str2double(r.sequence.branch{1}.branch{1}.DATakeMovie.length.Text);
    end
end

scan_fps = 5;
isBleach = StringFind(allMoviePars,'Bleach');
isScan = StringFind(allMoviePars,'scan');
bleachTime = sum(allMovieFrames(isBleach));
scanTime = sum(allMovieFrames(isScan))/scan_fps;

movieTime = (bleachTime + scanTime);

% parse kilroy
kilroyXML = xml2struct('E:\Leslie\2017-05-18_WT_14-16hrs_Bxc_10kb\Settings\Kilroy_Configuration_15minHybe.xml')
kiloryProtocolNames = cellfun(@(x) x.Attributes.name,kilroyXML.kilroy_configuration.kilroy_protocols.protocol,'UniformOutput',false);

% usedKPs = StringFind(kiloryProtocolNames,allValveCommands)
numKPs = length(kiloryProtocolNames);
usedKPs = StringFind(allValveCommands,kiloryProtocolNames,'exactly',true);
KPtime = zeros(numKPs,1);
for i=1:numKPs
    % flatten valves
    valve = kilroyXML.kilroy_configuration.kilroy_protocols.protocol{i}.valve;
    if iscell(valve)
        for j=1:length(valve)
            KPtime(i) = KPtime(i)+ str2double(valve{j}.Attributes.duration);
        end
    else
        KPtime(i) = KPtime(i)+ str2double(valve.Attributes.duration);
    end
    pump  = kilroyXML.kilroy_configuration.kilroy_protocols.protocol{i}.pump;
    if iscell(pump)
        for j=1:length(pump)
            KPtime(i) = KPtime(i)+ str2double(pump{j}.Attributes.duration);
        end
    else
        KPtime(i) = KPtime(i)+ str2double(pump.Attributes.duration);
    end
    KPtime(i) = KPtime(i)*length(usedKPs{i});
end
    

fluidTime = sum(KPtime);
duration(seconds(fluidTime),'Format','hh:mm:ss')
duration(seconds(movieTime+fluidTime),'Format','hh:mm:ss')
disp(['Expected completion: ',datestr(now+seconds(movieTime),'mmmm dd, yyyy HH:MM AM')]);






