function fluidTime = ParseKilroyTime(xml,allValveCommands,varargin)
% Inputs: xml- a kilroy xml configuration file 
%        allValveCommands - a cell array containing the names of all valve
%        commands used, repeated if used multiple times.
% Outputs: fluidTime - total time in seconds required to run all kilroy
%        commands.

kilroyXML = xml2struct(xml);
kiloryProtocolNames = cellfun(@(x) x.Attributes.name,kilroyXML.kilroy_configuration.kilroy_protocols.protocol,'UniformOutput',false);
% kilroyXML.

numKPs = length(kiloryProtocolNames);
usedKPs = StringFind(allValveCommands,kiloryProtocolNames,'exactly',true,'cellOutput',true);
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
display('Fluid time:');
display( duration(seconds(fluidTime),'Format','hh:mm:ss') )