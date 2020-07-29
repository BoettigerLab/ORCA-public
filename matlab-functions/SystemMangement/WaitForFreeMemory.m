function WaitForFreeMemory(varargin)
defaults = cell(0,3);
defaults(end+1,:) = {'minFreePhysicalMemory',nonnegative',2E6} % in Bytes
defaults(end+1,:) = {'refreshTime',nonnegative',1} % in Seconds
defaults(end+1,:) = {'verbose',boolean',true} % 


[~,wmicOut] = system('wmic os');
sysData = strsplit(wmicOut)
freePhysicalMemory = str2double(sysData(89));
while freePhysicalMemory < parameters.minFreePhysicalMemory
	if parameters.verbose
		disp('waiting for free RAM...')
		end
		gotpaused = true;
		pause(parameters.refreshTime);
		[~,wmicOut] = system('wmic os');
		sysData = strsplit(wmicOut)
		freePhysicalMemory = str2double(sysData(89));
	end
	if parameters.verbose && gotpaused
		disp('now running...');
	end
else
	warning('WaitForFreeMemory only works for windows');
end