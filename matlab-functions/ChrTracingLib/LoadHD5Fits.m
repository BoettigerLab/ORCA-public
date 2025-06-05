function data = LoadHD5Fits(binName,varargin)

% Parse default parameters
defaults = cell(0,3);
defaults(end+1,:) = {'method',{'h5read','csv'},'h5read'}; 
defaults(end+1,:) = {'verbose','boolean',true}; 
defaults(end+1,:) = {'model',{'2d','2dfixed','3d','Z'},'2d'}; % 0 uses default in xml file  
defaults(end+1,:) = {'nFrames','integer',0}; % 0 will auto detect, may be slow.  
defaults(end+1,:) = {'frames','integer',0}; % allow user to specify some frames to load 
defaults(end+1,:) = {'parallel','boolean',false}; % allow user to specify some frames to load 
pars = ParseVariableArguments(varargin,defaults,mfilename);

%%

if strcmp(pars.method,'csv')
    data = LoadDaoFits(binName);
else
    
        %%
    try
        a = h5info(binName,'/fr_0');
        dataFields =  {a.Datasets.Name};
    catch
        a = h5info(binName);
        dataFields =  {a.Datasets.Name};
    end
    
    % decide which frames to load
    if pars.frames(1) ~= 0 % if specific frames were passed, we just load those.
        frameNums = pars.frames; % these are the (frame) numbers written in the hdf5 file that we need to load. They should be treated like text.  Except note they are "1" indexed, not "0" indexed    
        nFrames = length(frameNums); 
        selFrames = 1:nFrames;  % these are the indices of from the list above that we are to use in the analysis; 
    elseif pars.frames(1) == 0 && pars.nFrames ~=0
        nFrames = pars.nFrames; 
        selFrames = 1:nFrames; 
        frameNums =  1:nFrames; 
    else % both are 0, we have to go fishing for the frame number. 
        % reading from the bin file is slow cause h5info and h5read are slow
        % It would be faster to get the corresponding dax name and read it's
        % info file.  
        a = h5info(binName);
        frameNames = {a.Groups.Name};
        isFrame = contains(frameNames,'/fr');
        frameNames = frameNames(isFrame);
        frameNums = cellfun(@(x) str2double(x(5:end)),frameNames)+1; % convert to 1-indexed 
        frameNums = sort(frameNums,'ascend'); % don't load in alphabetic order !!
        nFrames = length(frameNums);
        selFrames = 1:nFrames;
    end
    
    
    
    
    data(nFrames).x = [];
    data(nFrames).y = [];
    data(nFrames).height = [];
    data(nFrames).background = [];
    if ~pars.parallel
        
        try
            if pars.verbose
                disp('loading data...')
            end
            for f=selFrames
                try
                    ff = frameNums(f)-1; % frame number is 0 indexed.
                    data(f).frame = ff+1; % matlab will keep things 1 indexed
                    data(f).x = h5read(binName,['/fr_',num2str(ff),'/x']);
                    data(f).y = h5read(binName,['/fr_',num2str(ff),'/y']);
                    data(f).height = h5read(binName,['/fr_',num2str(ff),'/height']);
                    data(f).background = h5read(binName,['/fr_',num2str(ff),'/background']);
                    data(f).sigma = h5read(binName,['/fr_',num2str(ff),'/xsigma']);
                
                if any(contains(dataFields,'z')) % strcmp(pars.model,'Z') || strcmp(pars.model,'3d')
                    data(f).z = h5read(binName,['/fr_',num2str(ff),'/z']);
                end
                % if any(contains(dataFields,'xsigma'))
                %     data(f).sigma = h5read(binName,['/fr_',num2str(ff),'/xsigma']);
                % end
        
                catch er
                    warning(er.message);
                    disp(er.getReport);
                    disp('place debug here');
                end
        
        
                if pars.verbose
                    if rem(f,1000)==0
                        disp([num2str(f/nFrames*100,3),'% complete'])
                    end
                end
        
            end
        catch er
            disp(h5info(binName));
            warning(er.message);
            warning(er.getReport);
        end
    
    else % run parallel
        try
            if pars.verbose
                disp('loading data...')
            end
            parfor f=selFrames
                try
                    ff = frameNums(f)-1; % frame number is 0 indexed.
                    data(f).frame = ff+1; % matlab will keep things 1 indexed
                    data(f).x = h5read(binName,['/fr_',num2str(ff),'/x']);
                    data(f).y = h5read(binName,['/fr_',num2str(ff),'/y']);
                    data(f).height = h5read(binName,['/fr_',num2str(ff),'/height']);
                    data(f).background = h5read(binName,['/fr_',num2str(ff),'/background']);
                    data(f).sigma = h5read(binName,['/fr_',num2str(ff),'/xsigma']);
                
                if any(contains(dataFields,'z')) % strcmp(pars.model,'Z') || strcmp(pars.model,'3d')
                    data(f).z = h5read(binName,['/fr_',num2str(ff),'/z']);
                end
                % if any(contains(dataFields,'xsigma'))
                %     data(f).sigma = h5read(binName,['/fr_',num2str(ff),'/xsigma']);
                % end
        
                catch er
                    warning(er.message);
                    disp(er.getReport);
                    disp('place debug here');
                end
        
        
                if pars.verbose
                    if rem(f,1000)==0
                        disp([num2str(f/nFrames*100,3),'% complete'])
                    end
                end
        
            end
        catch er
            disp(h5info(binName));
            warning(er.message);
            warning(er.getReport);
        end
    end
end