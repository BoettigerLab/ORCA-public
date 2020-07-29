
% Run splitQVdax on single channel and mutlichannel data
minSize = 5E2; % Mb

folders ={    
    'L:\140310_L4E4'
} 
gain  = 1; 
        
% Run bead averaging on all data
for i = 1:length(folders); 
    beaddax = dir([folders{i},'\','561quad_*.dax']);
    idxBeads = [beaddax.bytes]/1E6 > minSize;
    beaddaxNames = strcat([folders{i},'\'],{beaddax(idxBeads).name})';

    for n=1:length(beaddaxNames)
        if isempty(strfind(beaddaxNames{n},'ds'))
            disp(['Averaging movie ',beaddaxNames{n}]);
            try 
                TimeAverageMovie(beaddaxNames{n},'gain',gain); 
            catch er
                disp(er.getReport);
                disp(['TimeAverageMovie Failed on ',beaddaxNames{n}]);
                beaddaxNames{n} = [];
            end
        else
            disp(['skipping movie ',beaddaxNames{n}]);
        end
    end
    
%     pause(1);
%     disp('deleting old bead-files...')
%     for n=1:length(beaddaxNames)
%         try
%          if isempty(strfind(beaddaxNames{n},'ds'))
%              downsampledName = regexprep(beaddaxNames{n},'\.dax','_ds60\.dax');
%              if exist(downsampledName,'file') > 0
%                 delCmd = ['del ',beaddaxNames{n}];
%                 disp(delCmd);
%                 system(delCmd); 
%              end
%          else
%             disp(['not deleting ds movie ',beaddaxNames{n}]);
%          end
%         catch
%              disp(['not deleting ds movie ',beaddaxNames{n}]);
%         end
%    end
    
end

