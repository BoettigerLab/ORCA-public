function LaunchOdysseyOligoArray(slurm_array_idx)

folder = '/n/home05/boettiger/OdysseyTSTORM/Lib9prep/14-09-05_Probes_30mers_t70Cx72Cs76C/';
load([folder,'oligoArrayCommands.mat'])
disp('recieved');
disp(slurm_array_idx);
disp(['running ',oligoArrayCommands{ slurm_array_idx } ]); %#ok<USENS>  loaded from file above 
system(oligoArrayCommands{ slurm_array_idx } );

