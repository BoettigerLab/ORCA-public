function ChrTracer_FixChromatic()
% Find Align Hybes (DataType 'A')

% Register align hybes first (this is the same first step)

% Actually, it will be better to apply this to the molecule lists. 


% 
% % Fix chromatic ------------------------------------------------% 
% defaults = cell(0,3);
% defaults(end+1,:) = {};
% CT{handles.id}.parsFixChromatic = ParseVariableArguments([],defaults,'ChrTracer_FixChromatic'); 
% 
% 
% %---------------------------- SAVE FOV OVERLAYS --------------------------%    
% elseif strcmp(currStepName,'Validate Overlay')
%     pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFovOverlay);
%     ChrTracer2_FOVsummaryPlots(CT{handles.id}.fiducialFrames,'parameters',pars);
%     
%     if CT{handles.id}.numDataChns ==1 % no need for chromatic correction
%         CT{handles.id}.currStep = 'Save Aligned Data';
%         disp('Select Save Aligned Data if the alignment looks good');
%         disp('Otherwise select Back Step, change Parameters, and Reload Data');
%     else
%          CT{handles.id}.currStep = 'Fix Chromatic';
%     end
% %-------------------------- FIX CHROMATIC ------------------------------% 
% elseif strcmp(currStepName,'Fix Chromatic')
%     pars = CT{handles.id}.parsChromatic;
%     ChrTracer2_FixChromatic();