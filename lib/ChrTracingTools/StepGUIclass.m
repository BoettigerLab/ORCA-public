classdef StepGUIclass < handle

% MosaicBuilder specific dynamic properties
    properties % editable parameters
        gui_h;
% %  The following properties must be defined by the child function:
%         currStep = 'Step 1';
%         saveFolder = '';
%        allStepNames = {'Step 1',...
%                        'Step 2',...
%                        'Step 3'}
%        stepDirections = {'Step 1: load stuff';
%                         'Step 2: do more stuff';
%                         'Step 3: finish stuff'}
%        allFxnHandles = {};
%        stepParameters = {};
     end

    
    methods            
        function self = StepGUIclass
            % make GUI handle to the fig.m file and store it locally
            self.gui_h = guihandles(MultistepAnalysisGUI);
        
            % set the callback function for the "Run Step" button.
            set(self.gui_h.ButtonRun,'callback',@(src,event) ButtonRun_callback(self,src,event));
            
            % set the callback function for the "Next Step" button.
            set(self.gui_h.ButtonNext,'callback',@(src,event) ButtonNext_callback(self,src,event));
            
            % set the callback function for the "Back Step" button.
            set(self.gui_h.ButtonBack,'callback',@(src,event) ButtonBack_callback(self,src,event));
            
            % set the callback function for the "Parameters" button.
            set(self.gui_h.ButtonPars,'callback',@(src,event) ButtonPars_callback(self,src,event));
            
            
            %sets the figure close function. This lets the class know that
            %the figure wants to close and thus the class should cleanup in
            %memory as well
            set(self.gui_h.figure1,  'closerequestfcn', @(src,event) Close_fcn(self,src,event));
            
            % setup initial function
            set(self.gui_h.TextDir,'String',self.stepDirections{1});
        end
        
    end
    % ==== These methods could be part of a MultistepAnalysisClass =======
    methods 
        % ---- Executes when "Run Step" button is clicked
        function ButtonRun_callback(self,src,event) %#ok<*INUSD>
            % calls the current step behavior using the current step
            % parameters
            %    see protected functions below
            for s=1:length(self.allStepNames)
                if strcmp(self.currStep,self.allStepNames(s))
                    pars = self.stepParameters{s};
                    self.allFxnHandles{s}(self,pars);
                end
            end
        end
        
        % ---- Executes when "Parameters" button is clicked
        function ButtonPars_callback(self,src,event)
            % Opens a GUI to edit the current step parameters
            for s=1:length(self.allStepNames)
                if strcmp(self.currStep,self.allStepNames(s))
                    self.stepParameters{s} = SimpleParameterGUI(self.stepParameters{s});
                end
            end
        end
        
        % ---- Executes when "Next Step" button is clicked
        function ButtonNext_callback(self,src,event)
            
            step = find(strcmp(self.allStepNames, self.currStep));
            numSteps = length(self.allStepNames);
            if step < numSteps
                 step = step + 1;
                 self.currStep = self.allStepNames{step};
                 set(self.gui_h.TextDir,'String',self.stepDirections{step});
            else
                warning('Already reached last step'); 
            end
        end
         
        % ---- Executes when "Back Step" button is clicked
        function ButtonBack_callback(self,src,event)
            step = find(strcmp(self.allStepNames, self.currStep));
            if step > 1
                 step = step - 1;
                 self.currStep = self.allStepNames{step};
                 set(self.gui_h.TextDir,'String',self.stepDirections{step});
            else
                warning('Currently on step 1, can not go back.'); 
            end
        end

        
        % Close Figure
        %this is the closerequestfcn of the figure. All it does here is
        %call the class delete function (presented below)
        function self = Close_fcn(self, src, event)
            % disp('exiting StepGUI');
            delete(self);
        end
        
        %  This (intentionally?) overloads the delete function
        %class deconstructor - handles the cleaning up of the class &
        %figure. Either the class or the figure can initiate the closing
        %condition, this function makes sure both are cleaned up
        function delete(self)
            %remove the closerequestfcn from the figure, this prevents an
            %infitie loop with the following delete command
            set(self.gui_h.figure1,  'closerequestfcn', '');
            %delete the figure
            delete(self.gui_h.figure1);
            %clear out the pointer to the figure - prevents memory leaks
            self.gui_h = [];
        end
        %function - Close_fcn
        %

        

    end
end