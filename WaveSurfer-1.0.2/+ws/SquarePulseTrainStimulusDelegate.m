classdef SquarePulseTrainStimulusDelegate < ws.StimulusDelegate
    
    properties (Constant)
        TypeString='SquarePulseTrain'
        AdditionalParameterNames={'Delay' 'Duration' 'Amplitude' 'DCOffset' 'Period' 'PulseDuration'};
        AdditionalParameterDisplayNames={'Delay' 'Duration' 'Amplitude' 'DC Offset' 'Period' 'Pulse Duration'};        
        AdditionalParameterDisplayUnitses={'s' 's' '' '' 's' 's'};        
    end
    
    properties (Dependent=true)
        Delay  % this and all below are string sweep expressions, in seconds
        Duration % s
            % this is the duration of the "support" of the stimulus.  I.e.
            % the stim is zero for the delay period, generally nonzero for
            % the Duration, then zero for some time after, with the total
            % duration being specified elsewhere.
        Amplitude
        DCOffset
        Period 
        PulseDuration
    end
    
    properties (Access=protected)
        Delay_ = '0.25'  % sec
        Duration_ = '0.5'  % sec
        Amplitude_ = '5'
        DCOffset_ = '0'
        Period_ = '0.1'  % s
        PulseDuration_ = '0.05'  % s
    end
    
    properties (Dependent = true, Transient=true)
        EndTime  % Delay + Duration
    end
    
    methods
        function val = get.EndTime(self)
            val = ws.Stimulus.evaluateSweepExpression(self.Delay,1) + ws.Stimulus.evaluateSweepExpression(self.Duration,1) ;
        end
    end        
    
    methods
        function self = SquarePulseTrainStimulusDelegate()
            %self=self@ws.StimulusDelegate();
        end
        
        function set.Delay(self, value)
            test = ws.Stimulus.evaluateSweepExpression(value,1) ;
            if ~isempty(test) && isnumeric(test) && isscalar(test) && isfinite(test) && isreal(test) && test>=0 ,
                % if we get here without error, safe to set
                self.Delay_ = value ;
            end                    
        end  % function
        
        function set.Duration(self, value)
            test = ws.Stimulus.evaluateSweepExpression(value,1) ;
            if ~isempty(test) && isnumeric(test) && isscalar(test) && isreal(test) && isfinite(test) && test>=0 ,
                % if we get here without error, safe to set
                self.Duration_ = value;
            end                    
        end  % function
        
        function set.Amplitude(self, value)
            test = ws.Stimulus.evaluateSweepExpression(value,1) ;
            if ~isempty(test) && (isnumeric(test) || islogical(test)) && isscalar(test) && isfinite(test) && isreal(test) ,
                % if we get here without error, safe to set
                self.Amplitude_ = value;
            end                
        end
        
        function set.DCOffset(self, value)
            test = ws.Stimulus.evaluateSweepExpression(value,1) ;
            if ~isempty(test) && isnumeric(test) && isscalar(test) && isfinite(test) && isreal(test) ,
                % if we get here without error, safe to set
                self.DCOffset_ = value;
            end                
        end

        function out = get.Delay(self)
            out = self.Delay_;
        end   % function

        function out = get.Duration(self)
            out = self.Duration_;
        end   % function

        function out = get.Amplitude(self)
            out = self.Amplitude_;
        end   % function

        function out = get.DCOffset(self)
            out = self.DCOffset_;
        end   % function        
        
        function set.Period(self, value)
            test = ws.Stimulus.evaluateSweepExpression(value,1) ;
            if ~isempty(test) && isnumeric(test) && isscalar(test) && isfinite(test) && isreal(test) && test>0 ,
                % if we get here without error, safe to set
                self.Period_ = value;
            end                    
        end  % function
        
        function set.PulseDuration(self, value)
            test = ws.Stimulus.evaluateSweepExpression(value,1) ;
            if ~isempty(test) && isnumeric(test) && isscalar(test) && isfinite(test) && isreal(test) && test>0 ,
                % if we get here without error, safe to set
                self.PulseDuration_ = value;
            end                    
        end  % function

        function out = get.Period(self)
            out = self.Period_;
        end
        
        function out = get.PulseDuration(self)
            out = self.PulseDuration_;
        end
    end
    
    methods
        function data = calculateSignal(self, t, sweepIndexWithinSet)
            % Process args
            if ~exist('sweepIndexWithinSet','var') || isempty(sweepIndexWithinSet) ,
                sweepIndexWithinSet=1;
            end
                        
            % Compute the delay from the expression for it
            delay = ws.Stimulus.evaluateSweepExpression(self.Delay,sweepIndexWithinSet) ;
            % Screen for illegal values
            if isempty(delay) || ~(isnumeric(delay)||islogical(delay)) || ~isscalar(delay) || ~isreal(delay) || ~isfinite(delay) || delay<0 ,
                data=zeros(size(t));
                return
            end
            
            % Shift the timeline to account for the delay
            tShiftedByDelay=t-delay;
            
            % Call a likely-overloaded method to generate the raw output data
            data = self.calculateCoreSignal(tShiftedByDelay, sweepIndexWithinSet) ;
                % data should be same size as t at this point
            
            % Compute the amplitude from the expression for it
            amplitude = ws.Stimulus.evaluateSweepExpression(self.Amplitude,sweepIndexWithinSet) ;
            % Screen for illegal values
            if isempty(amplitude) || ~(isnumeric(amplitude)||islogical(amplitude)) || ~isscalar(amplitude) || ~isreal(amplitude) || ~isfinite(amplitude) ,
                data=zeros(size(t));
                return
            end

            % Compute the delay from the expression for it
            dcOffset = ws.Stimulus.evaluateSweepExpression(self.DCOffset,sweepIndexWithinSet) ;
            % Screen for illegal values
            if isempty(dcOffset) || ~(isnumeric(dcOffset)||islogical(dcOffset)) || ~isscalar(dcOffset) || ~isreal(dcOffset) || ~isfinite(dcOffset) ,
                data=zeros(size(t));
                return
            end
            
            % Scale by the amplitude, and add the DC offset
            data = amplitude*data + dcOffset;

            % Compute the duration from the expression for it
            duration = ws.Stimulus.evaluateSweepExpression(self.Duration,sweepIndexWithinSet) ;
            % Screen for illegal values
            if isempty(duration) || ~(isnumeric(duration)||islogical(duration)) || ~isscalar(duration) || ~isreal(duration) || ~isfinite(duration) || duration<0 ,
                data=zeros(size(t));
                return
            end
            
            % Zero the data outside the support
            % Yes, this is supposed to "override" the DC offset outside the
            % support.
            isOnSupport=(0<=tShiftedByDelay)&(tShiftedByDelay<duration);
            data(~isOnSupport,:)=0;
            
            if size(data,1)>0 ,
                data(end,:)=0;  % don't want to leave the DACs on when we're done
            end
        end        
        
        function data = calculateCoreSignal(self, t, sweepIndexWithinSet)
            % Compute the period from the expression for it
            period = ws.Stimulus.evaluateSweepExpression(self.Period,sweepIndexWithinSet) ;
            if isempty(period) || ~isnumeric(period) || ~isscalar(period) || ~isreal(period) || ~isfinite(period) || period<=0 ,
                period=nan;  % s
            end
            
            % Compute the period from the expression for it
            pulseDuration = ws.Stimulus.evaluateSweepExpression(self.PulseDuration,sweepIndexWithinSet) ;
            if isempty(pulseDuration) || ~isnumeric(pulseDuration) || ~isscalar(pulseDuration) || ~isreal(pulseDuration) || ~isfinite(pulseDuration) || ...
               pulseDuration<=0 ,
                pulseDuration=nan;  % s
            end
            
            % compute the data
            if isfinite(period) && isfinite(pulseDuration) ,
                timeInCycle=mod(t,period);
                data=double(timeInCycle<pulseDuration);
            else
                % invalid params causes signal to be all-zero
                data=zeros(size(t));
            end
        end        
    end
    
%     %
%     % Implementations of methods needed to be a ws.ValueComparable
%     %
%     methods
%         function value=isequal(self,other)
%             % Custom isequal.  Doesn't work for 3D, 4D, etc arrays.
%             value=isequalHelper(self,other,'ws.SquarePulseTrainStimulusDelegate');
%         end                            
%     end
%     
%     methods (Access=protected)
%        function value=isequalElement(self,other)
%             propertyNamesToCompare={'Delay' 'Duration' 'Amplitude' 'DCOffset' 'Period' 'PulseDuration'};
%             value=isequalElementHelper(self,other,propertyNamesToCompare);
%        end
%     end
    
    methods 
        function out = getPropertyValue_(self, name)
            out = self.(name);
        end  % function
        
        % Allows access to protected and protected variables from ws.Encodable.
        function setPropertyValue_(self, name, value)
            self.(name) = value;
        end  % function
    end
    
    methods
        function mimic(self, other)
            ws.mimicBang(self, other) ;
        end
    end        
    
    methods
        % These are intended for getting/setting *public* properties.
        % I.e. they are for general use, not restricted to special cases like
        % encoding or ugly hacks.
        function result = get(self, propertyName) 
            result = self.(propertyName) ;
        end
        
        function set(self, propertyName, newValue)
            self.(propertyName) = newValue ;
        end           
    end  % public methods block        
    
end
