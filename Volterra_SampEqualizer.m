classdef Volterra_SampEqualizer
    % dual-pol Volterra - decision direct Summary of this class goes here
    % two by one FFE, use for single-pol modulation but dual pol detection
    % created by Zhixin Liu @ 2021 zhixin.liu@ucl.ac.uk
    % edited by Jinlong Wei, August 2021: moving NLtaps and step size into properties
    
    properties
        name;               % name
        trainlen;		    % number of train symbols
        NLtap;              % memory length of 2nd and 3rd order FFE added by JWei
        NLmu;               % step size of 2nd and 3rd order terms added by JWei
    end
    
    methods
        
        function obj = Volterra_SampEqualizer(varargin)
            % Check number of input arguments
            narginchk(0,4);
            % Set default values, default is a linear equalizer
            fieldArray = {'NLtap', 'NLmu','trainlen', 'name'};   
            defaltVal = containers.Map( ...
                fieldArray, ...            
                {9 , 1e-4*[5 1 0.5], 1e4, 'Volterra single equalizer'});
            % Set the fields
            %obj = OOPUtility.ComFun.Init(obj, fieldArray, varargin, defaltVal);
        end
        
        function [symbolEst, symbolDet] = Equalize(obj, samp, refsym, const, nSampPerSym)
            % check the polariztion of the input signal
%             if size(samp,2) ~=2
%                 warning('the input signal is single-polarized');
%                 samp = [samp zeros(numel(samp),1)];
%             end
            %% preparation work
            % adjust the input sample sequence when necessary
            if (nSampPerSym) == 2 && bitand(obj.NLtap(1), 2) % when the tap number is odd, adjust the tap to ensure the center corresponds to the output
                samp = circshift(samp, 1);
            end
            
            if numel(refsym) < obj.trainlen
               error('reference symbol length smaller than training signal length: not enough training symbols');
            else
               trainSig = refsym(1:obj.trainlen);
            end
            
            % does window needed here to prevent the frequency leakage?
            %% the FFE equalizer            
            % create the butterfly cma equalizer
           
            hxx = lineareq(obj.NLtap(1), lms(obj.NLmu.scale1st(1)), const, nSampPerSym);
           %%% hxx = comm.LinearEqualizer('Algorithm','LMS','NumTaps',obj.NLtap(1),'StepSize',obj.NLmu.scale1st(1));
            % and adjust its parameters
           
            delay = round(((1+obj.NLtap(1))/2-1)/nSampPerSym); % number of taps delay for output
%             hxx.Weights(1:end) = 0; hxx.Weights(nSampPerSym*delay+1) = 1+1i; % center tap 1, others zero
      
           hxx.ResetBeforeFiltering = 0;
           hxx.RefTap = floor((obj.NLtap(1)+1)/2/nSampPerSym)*nSampPerSym + 1;
          
            % and adjust its parameters
             
            

          
            % define the cross connect filter hyy, note it needs a new definition, cannot simply
            % pass hxx=hyy (they will use the same property if define this way)
%             hxy = lineareq(obj.NLtap(1), lms(obj.NLmu.scale1st(1)), const, nSampPerSym);
%             hxy.ResetBeforeFiltering = 0;
%             hxy.RefTap = hxx.RefTap;
             
            %%%%%%%%%%%%%%%% nonlinear terms %%%%%%%%%%%%%
            NLobj.NLtap = obj.NLtap;
            if obj.NLtap(2) ~= 0
                NLobj.Hxx_NL  = zeros(NLobj.NLtap(2),NLobj.NLtap(2));
%                 NLobj.Hyx_NL  = zeros(NLobj.NLtap(2),NLobj.NLtap(2));
%                 NLobj.Hyyxx_NL  = zeros(NLobj.NLtap(2),NLobj.NLtap(2));
            else
                NLobj.Hxx_NL = 0; 
%                 NLobj.Hyx_NL = 0;  
%                 NLobj.Hyyxx_NL = 0;
            end
            if obj.NLtap(3) ~= 0
                NLobj.Hxx_3NL  = zeros(NLobj.NLtap(3),NLobj.NLtap(3),NLobj.NLtap(3));
%                 NLobj.Hyx_3NL  = zeros(NLobj.NLtap(3),NLobj.NLtap(3),NLobj.NLtap(3));
            else
                NLobj.Hxx_3NL = 0; 
%                 NLobj.Hyx_3NL = 0;
            end
            %%%%%%%%%%%%%%%% add ended %%%%%%%%%%%%%%%%%%%
            
            % defein hyy and hyx if need dual pol output. 
            %% equalization
            % training the taps before equalization
            for idx = 1:numel(obj.NLmu.scale1st)-1
                % change the stepsize of the equalizer
                hxx.StepSize = obj.NLmu.scale1st(idx);
                NLobj.Step2nd = obj.NLmu.scale2nd(idx);    % added by JWei
                NLobj.Step3rd = obj.NLmu.scale3rd(idx);    % added by JWei
                % train the equalizer
                [~,~,~,wxx,NLobj] = Equalizer.Volterra_SampEqualizer.MyEqualizeNL(hxx, NLobj, samp, trainSig(1:obj.trainlen));
%                 [~,~,~,~,~,NLobj] = Equalizer.Volterra_SampEqualizer.MyEqualizeNL(hxx, hxy, NLobj, samp, trainSig(1:obj.trainlen));
            end
            % perform equalization
            hxx.StepSize = obj.NLmu.scale1st(end);       % linear terms
            NLobj.Step2nd = obj.NLmu.scale2nd(end);    % 2nd order terms
            NLobj.Step3rd = obj.NLmu.scale3rd(end);    % 3rd order terms
            hxx.WeightInputs = wxx.';
            [symbolEst, symbolDet, symbolErr] =  Equalizer.Volterra_SampEqualizer.MyEqualizeNL(hxx, NLobj, samp, trainSig(1:obj.trainlen)); 
            % correct the delay
            symbolDet = symbolDet((obj.trainlen+1+delay):end);
            symbolEst = symbolEst((obj.trainlen+1+delay):end);
        end
        
        function cellOfproperties = Property2cell(obj)
            cellOfproperties = {...
                'Number of FFE taps', obj.nFFETaps, '';
                'Step size vector', num2str(obj.stepSizeVec), '';
                'NLtap', obj.NLtap, '';
                'NLmu', obj.NLmu, '';
                };
        end
    end
    
    methods (Static = true)
%         [symbolEst, symbolDet, symbolErr, wxx, wxy, NLobj] = MyEqualizeNL(deqObj, xeqObj, NLobj, inputSig, trainSig);
%         [symbolEst, symbolDet, symbolErr, wxx, wxy, NLobj] = MyThisEqualizeNL(deqObj, xeqObj, NLobj, inputSig, trainSig);     
        [symbolEst, symbolDet, symbolErr, wxx, NLobj] = MyEqualizeNL(deqObj, NLobj, inputSig, trainSig);
        [symbolEst, symbolDet, symbolErr, wxx, NLobj] = MyThisEqualizeNL(deqObj, NLobj, inputSig, trainSig);
    end
    
end