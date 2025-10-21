function [symbolEst, symbolDet, symbolErr, wxx, NLobj] = MyThisEqualizeNL(deqObj, NLobj,inputSig, trainSig)
%THISEQUALIZE   Equalize a signal with an equalizer object.
%  THISEQUALIZE is a method for an equalizer.lineareq object.
%  For more information, type 'help equalize' in MATLAB.

% Copyright 2003 The MathWorks, Inc.

% Most error checking is inherited from baseclass.

% Note that this function is inherited by the DFE as well.  The key
% difference is the external function eqindices.m, which is overloaded for
% the DFE.

%--------------------------------------------------------------------------
% Extract equalizer object parameters.

alg = deqObj.AdaptAlg;
nSampPerSym = deqObj.nSampPerSym;
sigConst = deqObj.SigConst;
blindMode = deqObj.AdaptAlg.BlindMode;
nWeights = deqObj.AdaptAlg.nWeights;


% Derive signal dimensions.
nSamples = length(inputSig);
nSymbols = ceil(nSamples/nSampPerSym);
nTrain = length(trainSig);
delay = floor(deqObj.RefTap-1)/nSampPerSym; % delay for the output symbol
xin = flipud(reshape(inputSig(:,1), [nSampPerSym nSymbols]));% polarisation, 
% yin = flipud(reshape(inputSig(:,2), [nSampPerSym nSymbols]));


%%%%%%%%%%%%%%% added by JWei %%%%%%%%%%%%%%
EQ.tap = [deqObj.nWeights, NLobj.NLtap(2), NLobj.NLtap(3)];             % memory length of the linear, 2nd and 3rd order terms. Note the latter two values must be not larger than the first one
EQ.mu  = [deqObj.StepSize, NLobj.Step2nd, NLobj.Step3rd];                % step siye of the linear, 2nd and 3rd order terms
% if EQ.tap(2) ~= 0
%     Hxx_NL  = zeros(EQ.tap(2),EQ.tap(2));
%     Hyx_NL  = zeros(EQ.tap(2),EQ.tap(2));
%     Hyyxx_NL  = zeros(EQ.tap(2),EQ.tap(2));
% else Hxx_NL = 0; Hyx_NL = 0;  Hyyxx_NL = 0; 
% end
% if EQ.tap(3) ~= 0
%     Hxx_3NL  = zeros(EQ.tap(3),EQ.tap(3),EQ.tap(3));
%     Hyx_3NL  = zeros(EQ.tap(3),EQ.tap(3),EQ.tap(3));
% else Hxx_3NL = 0; Hyx_3NL = 0;  
% end
%%%%%%%%%% add ended %%%%%%%%%%%%%%

% nTrain should be shorter than nSymbols -
% ceil(delay/nSampPerSym), else it will be truncated
if nTrain > nSymbols - delay
    warning('comm:equalize:trainSigLength', ...
          ['TRAINSIG will be truncated to a length equal to '...
           'length(INPUTSIG)- delay, where delay = '...
           '(EQOBJ.RefTap-1)/EQOBJ.nSampPerSym.']);
end

% Initialize symbol estimates and errors.
symbolEst = zeros(nSymbols, 1);    % symbol estimate vector
symbolDet = zeros(nSymbols, 1);    % detected symbols vector
symbolErr = zeros(nSymbols, 1);    % symbol error vector

% Reset channel if required.
if deqObj.ResetBeforeFiltering
    reset(deqObj);
end
% if xeqObj.ResetBeforeFiltering
%     reset(xeqObj);
% end

% Indices for input buffer, u, and decision sample, d.
% - i1 is the input buffer index for the current input signal sample:
%      u(i1)=sample; see below
% - nn1 and nn2 are input buffer indices used to perform the 
%   sample shifting operation each equalizer iteration: 
%      u(nn2)=u(nn1); see below.
% - id is an "index" into the decision value: 
%      decision_feedback_value=d(id); id=[] for a linear equalizer; 
%      id=1 for DFE.
[i1, nn1, nn2, id] = eqindices(deqObj);

% Initialize input buffer, accounting for reference tap delay.
% Get current equalizer weights.
wxx = deqObj.Weights.';  % Use hermitian transpose for mathematical simplicity. 
uxx = deqObj.WeightInputs.';
% wxy = xeqObj.Weights.';  
% uxy = xeqObj.WeightInputs.';

%%%%%%%%% added by JWei %%%%%%%%%
Hxx_NL    =  NLobj.Hxx_NL;
% Hyx_NL    =  NLobj.Hyx_NL;
% Hyyxx_NL  = NLobj.Hyyxx_NL;
Hxx_3NL   = NLobj.Hxx_3NL;
% Hyx_3NL   = NLobj.Hyx_3NL; 

% Constellation.
if nTrain >= 1
    constParam = sigConst * sqrt(mean(abs(trainSig).^2)/mean(abs(sigConst).^2));
else
    constParam = sigConst;
    warning('no training signal used, should use CMA mode');
end
nconst = numel(constParam);

% Initial decision value (required for compatibility with DFE).
dref = deqObj.dref;

if blindMode
    % if var(abs(sigConst))>eps
    %     warning('comm:equalize:CMAsigconst',...
    %         ['Signal constellation with non-constant magnitude '...
    %          'may cause CMA equalizer to be unstable.'])
    % end
    % Set CMA constant
    % R2 = mean(abs(sigConst).^4)/mean(abs(sigConst).^2);
    R2Vec = unique(real(sigConst).^2+imag(sigConst).^2); % the number of levels
end

%  Main loop.
for n = 1:nSymbols
    % the direct filter
    uxx(nn2) = uxx(nn1);             % shift input signal buffer samples
    uxx(i1) = [xin(:, n); dref(id)]; % current input signal vector
%     % the cross filter
%     uxy(nn2) = uxy(nn1);             % shift input signal buffer samples
%     uxy(i1) = [yin(:, n); dref(id)]; % current input signal vector
    %%%%%%%%%%%%%%%%%%%% added by JWei %%%%%%%%%%%%%%%%%%%%%
    if EQ.tap(2) ~= 0
        Ex_reg_NL(1:EQ.tap(2),1)  = uxx((EQ.tap(1)-EQ.tap(2))/2+1:end - (EQ.tap(1)-EQ.tap(2))/2, 1);  
        Ex_reg_NL                 = tril(Ex_reg_NL*Ex_reg_NL.',-1);
%         Ey_reg_NL(1:EQ.tap(2),1)  = uxy((EQ.tap(1)-EQ.tap(2))/2+1:end - (EQ.tap(1)-EQ.tap(2))/2, 1);  
%         Ey_reg_NL                 = tril(Ey_reg_NL*Ey_reg_NL.',-1);
%         Exy_reg_NL                = tril(Ex_reg_NL*Ey_reg_NL.',0);
    else
        Ex_reg_NL = 0; 
%         Ey_reg_NL = 0;  
%         Exy_reg_NL = 0;
    end
    if EQ.tap(3) ~= 0
        Ex_reg_NL3(1:EQ.tap(3),1)  = uxx((EQ.tap(1)-EQ.tap(3))/2+1:end - (EQ.tap(1)-EQ.tap(3))/2, 1);  
        Ex_reg_3NL                 = zeros(EQ.tap(3),EQ.tap(3),EQ.tap(3));
%         Ey_reg_NL3(1:EQ.tap(3),1)  = uxy((EQ.tap(1)-EQ.tap(3))/2+1:end - (EQ.tap(1)-EQ.tap(3))/2, 1);  
%         Ey_reg_3NL                 = zeros(EQ.tap(3),EQ.tap(3),EQ.tap(3));
        for k1 = 1:length(Ex_reg_NL3)
            Ex_reg_3NL(:,:,k1)     = tril(Ex_reg_NL3*Ex_reg_NL3.')*conj(Ex_reg_NL3(k1));
%             Ey_reg_3NL(:,:,k1)     = tril(Ey_reg_NL3*Ey_reg_NL3.')*conj(Ey_reg_NL3(k1));
        end
    else
        Ex_reg_NL3 = 0; 
%         Ey_reg_NL3 = 0;
    end
    % output
%     out = wxx'*uxx + wxy'*uxy + sum(sum(Ex_reg_NL.*Hxx_NL)) + sum(sum(Ey_reg_NL.*Hyx_NL)) + sum(sum(Exy_reg_NL.*Hyyxx_NL)) + ...
%           sum(sum(sum(Hxx_3NL.*Ex_reg_3NL))) + sum(sum(sum(Hyx_3NL.*Ey_reg_3NL))); 
      out = wxx'*uxx + sum(sum(Ex_reg_NL.*Hxx_NL)) + sum(sum(sum(Hxx_3NL.*Ex_reg_3NL))); 
      %%%%%%%%%%%%%%%%%%%% add ended %%%%%%%%%%%%%%%%%%%%%
    
    % the linear term output 
%     out = wxx'*uxx + wxy'*uxy;       % current output symbol (note hermitian transpose)
    
    % decision: find closest signal constellation point (detector/slicer)
    [~, idx] = min(abs(constParam - repmat(out,1,nconst)));
%     [minMetric, idx] = min(constParam - abs(out)^2); % this for CMA
    d = constParam(idx);

    % Symbol error for training/blind mode or decision-directed mode.
    if ~blindMode
        m = n - delay;
        if m >= 1 && m <= nTrain
            dref = trainSig(m);  % training mode
        else
            dref = d;  % decision-directed mode
        end
        e = dref - out; % error 
    else % blind mode
        dref = d;  % required for DFE
        rawErrVec = R2Vec - abs(out).^2;
        [~, errIdx] = min(abs(rawErrVec));
        e = out .* rawErrVec(errIdx); % blind mode (CMA)
    end
    if isinf(e) || isnan(e)
        warning('comm:equalize:error',...
           ['Equalizer unstable. ' ...
            'Try adjusting step size or forgetting factor.']);
    end

    % Update weights until exhausted input samples. The weight update
    % algorithm, update_weights, is a UDD method. That is, the method
    % called will depend on the class of adaptive algorithm object.
    wxx = update_weights(alg, wxx, uxx, e);
%     wxy = update_weights(alg, wxy, uxy, e);
    %%%%%%%%%%%%%%%%%%% added by JWei %%%%%%%%%%%%%%
    M = 1; % LeakageFactor, default 1 means no leakage
    Hxx_NL    = M*Hxx_NL + EQ.mu(2)*e*conj(Ex_reg_NL);
%     Hyx_NL    = M*Hyx_NL + EQ.mu(2)*e*conj(Ey_reg_NL);
%     Hyyxx_NL  = M*Hyyxx_NL + EQ.mu(2)*e*conj(Exy_reg_NL);
    Hxx_3NL   = M*Hxx_3NL + EQ.mu(3)*e*conj(Ex_reg_3NL);
%     Hyx_3NL   = M*Hyx_3NL + EQ.mu(3)*e*conj(Ey_reg_3NL);
    %%%%%%%%%%%%%%%%% add ended %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % debug --- evolve of equalizer tap coefficients
%     plot(1:nWeights, real(wxx),'bo', 1:nWeights, real(wxy),'r^'); drawnow;
%     pause(0.01);
%     scatter(real(e),imag(e),'b');axis([-1.5 1.5 -1.5 1.5]); hold on; pause(0.01)
    % Assign outputs and symbol error.l
    symbolDet(n) = d;
    symbolErr(n) = e;
    symbolEst(n) = out;
end

%--------------------------------------------------------------------------
% Save equalizer weights, weight inputs, and reference sample.
deqObj.Weights = wxx.';  % note hermitian transpose
% deqObj.WeightInputs = uxx.';   % note physical transpose
% eqObj.dref = dref;
% xeqObj.Weights = wxy.';  % note hermitian transpose
% xeqObj.WeightInputs = uxy.';   % note physical transpose
%%%%%%%%%%% added by JWei %%%%%%%%%%%%%%%%5
NLobj.Hxx_NL   = Hxx_NL;
% NLobj.Hyx_NL   = Hyx_NL;
% NLobj.Hyyxx_NL = Hyyxx_NL;
NLobj.Hxx_3NL  = Hxx_3NL;
% NLobj.Hyx_3NL  = Hyx_3NL;
%%%%%%%%%%% add ended %%%%%%%%%%%5

% Update number of samples processed by equalizer.
% eqObj.NumSamplesProcessed = eqObj.NumSamplesProcessed + nSamples;
