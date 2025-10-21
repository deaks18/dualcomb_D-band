function [symbolEst, symbolDet, symbolErr, wxx, NLobj] = ...
    MyEqualizeNL(deqObj, NLobj, inputSig, trainSig)
%EQUALIZE   Equalize a signal with an equalizer object.
%  EQUALIZE is a method for an equalizer object.
%  For more information, type 'help equalize' in MATLAB.
% deqObj stands for direct equalizer object hxx, hyy
% xeqObj stands for cross equalizer object hxy, hyx
% eqObj must be an equalizer object for this method to be invoked.

%   Copyright 2003-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2021/07/14 $

%--------------------------------------------------------------------------
% Check number of arguments; if only two, set trainSig to maximum length.
% narginchk(3, 4);

% Check inputSig is numeric.
if isempty(inputSig) || ~isnumeric(inputSig)
    error('comm:equalize:inputSigNumeric', 'INPUTSIG must be numeric.');
end

% Force inputSig to be N by 2 matrix.
sizeSig = size(inputSig);
if sizeSig(2) == 2
    nSamples = size(inputSig,1);
elseif sizeSig(1) == 2
    inputSig = inputSig.';  % transposition, noted should use .' for none conjugate transpose
    nSamples = size(inputSig,1);
else
%     error('comm:equalize:inputSig', 'INPUTSIG stands for duo pol singal, it must be a vector.')
    nSamples = size(inputSig,1);
end

% set trainSig to maximum length, i.e. use all symbols for training
if nargin<4, trainSig=[]; end

%--------------------------------------------------------------------------
% Error checking for inputSig and trainSig.

% Determine number of symbols to estimate.


nSampPerSym = deqObj.nSampPerSym;

% % % nSampPerSym = deqObj.InputSamplesPerSymbol;


nSymbols = ceil(nSamples/nSampPerSym);

% nSamples must be an integer multiple of nSampPerSym.
if nSymbols*nSampPerSym ~= nSamples
   error('comm:equalize:inputSigLength', ...
       ['Number of elements in INPUTSIG must be an integer multiple of '...
        'EQOBJ.nSampPerSym.']);
end

% Force trainSig to be a column vector.
if size(trainSig,1) < size(trainSig,2)
    trainSig = trainSig.';
end
nTrain = size(trainSig,1);   

    % Check trainSig is numeric and shorter than number of output symbols.
if ~isnumeric(trainSig) || nTrain > nSymbols
    error('comm:equalize:trainSigNumeric', ...
          ['TRAINSIG must be numeric vector, '...
            'with length less than or equal to ', ...
            'length(INPUTSIG)/EQOBJ.nSampPerSym.']);
end

% Warn if CMA equalizer is used with training signal.
% if (deqObj.AdaptAlg.BlindMode && nTrain>0) || (xeqObj.AdaptAlg.BlindMode && nTrain>0)
%     warning('comm:equalize:blindmode', ...
%         ['Training signal is ignored for constant modulus algorithm.']);
% end

%--------------------------------------------------------------------------
[symbolEst, symbolDet, symbolErr, wxx, NLobj] = Equalizer.Volterra_SampEqualizer.MyThisEqualizeNL(deqObj, NLobj,inputSig, trainSig);

if sizeSig(2)>sizeSig(1)
    symbolEst = symbolEst.';
    symbolDet = symbolDet.';
    symbolErr = symbolErr.';
end
    