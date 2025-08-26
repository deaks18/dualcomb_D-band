% ****************************************************************************
% MODULE:	  IQ2SmuWv.m
% COPYRIGHT:  (c) 2004 Rohde & Schwarz GmbH & Co. KG, Munich
% PROJECT:    Demo Scripts
% LANGUAGE:   MATLAB Version 6.1.0.450 (R12.1)
% DATE:       22-SEP-2004
% AUTHOR:     Johannes Ganzert
%
% HISTORY:
%             2004-09-04: Initial revision
% ****************************************************************************
%
% 
% ****************************************************************************
% FUNCTION DEFINITIONS
% ****************************************************************************
%
function rs_smu_iq2wv(data, filename, SampleRate, varargin)
%
% function IQ2SmuWv(data,filename,[SampleRate])
%
% SPECIFICATION: The function IQ2SmuWv creates a waveform file for the
%                SMU/SMW/SGT from the complex vector IQ. The sample rate
%                is an optional argument.
%
% PARAMETERS:	 INPUT:
%                   data         : I/Q samples (complex)
%                   filename     : Name der WV-Datei, die erzeugt
%                                  werden soll.
%                   SampleRate   : sample rate for IQ waveform in baseband ARB
%                   [Tburst]     : burst period in samples for Crest Factor
%                                  calculation
%                   
%
%                OUTPUT: None
%
% PRECONDITIONS: None
% 
% SIDE_EFFECTS:	 None
% 
% RETURN VALUES: None
% 
% EXCEPTIONS: 	 None
%
% AUTHOR:        
%
%
% ****************************************************************************

% Check number of arguments
if ~(nargin == 3 | nargin == 4)
   error('Wrong number of input parameters.');
end

if nargin == 4
    Tburst = varargin{1};
else
    Tburst = [1 length(data)];
end

% check FILENAME for extension
if ~strcmp( lower(filename(end-2:end)), '.wv' )
   filename = [filename '.wv'];
end

% check size of DATA
data_length = length(data);
if data_length < 512
   error('The waveform must contain a minimum of 512 I/Q samples.');
end
if data_length > 256e6
   error('The maximum waveform length is 256 M samples');
end
if mod(data_length,4) ~= 0
   error('Die Anzahl der I/Q-Samples muﬂ durch 4 teilbar sein.');
end

% fill in imaginary part if DATA is not complex
if isreal(data)
   data = complex(data);
end

% reorganize real and imag part of DATA as columns
if size(data,1) < size(data,2)
   data = data.';
end

% decompose DATA in real and imaginary part
data_real = real(data);
data_imag = imag(data);

%Compute the Crest Factor 
%for bursted signals Tburst defines the burst period
%for continuous signals Tburst = full signal duration
idx_nz = Tburst(1):Tburst(2); % burst period, magnitude > 0
fPeak = max(abs(data(idx_nz)));
fRMS  = sqrt(mean(data_real(idx_nz).^2 + data_imag(idx_nz).^2));
fCrestFactorLog = 20*log10(fPeak/fRMS);   
clear data;

%normalize the peak to -1,+1 in order to get the correct signal level.
%  I.e. the peak must be the fullscale of the D/A converter.
scaling = max(abs([data_real;data_imag]));
data_real = data_real/scaling;
data_imag = data_imag/scaling;
clear scaling;


if 0
  % IMPORTANT NOTE:
  %   Normalization must be done all the time in order to get the
  %   SMU Signal Level right. I.e. the peak must be the 
  %   fullscale of the D/A converter.
  % normalize or limit DATA to - 1 ... +1
  if max(abs(data_real)) > 1 | max(abs(data_imag)) > 1
     resp = [];
     while ~(strcmp(lower(resp),'n') | strcmp(lower(resp),'l'))
        clc;
        disp('The I and Q samples must be in the range -1 ... +1');
        disp('This data contains samples that are out of range');
        resp = input('Do you want to (n)ormalize or (l)imit the data?  n/l [n]: ','s');
        if isempty(resp)
           resp = 'n';
        end
     end
     switch lower(resp)
        case 'n'
           scaling = max(abs([data_real;data_imag]));
           data_real = data_real/scaling;
           data_imag = data_imag/scaling;
           clear scaling;
        case 'l'
           data_real(find(data_real > 1)) = 1;
           data_real(find(data_real < -1)) = -1;
           data_imag(find(data_imag > 1)) = 1;
           data_imag(find(data_imag < -1)) = -1;
     end
  end
  clear resp;
end;

% Quantization, cf. SMU manual, "Waveform File Format".
%
%   IQ complex     IQ_quantized
%   ----------------------------
%       -1.0         -32767
%        0.0              0
%       +1.0          32767
data_real = round(32767*data_real);
data_imag = round(32767*data_imag);
data_IQ = reshape([data_real data_imag]',1,2*data_length);
clear data_real data_imag;

% Write waveform file
file_id = fopen(filename,'w');
fprintf(file_id,'%s','{TYPE: SMU-WV,0}');
fprintf(file_id,'%s',['{CRESTFACTOR: ' num2str(fCrestFactorLog) '}']);
fprintf(file_id,'%s','{COMMENT: No Comment}');
fprintf(file_id,'%s',['{DATE: ' datestr(clock,0) '}']);
fprintf(file_id,'%s',['{LEVEL OFFS: ' num2str(fCrestFactorLog) ', ' num2str(0.0) '}']);
fprintf(file_id,'%s',['{CLOCK: ' num2str(SampleRate) '}']);
fprintf(file_id,'%s',['{SAMPLES: ' num2str(data_length) '}']);
fprintf(file_id,'%s',['{WAVEFORM-' num2str(4*data_length+1) ':#']);
fwrite(file_id,data_IQ,'int16');
fprintf(file_id,'%s','}');
fclose(file_id);

% end of function