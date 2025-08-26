%(c)2021 UCL employee of the year Dhecha Nopchinda d.nopchinda@ucl.ac.uk dhecha@ieee.org
function sendToSMW_DN_241121(iq,fs,fc,RMSin,ip_SMW)
%% Note
%This is a function to send I and Q vectors to the R&S SMW Signal Generator 
%iq is the IQ signal vector. fs is the sampling rate (Hz), fc is the center
%frequency (Hz), RMSin is the RF powers setting.
%The IQ signal vector will be normalized such that the dynamic range of the
%DAC is maximized and clipping is avoided.

%If error occurred whilst connected to the instrument. You must force a
%disconnect by running the following command: instrreset

%ver0.1 241121
%% PA
%warning('"Remember that it is the peak power, not the average, that kills an instrument." - Dhecha the voice of reason')
%% Check
if (fs>2.4e9)||(fs<400)
    error('Dhecha the voice of reason: sampling rate out of range. Operation terminated.')
end
if (fc>44e9)||(fc<100e3)
    error('Dhecha the voice of reason: center frequency out of range. Operation terminated.')
end
if (RMSin>30)||(RMSin<-145)
    error('Dhecha the voice of reason: assigned power out of range. Operation terminated.') 
end
N=length(iq);
if mod(N,4)~=0
    error('Dhecha the voice of reason: Length of the IQ signal must be divisible by 4. Operation terminated.')
end
if (N>256e6)||(N<512)
    error('Dhecha the voice of reason: Length of the IQ signal out of range. Operation terminated.') 
end
if ~isvector(iq)
    error('Dhecha the voice of reason: The IQ signal is not a vector. Operation terminated.') 
end
%% Param
if size(iq)>1
    iq=iq.';%convert to row vector
end
% fs=100e6;%400 Hz to 2.4 GHz
% fc=35e9; %RF frequency
% RMSin=-10;%RF output power in dBm
% N=length(iq);%Must be divisible by 4, max length 56M (by rs_smu_iq2wv not the instrument), min length 512
% timeVect=0:1/fs:(N-1)/fs;
% iq=exp(1i*2*pi*3e6*timeVect);
wv_name='tempWaveFile';%No need to modify, this is just a temp handler.
%% Prep IQ vectors in a format recognized by the instrument
rs_smu_iq2wv(iq, wv_name, fs) %convert the IQ samples into a .wv fileformat (.wv is the fileformat understood by the SMW ARB)
%% Open connection
if N<65e6
    SMW = visa('rs',['TCPIP::' ip_SMW '::hislip0'],'OutputBufferSize',100e6,'Timeout',20);
else 
    SMW = visa('rs',['TCPIP::' ip_SMW '::hislip0'],'OutputBufferSize',100e6,'Timeout',90);
end
disp('Contacting SMW, please wait...')
%%The OutputBufferSize property specifies the maximum number of bytes 
%%that can be written to the instrument at once. 
%%By default, OutputBufferSize is 512 bytes.
pause(1)
fopen(SMW);
pause(10)
%Turn off RF - DN
query(SMW,':OUTPut1:STATe 0; *OPC?');
%%append *OPC? to use query instead of write. This effectively forces the
%%command to finish before proceeding.
%% Baseband
disp('Agreeing on baseband settings...')
%Upload signal to instrument
upload_iq_tar(SMW,[wv_name '.wv'],[wv_name '.wv'])
pause(5)
%Copy the uploaded signal to a local directory on SMW / Set waveform source
%to the uploaded file.
query(SMW,sprintf(':SOURce1:BB:ARBitrary:WAVeform:SELect ''/var/user/%s''; *OPC?',[wv_name '.wv']));
%Set sampling rate of the ARB
query(SMW,[':SOURce1:BB:ARBitrary:CLOCk ' num2str(fs) '; *OPC?']);
% Output auxillary trigger with Marker 1 Output. Triggered each time the
% ARB restarts (uncomment to enable)
% query(SMW,':SOURce1:BB:ARBitrary:TRIGger:OUTPut1:MODE REST; *OPC?');
%Switch BB mode to ARB and turn on BB component
query(SMW,':SOURce1:BB:ARBitrary:STATe 1; *OPC?');
%% RF
disp('Agreeing on RF settings...')
%Set Reference CLK
query(SMW,':SOURce1:ROSCillator:SOURce INT; *OPC?');% Replace EXT with INT for internal.
query(SMW,'SOURce1:ROSCillator:EXTernal:FREQuency 10MHZ; *OPC?');
%Set center frequency
query(SMW,[':SOURce1:FREQuency:CW ' num2str(fc) '; *OPC?']);
%Set output power
query(SMW,[':SOURce1:POWer:LEVel:IMMediate:AMPLitude ' num2str(RMSin) '; *OPC?']);
%Turn on RF - DN
query(SMW,':OUTPut1:STATe 1; *OPC?');
%% Close connection
%terminate connection
pause(1)
fclose(SMW);
delete(SMW);
disp('Handover to SMW completed. RF stage engaged.')
% delete(SMW);