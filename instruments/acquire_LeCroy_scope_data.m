function output = acquire_LeCroy_scope_data(ip_address, N)

% INPUTS
% ip_address -> static ipv4 address of scope
% N -> number of points to capture

% OUTPUTS
% output -> Nx8 array of waveforms, from all 8 scope channels.


% This script uses a modified version of the lecroy driver "lecroy_basic_driver_8Ch.mdd"
%
% The original driver can be downloaded from 
% https://uk.mathworks.com/matlabcentral/fileexchange/12862-lecroy-oscilloscopes-matlab-instrument-driver
% which only works up to 4 channels.

% Note that only the "readwaveform" function has been modiefied to allow
% for 8 channels. Other driver functions, that are not used in this script,
% may not work when calling channel numbers > 4.

% Callum Deakin 29/11/2022

interfaceObj = tcpip(ip_address, 1861);

% Create a device object. 
deviceObj = icdevice('lecroy_basic_driver_8Ch.mdd', interfaceObj);

% Connect device object to hardware.
connect(deviceObj);

% Execute device object function(s).
% groupObj = get(deviceObj, 'Acquisition');
% groupObj = groupObj(1);
% set(groupObj, 'State', 'single')

groupObj = get(deviceObj, 'Trigger');
groupObj = groupObj(1);
set(groupObj, 'Source', 'acline'); % trigger off AC rising edge cos I can't find the function for instant capture
invoke(groupObj, 'trigger');

groupObj = get(deviceObj, 'Waveform');
groupObj = groupObj(1);
set(groupObj, 'MaxNumberPoint',N-1);
set(groupObj, 'Precision','int16');
output = zeros(8,N);

for i = 1:8
    [output(i,:), ~, ~, ~, ~] = invoke(groupObj, 'readwaveform', strcat('channel',string(i)));    
end
output = output.';

% Disconnect device object from hardware.
disconnect(deviceObj);

% Delete objects.
delete([deviceObj interfaceObj]);
