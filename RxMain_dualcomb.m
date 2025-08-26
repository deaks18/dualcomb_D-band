clear; close all;

folder = fileparts(which(mfilename)); 
addpath(genpath(strcat(folder,'\..')));

ParaSetting_DualComb;

%% Parameters
% computer is set as 192.168.0.1
ip_lecroy = '192.168.0.100';
ip_SMW = '192.168.0.2';  
bit_sequence_length = 2^14;
constellation_points = 16;
fs_DAC = 2.4e9;
baud_rate = SC.SymbolRates;
center_freq = 60e6;
power = -11;

mixer_freq = 25.97e9*6;
comb_freq = 26e9*6;

comb_diff = comb_freq - mixer_freq;

fs_scope = 10e9;
comment = "_4_no_eravant";

file_location = strcat("C:\\Users\\",getenv('username'),"\\OneDrive - Nokia\\ExperimentData\\2025-08-25_UCL_dualcomb\\");
if ~exist(file_location, 'dir')
   mkdir(file_location)
end
filename = strcat('DualComb_',....
                  string(center_freq*1e-6),'MHz_',...
                  string(power),'dBm_',...
                  string(int64(baud_rate/1e6)),'MBd',...
                  string(constellation_points),'-QAM',...
                  comment);
              
%% TxSig
bitSeqVal = nrPRBS(1,bit_sequence_length);
constSeq = qammod(int8(bitSeqVal),constellation_points,'InputType','bit','UnitAveragePower',true);

if FLAG.if_CreateDACFiles 
    filter_span_symbols = 50;
    txfilter = comm.RaisedCosineTransmitFilter('Shape','Square root', ... % RRC
                                               'RolloffFactor',0.1, ...
                                               'FilterSpanInSymbols',filter_span_symbols, ...
                                               'OutputSamplesPerSymbol',2);
    txSig = txfilter(constSeq);
    % txSig = txSig((2*filter_span_symbols):end);
    % constSeq = constSeq(filter_span_symbols:end);
    txSig = resample(txSig,fs_DAC,2*baud_rate);
    
    writematrix(txSig,strcat(file_location,"Tx_",filename));
    sendToSMW200A(txSig,fs_DAC,center_freq,power,ip_SMW);
else
    txSig = readmatrix(strcat(file_location,"Tx_",filename));
end

%% Scope
if FLAG.if_CaptureNew
    output = acquire_LeCroy_scope_data(ip_lecroy, 2*length(txSig)*fs_scope/fs_DAC);
    writematrix(output,strcat(file_location,"Rx_",filename));

else
    output = readmatrix(strcat(file_location,"Rx_",filename));
end

rxSig = output(:,2);
spectrum_plot(rxSig,fs_scope);

%% shift to baseband
T = (0:(length(rxSig)-1)).*(1/fs_scope);
freq_shift = exp(1i*2*pi*(-center_freq-comb_diff)*T).';
rxSig = freq_shift.*rxSig;

% rxSig = conj(rxSig);
spectrum_plot(rxSig,fs_scope);

txSig_resampled = resample(txSig,fs_scope,fs_DAC);
synced_rx = SignalSync(rxSig, txSig_resampled);

% SNR filter
filter_BW = 28e6;
N = length(synced_rx);
resolution = fs_scope/N;
freqs = (-fs_scope/2:resolution:(fs_scope/2-resolution)).';
SNR_filter = freqs > -filter_BW & freqs < filter_BW;
SNR_filter = fftshift(SNR_filter);

synced_rx_cut_snr = ifft(SNR_filter.*fft(synced_rx));
txSig_resampled_snr = ifft(SNR_filter.*fft(txSig_resampled));

synced_rx_cut_snr = synced_rx_cut_snr/sqrt(sum(abs(synced_rx_cut_snr.^2)));
txSig_resampled_snr = txSig_resampled_snr/sqrt(sum(abs(txSig_resampled_snr.^2)));

% mean_angle = 1;
% while abs(mean_angle) > 1e-15
%     [synced_rx_cut_snr,mean_angle] = mean_phase_remove(synced_rx_cut_snr,txSig_resampled_snr);
% end

[freqs,PSD] = doublesided_PSD(synced_rx_cut_snr,fs_scope);
figure; plot(freqs,PSD); hold on;
[freqs,PSD] = doublesided_PSD(txSig_resampled_snr,fs_scope);
plot(freqs,PSD); 

SNR_signal = calculateSNR(synced_rx_cut_snr,txSig_resampled_snr);

[r,lags] = xcorr(synced_rx,txSig_resampled);
figure; plot(lags,abs(r));

%% Go to regular DSP
rxSig_resampled = resample(rxSig,fs_DAC,fs_scope);
spectrum_plot(rxSig_resampled,fs_DAC);

indexes_xcorr = FrameSync(rxSig_resampled, txSig, FLAG, EQUIPMENT);
DATA_MOD = SC;
DATA_MOD.indexes = indexes_xcorr;
DATA_MOD.DAC_FrameSize = length(txSig);

[Report]=RX_DSP_SC_Main(rxSig_resampled,constSeq,DATA_MOD,EQUIPMENT,FLAG);
