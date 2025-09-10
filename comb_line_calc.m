clear; close all;

sig_comb_freq = 25e9;
LO_comb_freq = 27e9;

mixer_freq = 110e9:0.1e9:170e9;
LO_channel = round(abs(mixer_freq - sig_comb_freq*round(mixer_freq/sig_comb_freq))/(LO_comb_freq-sig_comb_freq));
sig_channel_pos = abs(LO_channel + (-1).^(round(mixer_freq/(2.*sig_comb_freq))) .* round(mixer_freq/sig_comb_freq));
sig_channel_neg = abs(LO_channel - (-1).^(round(mixer_freq/(2.*sig_comb_freq))) .* round(mixer_freq/sig_comb_freq));

figure; plot(mixer_freq/1e9,LO_channel);
hold on; plot(mixer_freq/1e9,sig_channel_pos);
hold on; plot(mixer_freq/1e9,sig_channel_neg);

legend('LO','Sig pos','Sig neg');

figure; plot(mixer_freq/6e9,LO_channel);
hold on; plot(mixer_freq/6e9,sig_channel_pos);
hold on; plot(mixer_freq/6e9,sig_channel_neg);

legend('LO','Sig pos','Sig neg');

%% write to csv for paper

T = table(mixer_freq.'/1e9, LO_channel.', sig_channel_pos.', sig_channel_neg.', 'VariableNames', {'freq', 'lo', 'pos', 'neg'});

writetable(T, 'channel_assignment.csv');