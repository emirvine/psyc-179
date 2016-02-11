cd('C:\Users\Emily\Desktop\R042-2013-08-18');
cfg_csc.fc = {'R042-2013-08-18-CSC03a.ncs'};
csc = LoadCSC(cfg_csc);

fs = csc.cfg.hdr{1}.SamplingFrequency;

%%


csc_restricted = restrict(csc, 3282, 3286);
plot(csc_restricted.tvec, csc_restricted.data);

% [s, freq, time, power] = spectrogram(csc_restricted.data, rectwin(512), 384, 1:200:0.01, fs);

[s, freq, time, power] = spectrogram(csc_restricted.data, hanning(256), 128, 1:200, fs);
subplot(211);
imagesc(time, freq, 10*log10(power));
[s, freq, time, power] = spectrogram(csc_restricted.data, hanning(1024), 384, 1:200:0.002, fs);
subplot(212);
imagesc(time, freq, 10*log10(power));
set(gca, 'Fontsize', 14);
axis xy; xlabel('Time (s)');
ylabel('Frequency (Hz)');
hold on;

lfp_minmax = 25;
lfp_location = 125;
t0 = csc_restricted.tvec - csc_restricted.tvec(1);
data = rescale(csc_restricted.data, -lfp_minmax, lfp_minmax);
data = data + lfp_location;

lfp_h = plot(t0, data, 'k');

%% Pitfalls: gaps
cscR = restrict(csc, 3300, 3340);
 
[S,F,T,P] = spectrogram(cscR.data, rectwin(256), 128 ,1:200, fs);
imagesc(T,F,10*log10(P)); % converting to dB as usual
set(gca,'FontSize',20);
axis xy; xlabel('time (s)'); ylabel('Frequency (Hz)');  
 
hold on;
 
lfp_minmax = 25; lfp_cent = 125; % range and mean of LFP plotting
tvec0 = cscR.tvec - cscR.tvec(1); % align LFP with spectrogram
data = rescale(cscR.data,-lfp_minmax,lfp_minmax); data = data+lfp_cent;
 
lfp_h = plot(tvec0,data,'k');
xlim([tvec0(1) tvec0(end)]);
