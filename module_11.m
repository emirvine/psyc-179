fs = 500;
dt = 1 ./ fs;
time_lim = [0, 2];
time = time_lim(1):dt:time_lim(2)-dt;

freq1 = 8;
data1 = sin(2*pi*freq1*time)+0.1*randn(size(time));

[auto, lag] = xcorr(data1, 100, 'coeff');
lag = lag .* (1/fs);
plot(lag, auto); grid on;

freq2 = 8;
data2 = sin(2*pi*freq2*time+pi/4) + 0.1*randn(size(time));

[cross, lag] = xcorr(data1, data2, 100, 'coeff');
lag = lag .* (1./fs);
plot(lag, cross); grid on;

figure;
subplot(221);
plot(time, data1, 'r', time, data2, 'b');
legend({'Signal 1', 'Signal 2'});
title('Raw signals');

[power_freq1, norm_freq] = pwelch(data1, hanning(250), 125, length(data1), fs);
[power_freq2, norm_freq] = pwelch(data2, hanning(250), 125, length(data1), fs);
subplot(222);
plot(norm_freq, abs(power_freq1),'r', norm_freq, abs(power_freq2), 'b');
xlim([0, 100]);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Cross-spectrum');

[power_both, norm_freq] = cpsd(data1, data2, hanning(250), 125, length(data1), fs);
subplot(223);
plot(norm_freq, abs(power_both));
xlim([0,100]);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Cross-spectrum');

[cross_freq, lag] = xcorr(data1, data2, 100, 'coeff');
lag = lag .*(1/fs);
subplot(224);
plot(lag, cross_freq); grid on;
xlabel('Time lag (s)');
ylabel('Correlation ({\itr})');
title('Xcorr')


%% Properties of the coherence measure
freq = 0.5;
freq2 = 8;
modulation = 4;
window = 250;

subplot(421);
signal2 = data2;
plot(time, signal2, time, data1);
title('Signal 1 = constant phase');

subplot(422);
signal3 = sin(2*pi*freq2*time + modulation .* sin(2*pi*freq*time - pi/2)) + 0.1*randn(size(time));
plot(time, signal3, time, data1);
title('Signal 2 - varying phase');

subplot(423);
[power2, norm_freq] = pwelch(signal2, hanning(window), window/2, length(data2),fs);
plot(norm_freq, abs(power2));
title('PSD signal2');

subplot(424);
[power3, norm_freq] = pwelch(signal3, hanning(window), window/2, length(data2),fs);
plot(norm_freq, abs(power3));
title('PSD signal3');

subplot(425);
[coherence1, norm_freq] = mscohere(data1, signal2, hanning(window), window/2, length(data1), fs);
plot(norm_freq, coherence1);
title('Coherence signal2');
xlabel('Frequency (Hz)');

subplot(426);
[coherence2, norm_freq] = mscohere(data1, signal3, hanning(window), window/2, length(data1), fs);
plot(norm_freq, coherence2);
title('Coherence signal3');
xlabel('Frequency (Hz)');

[auto1, lag] = xcorr(data1, signal2, 100, 'coeff');
lag = lag .* (1/fs);
subplot(427);
plot(lag, auto1); grid on;
xlabel('Time lag (s)');
ylabel('Correlation ({\itr})')
title('Xcorr signal2');

[auto2, lag] = xcorr(data1, signal3, 100, 'coeff');
lag = lag .* (1/fs);
subplot(428);
plot(lag, auto2); grid on;
xlabel('Time lag (s)');
ylabel('Correlation ({\itr})')
title('Xcorr signal2');


%% another example.

window = 50;
fs = 500;
dt = 1./fs;
time_lim = [0, 2];
time = time_lim(1):dt:time_lim(2)-dt;
freq1 = 40;
freq2 = 40;

mod1 = square(2*pi*4*time, 20);
mod1(mod1 < 0) = 0;

mod2 = square(2*pi*4*time+pi, 20);
mod2(mod2 < 0) = 0;

data1 = sin(2*pi*freq1*time);
data1 = data1 .* mod1 + 0.01*randn(size(time));

data2 = sin(2*pi*freq2*time);
data2 = data2 .* mod2 + 0.01*randn(size(time));

subplot(221);
plot(time, data1, 'r', time, data2, 'b');
legend({'Signal 1', 'Signal 2'});
title('Raw signals');

[power1, norm_freq] = pwelch(data1, hanning(window), window/2, length(data1), fs);
[power2, norm_freq] = pwelch(data2, hanning(window), window/2, length(data2), fs);
subplot(222);
plot(norm_freq, abs(power1), 'r', norm_freq, abs(power2), 'b');
title('PSD');

subplot(223);
[coherence, norm_freq] = mscohere(data1, data2, hanning(window), window/2, length(data1), fs);
plot(norm_freq);
title('Coherence');
xlabel('Frequency (Hz)');

[cross, lag] = xcorr(data1, data2, 100, 'coeff');
lag = lag .* (1/fs);
subplot(224);
plot(lag, cross); grid on;
xlabel('Time lag (s)');
ylabel('Correlation ({\itr})');
title('Xcorr');


%% real data application
cd('C:\Users\Emily\Desktop\R016-2012-10-03_promoted');
LoadExpKeys;

cfg = [];
cfg.fc = cat(2, ExpKeys.goodGamma(1:2), ExpKeys.goodTheta(1));
cfg.label = {'vStr1', 'vStr2', 'HC'};
csc = LoadCSC(cfg);

csc = restrict(csc, ExpKeys.TimeOnTrack(1), ExpKeys.TimeOffTrack(2));

%% cont'd
fs = csc.cfg.hdr{1}.SamplingFrequency;
window = 2048;

num_samples = length(csc.label);

for sample = 1:num_samples
    [power{sample}, norm_freq{sample}] = pwelch(getd(csc, csc.label{sample}), hanning(window), window/2, 2*window, fs);
    for sample2 = sample+1:num_samples
        [cohere{sample, sample2}, norm_freqc{sample}] = mscohere(getd(csc, csc.label{sample}), getd(csc, csc.label{sample2}), hanning(window), window/2, 2*window, fs);
    end
end

subplot(121);
colors = 'kgm';
for sample = 1:num_samples
    h(sample) = plot(norm_freq{sample}, 10*log10(power{sample}), colors(sample), 'LineWidth', 2);
    hold on;
end
set(gca, 'XLim', [0, 150], 'XTick', 0:25:150, 'FontSize', 12); grid on;
legend(h, csc.label, 'Location', 'Northeast');
legend boxoff;
xlabel('Frequenct (Hz)');
ylabel('Power (dB)');

subplot(122); clear h;
h(1) = plot(norm_freqc{1}, cohere{1,2}, 'LineWidth', 2); 
hold on;
h(2) = plot(norm_freqc{1}, cohere{1,3}, 'LineWidth', 2);
set(gca,'XLim', [0,150], 'XTick', 0:25:150, 'FontSize', 12);
grid on;
legend(h, {'cStr1-vStr2', 'vStr1-HC'}, 'Location', 'Northeast');
legend boxoff;
xlabel('Frequency (Hz)');
ylabel('Coherence');

%% Directly from Module 11
cd('C:\Users\Emily\Desktop\R016-2012-10-03_promoted'); 
LoadExpKeys;
 
cfg.fc = cat(2,ExpKeys.goodGamma(1:2),ExpKeys.goodTheta(1));
data = ft_read_neuralynx_interp(cfg.fc);
data.label = {'vStr1','vStr2','HC'};
%% trialify (directly from Module 11)
data.hdr.Fs = data.fsample;
 
cfg = [];
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = 2.5; cfg.trialdef.post = 5; % define time window of interest
 
cfg.trialdef.eventtype = 'nosepoke'; % could be 'nosepoke', 'reward', 'cue'; this and what follows are all task-specific
cfg.trialdef.location = 'both'; % could be 'left', 'right', 'both'
cfg.trialdef.block = 'both'; % could be 'value', 'risk', 'both'
cfg.trialdef.cue = {'c1','c3','c5'}; % cell array with choice of elements {'c1','c3','c5','lo','hi'} (1, 3, 5 pellets; low and high risk)
 
[trl, event] = ft_trialfun_lineartracktone2(cfg);
cfg.trl = trl;
 
data_trl = ft_redefinetrial(cfg,data);

cfg              = [];
cfg.output       = 'powandcsd';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:100; % frequencies to use
cfg.t_ftimwin    = 20./cfg.foi;  % frequency-dependent, 20 cycles per time window
cfg.keeptrials   = 'yes';
cfg.channel      = {'vStr1', 'vStr2', 'HC1'};
cfg.channelcmb   = {'vStr2', 'HC1'; 'vStr2', 'vStr1'}; % channel pairs to compute csd for

cfg.toi          = -2:0.05:0; % pre-nosepoke baseline (time 0 is time of nosepoke)

TFR_pre = ft_freqanalysis(cfg, data_trl);



cfg            = [];
cfg.method     = 'coh'; % compute coherence; other measures of connectivity are also available
fd             = ft_connectivityanalysis(cfg,TFR_pre);


figure;
cols = 'rgb';
for iCmb = 1:size(fd.labelcmb,1)
    lbl{iCmb} = cat(2,fd.labelcmb{iCmb,1},'-',fd.labelcmb{iCmb,2});
 
    temp = nanmean(sq(fd.cohspctrm(iCmb,:,:)),2);
    h(iCmb) = plot(fd.freq,temp,cols(iCmb));
    hold on;
end
legend(h,lbl);


%% time-frequency coherence analysis (directly from Mondule 11)
iC = 1; % which signal pair to plot
lbl = [fd.labelcmb{1,:}]; % get the label of this pair
imagesc(fd.time,fd.freq,sq(fd.cohspctrm(iC,:,:))); axis xy; colorbar
xlabel('time (s)'); ylabel('Frequency (Hz)'); title(lbl);