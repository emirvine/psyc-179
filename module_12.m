signal = randn(100, 4);
signal(:, 4) = sum(signal, 2);
[cor, p_val] = corrcoef(signal);
imagesc(cor);

%% Correlation matrix on vStr data
cd('C:\Users\Emily\Desktop\R016-2012-10-03_promoted');
cfg = [];
cfg.fc = {'R016-2012-10-03-CSC04a.ncs'};
csc = LoadCSC(cfg);

csc = restrict(csc, csc.cfg.ExpKeys.TimeOnTrack(2), csc.cfg.ExpKeys.TimeOffTrack(2));
fs = csc.cfg.hdr{1}.SamplingFrequency;

[s, f, t, p] = spectrogram(csc.data, hanning(512), 384, 1:0.25:200, fs);
[cor, p_val] = corrcoef(10*log10(p'));

imagesc(f, f, cor);
caxis([-0.1, 0.5]);
axis xy;
colorbar;
grid on;
set(gca, 'XLim', [0, 150], 'FontSize', 14, 'XTick', 0:10:150, 'YTick', 0:10:150);

%% Autocorrelation for continuous signals
noise = randn(1000, 1);
[cor, p_val] = corrcoef(noise(1:end-1), noise(2:end));

[autocorr, lags] = xcorr(noise, 100, 'coeff');
plot(lags, autocorr, 'LineWidth', 2); grid on;
set(gca, 'FontSize', 18);
xlabel('Time lag (samples)');
ylabel('Correlation ({\itr})');

%% Correlated signal
fs = 500;
dt = 1 ./ fs;
time = 0:dt:2-dt;

noise = randn(size(time));
signal = sin(2*pi*10*time) + noise;

subplot(211);
plot(time, signal); grid on;
set(gca, 'FontSize', 18);
title('Signal');

subplot(212);
[autocorr, lags] = xcorr(signal, 100, 'coeff');
lags = lags .* dt;
plot(lags, autocorr); grid on;
set(gca, 'FontSize', 18);
xlabel('Time lag (s)');
ylabel('Correlation ({\itr})');
title('Autocorrelation');

%% Application to real data
fs = 500;
csc.data = decimate(csc.data, 4);
csc.tvec = downsample(csc.tvec, 4);
[s, f, t, p] = spectrogram(csc.data, hanning(125), 115, 1:200, fs);

%%
freq_idx = find(f > 70 & f < 100); % high gamma
freq_idx2 = find(f > 50 & f < 65); % low gamma

power = mean(p(freq_idx, :));
power2 = mean(p(freq_idx2, :));

[autocorr, lags] = xcorr(power2 - mean(power2), 50, 'coeff');
lags = lags .* mean(diff(t));

figure;
plot(lags, autocorr);
set(gca, 'FontSize', 18);
xlabel('Time lag (s)');
ylabel('Correlation ({\itr})');
title('Low-gamma autocorrelation');

%%
freq_idx = find(f > 70 & f < 100); % high gamma
freq_idx2 = find(f > 50 & f < 65); % low gamma

power = mean(p(freq_idx, :));
power2 = mean(p(freq_idx2, :));

[autocorr, lags] = xcorr(power - mean(power), power2 - mean(power2), 50, 'coeff');
lags = lags .* mean(diff(t));

figure;
plot(lags, autocorr);
set(gca, 'FontSize', 18);
xlabel('Time lag (s)');
ylabel('Correlation ({\itr})');
title({'Crosscorrelation', '(high-gamma relative to low-gamma at time=0)'});

%% Extracting signal phases
fs = 500;
dt = 1 ./ fs;
time = 0:dt:1-dt;

noise = randn(size(time));
signal = sin(2*pi*4*time) + 0.1*noise;

hilber = hilbert(signal);
phi = angle(hilber);

[ax_h, h1, h2] = plotyy(time, signal, time, phi);
set(ax_h(2), 'YColor', [1,0,0], 'YLim', [-pi, pi], 'Ytick', -pi:pi/2:pi, ...
    'YTickLabel', {'-p', '-p/2', '0', 'p/2', 'p'}, 'fontname', 'symbol', 'XTick', []);
set(h2, 'Color', [1,0,0]);
set(get(ax_h(2), 'YLabel'), 'String', 'Phase (radians)');
box off;

%% Relating signal phases to amplitude
fs = 500;
dt = 1 ./ fs;
time = 0:dt:1-dt;

freq1 = 8;
signal1 = sin(2*pi*freq1*time) + 0.05*randn(size(time));
signal1_phi = angle(hilbert(signal1));
dphi = exp(-abs(signal1_phi)/1.5);

freq2 = 80;
signal2 = 0.3 .* sin(2*pi*freq2*time) .* dphi;

signal = signal1 + signal2;
plot(time, signal);
xlim([0,0.5]);

%%
phi_edges = -pi:pi/8:pi;
pow_bin = average_x_by_y_bin(dphi, signal1_phi, phi_edges);

pow_bin(end-1) = pow_bin(end-1) + pow_bin(end);
pow_bin = pow_bin(1:end-1);

phi_centers = phi_edges(1:end-1) + pi/16;
plot(phi_centers, pow_bin);

%% Real data
addpath('C:\Users\Emily\Code\fieldtrip');
ft_defaults;
addpath('C:\Users\Emily\Code\neuraldata-w16\shared');
rmpath('C:\Users\Emily\Code\fieldtrip\external\signal');

cd('C:\Users\Emily\Desktop\R016-2012-10-03_promoted');

%%
twin = [-1 3]; % time window of interest, relative to reward delivery
f_lg = [50 65]; % low-gamma freq band, to get power
f_hg = [70 85]; % high-gamma freq band, to get power
f_x = [8 12]; % low-frequency band, to get phase
 
LoadExpKeys;
fc = ExpKeys.goodGamma(1);
data = ft_read_neuralynx_interp(fc);
%%
data2 = data; % temporary data variable
 
% low-gamma
data2f = ft_filterLFP(data2,f_lg,'ford',4); % a FieldTrip wrapper for the filtfilt() function
data.trial{1}(2,:) = data2f.trial{1};
data.label{2} = 'lg_pow';
data.trial{1}(2,:) = abs(hilbert(data.trial{1}(2,:))).^2;
 
% high-gamma
data2f = ft_filterLFP(data2,f_hg,'ford',4);
data.trial{1}(3,:) = data2f.trial{1};
data.label{3} = 'hg_pow';
data.trial{1}(3,:) = abs(hilbert(data.trial{1}(3,:))).^2;

%%
d = fdesign.bandpass('N,F3dB1,F3dB2',4,f_x(1),f_x(2),2e3);
Hd = design(d,'butter');
 
data2f = ft_filterLFP(data2,f_x,'B',Hd.sosMatrix,'A',Hd.scaleValues);
data.trial{1}(4,:) = data2f.trial{1};
data.label{4} = 'x_phase';
data.trial{1}(4,:) = angle(hilbert(data.trial{1}(4,:)));

%%
all_cues = {'c1','c3','c5','lo','hi'}; 

cfg = [];
cfg.trialdef.cue = all_cues;
cfg.trialdef.block = 'both';
cfg.trialdef.eventtype = 'nosepoke';
cfg.trialdef.location = 'both';
cfg.trialfun = 'ft_trialfun_lineartracktone2';
cfg.trialdef.hdr = data.hdr;
cfg.trialdef.pre = 2; cfg.trialdef.post = 4;
 
trl = ft_trialfun_lineartracktone2(cfg); cfg.trl = trl;
this_data = ft_redefinetrial(cfg,data);

%%
dpi = pi/18; % binsize for binning phase
pbin_edges = -pi:dpi:pi; pbin_centers = pbin_edges(1:end-1)+dpi/2;
maxlag = 2000; % samples for xcorr
nShuf = 100; % number of shuffles to perform
 
% restrict data
cfg = []; cfg.latency = twin;
temp = ft_selectdata(cfg,this_data);
trial_len = length(temp.trial{1}(1,:));

%%
% first get observed PAC value
trial_mat = cell2mat(temp.trial); % trials concatenated
 
lg_idx = strmatch('lg_pow',temp.label);
obs_lg_pwr = trial_mat(lg_idx,:); % note, ft_selectdata switches labels!!
 
hg_idx = strmatch('hg_pow',temp.label);
obs_hg_pwr = trial_mat(hg_idx,:);
 
xf_idx = strmatch('x_phase',temp.label);
obs_xf_phase = trial_mat(xf_idx,:);
 
% explain how to quantify PAC.. mean vector length
obs_pac_lg = abs(mean(obs_lg_pwr.*exp(1i*obs_xf_phase))); % mean vector length
obs_pac_hg = abs(mean(obs_hg_pwr.*exp(1i*obs_xf_phase)));
 
% get observed phase-frequency relationships
obs_xf_lg = averageXbyYbin(obs_lg_pwr,obs_xf_phase,pbin_edges);
obs_xf_hg = averageXbyYbin(obs_hg_pwr,obs_xf_phase,pbin_edges);

%%
shuf_pac_lg = zeros(nShuf,1); shuf_pac_hg = zeros(nShuf,1);
shuf_xf_lg = zeros(nShuf,length(pbin_edges)-1); shuf_xf_hg = zeros(nShuf,length(pbin_edges)-1);
 
for iShuf = 1:nShuf
    shuf_data = temp;
 
    % random amount to shift data
    perm_shift1 = ceil(length(shuf_data.trial{1}(lg_idx,:)-1)*rand(1));
    perm_shift2 = ceil(length(shuf_data.trial{1}(lg_idx,:)-1)*rand(1));
 
    % for each trial, shift the data
    for iT = 1:length(shuf_data.trial)
 
        piece1 = shuf_data.trial{iT}(lg_idx,perm_shift1+1:end); piece2 = shuf_data.trial{iT}(lg_idx,1:perm_shift1);
        shuf_data.trial{iT}(lg_idx,1:trial_len-perm_shift1) = piece1;
        shuf_data.trial{iT}(lg_idx,end-perm_shift1+1:end) = piece2;
 
        piece1 = shuf_data.trial{iT}(hg_idx,perm_shift2+1:end); piece2 = shuf_data.trial{iT}(hg_idx,1:perm_shift2);
        shuf_data.trial{iT}(hg_idx,1:trial_len-perm_shift2) = piece1;
        shuf_data.trial{iT}(hg_idx,end-perm_shift2+1:end) = piece2;
 
 
    end
 
    % get PAC for this shuffle
    trial_mat = cell2mat(shuf_data.trial); % trials concatenated
 
    shuf_lg_pwr = trial_mat(lg_idx,:);
    shuf_hg_pwr = trial_mat(hg_idx,:);
 
    shuf_xf_lg(iShuf,:) = averageXbyYbin(shuf_lg_pwr,obs_xf_phase,pbin_edges);
    shuf_xf_hg(iShuf,:) = averageXbyYbin(shuf_hg_pwr,obs_xf_phase,pbin_edges);
 
    shuf_pac_lg(iShuf) = abs(mean(shuf_lg_pwr.*exp(1i*obs_xf_phase)));
    shuf_pac_hg(iShuf) = abs(mean(shuf_hg_pwr.*exp(1i*obs_xf_phase)));
 
end

%%
pacz_lg = (obs_pac_lg-mean(shuf_pac_lg))/std(shuf_pac_lg);
pacz_hg = (obs_pac_hg-mean(shuf_pac_hg))/std(shuf_pac_hg);

