cd('C:\Users\Emily\Desktop\R016-2012-10-03_promoted');
cfg_csc.fc = {'R016-2012-10-03-CSC04a.ncs'};
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

%% Getting events
cfg_evt = [];
 
cfg_evt.eventList = {'Feeder 0 nosepoke','Feeder 1 nosepoke', ...
    '1 pellet cue','3 pellet cue','5 pellet cue', ...
    '1 pellet dispensed','3 pellet dispensed','5 pellet dispensed'};
 
cfg_evt.eventLabel = {'n0','n1', 'c1','c3','c5', 'd1','d3','d5'};
 
evt = LoadEvents(cfg_evt);

%% Loads neuralynx data into fieldtrip
cd('C:\Users\Emily\Desktop\R016-2012-10-03_promoted');
fc = {'R016-2012-10-03-CSC04a.ncs'};
data = ft_read_neuralynx_interp(fc);

%% trial data in fieldtrip
cfg = [];
cfg.t = cat(2,getd(evt,'n0'),getd(evt,'n1'));
cfg.mode = 'nlx';
cfg.hdr = data.hdr;
cfg.twin = [-1 4];
 
trial = ft_maketrl(cfg);
 
cfg = [];
cfg.trl = trial;
data_trial = ft_redefinetrial(cfg,data); 

%% event-triggered spectrogram
cfg = []; % start with empty cfg
cfg.output = 'pow';
cfg.channel = 'R016-2012-10-03-CSC04a';
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 10:2:100; % frequencies of interest
cfg.t_ftimwin = ones(size(cfg.foi)).*0.5;  % window size: fixed at 0.5s
cfg.toi = -0.5:0.05:3.5; % times of interest
 
TFR = ft_freqanalysis(cfg, data_trl);
 
figure
cfg = []; cfg.channel = 'R016-2012-10-03-CSC04a';
ft_singleplotTFR(cfg, TFR);

%% Baseline correction
figure
cfg = [];
cfg.baseline = [-2 0];
cfg.baselinetype = 'relative';
cfg.channel = 'R016-2012-10-03-CSC04a';
ft_singleplotTFR(cfg, TFR);

%% Statistical tests
cfg = [];
cfg.output = 'pow';
cfg.channel = 'R016-2012-10-03-CSC04a';
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 1:1:100;
cfg.keeptrials = 'yes'; % need this for stats later
cfg.t_ftimwin = 20./cfg.foi;  % 20 cycles per time window
 
cfg.toi = -1:0.05:0; % pre-nosepoke baseline
TFR_pre = ft_freqanalysis(cfg, data_trl);
 
cfg.toi = 0:0.05:1; % post-nosepoke
TFR_post = ft_freqanalysis(cfg, data_trl);
 
TFR_pre.time = TFR_post.time; % time should be identical for comparison

%% t-test
cfg = [];
cfg.channel = 'R016-2012-10-03-CSC04a';
cfg.latency = 'all';
cfg.trials = 'all';
cfg.frequency = 'all';
cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';
cfg.avgoverfreq = 'no';
cfg.parameter = 'powspctrm';
cfg.method = 'stats';
cfg.statistic = 'ttest2';
cfg.alpha = 0.05;
 
nTrials1 = size(TFR_pre.powspctrm,1); 
nTrials2 = size(TFR_post.powspctrm,1);
 
cfg.design = cat(2,ones(1,nTrials1),2*ones(1,nTrials2)); % two conditions
cfg.ivar = 1; % dim of design var with the independent variable (group)
 
stat = ft_freqstatistics(cfg,TFR_post,TFR_pre);
 
cfg.parameter = 'stat';
ft_singleplotTFR(cfg,stat); % plot the t-statistic

%% Multichannel data
cd('C:\Users\Emily\Desktop\R016-2012-10-03_promoted');
fc = FindFiles('*.ncs'); % get filenames of all LFPs recorded
data_all = ft_read_neuralynx_interp(fc); % load them all -- this will take a while
data_all.hdr.Fs = data_all.fsample; % for some reason this is missing from the header
 
% define layout for later plotting
cfg = [];
cfg.layout = 'ordered'; cfg.channel = data_all.label;
layout = ft_prepare_layout(cfg, data_all);

%%
cfg = [];
cfg.t = cat(2,getd(evt,'n0'),getd(evt,'n1'));
cfg.mode = 'nlx';
cfg.hdr = data_all.hdr;
cfg.twin = [-1 4];
 
trl = ft_maketrl(cfg);
 
cfg = [];
cfg.trl = trl;
data_trl = ft_redefinetrial(cfg,data_all); 

cfg = []; % start with empty cfg
cfg.output = 'pow';
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 1:100; % frequencies of interest
cfg.toi = -0.5:0.05:3.5; % times of interest
cfg.t_ftimwin = 20./cfg.foi;
 
TFR = ft_freqanalysis(cfg, data_trl);

%% Plotting
figure
cfg = [];
cfg.baseline = [-2 0];
cfg.baselinetype = 'relative';
cfg.layout = layout;
 
ft_multiplotTFR(cfg, TFR); % note this is now multiplot rather than singleplot