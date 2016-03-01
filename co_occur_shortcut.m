cd('C:\Users\Emily\Desktop\R063-2015-03-22_recording');
LoadExpKeys;

%% Get SWR events
cfg = [];
cfg.fc = {'R063-2015-03-22-CSC14d.ncs'};
csc = LoadCSC(cfg);

cfg = [];
cfg.load_questionable_cells = 1;
spikes = LoadSpikes(cfg);

load('R063-2015-03-22-emi-vt.mat');

pos_tsd.data(1,:) = pos_tsd.data(1,:)./ExpKeys.pxl_to_cm(1); 
pos_tsd.data(2,:) = pos_tsd.data(2,:)./ExpKeys.pxl_to_cm(2);

pos = pos_tsd;

cfg = [];
cfg.type = 'fdesign';
cfg.f = [140, 250];
csc_filtered = FilterLFP(cfg, csc);

envelope = abs(csc_filtered.data);

kernel = gausskernel(60, 20);
envelope = conv(envelope, kernel, 'same');

swr = tsd(csc.tvec, envelope);

cfg = [];
cfg.method = 'zscore';
cfg.threshold = 3;
cfg.minlen = 0.02;
cfg.merge_thr = 0;
swr_evt = TSDtoIV(cfg, swr);

clearvars -except spikes pos swr_evt ExpKeys


%% Estimating place fields
linspeed = getLinSpd([], pos);

cfg = [];
cfg.method = 'raw';
cfg.operation = '>';
cfg.threshold = 3.5;
iv_fast = TSDtoIV(cfg, linspeed);

pos_fast = restrict(pos, iv_fast); 
spikes_fast = restrict(spikes, iv_fast);

plot(getd(pos_fast, 'x'), getd(pos_fast, 'y'), 'b')


%% Extracting trials
pos_phase3 = restrict(pos_fast, ExpKeys.phase3(1), ExpKeys.phase3(2));
[u_times, shortcut_times, novel_times, ~, ~, ~, ~, ~, ~] = shortcut_trial_idx(pos_phase3, ExpKeys);


