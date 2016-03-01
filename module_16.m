cd('C:\Users\Emily\Desktop\R064-2015-04-22');
LoadExpKeys;
cfg = [];
cfg.fc = ExpKeys.goodSWR(1);
csc = LoadCSC(cfg);

cfg = [];
cfg.load_questionable_cells = 1;
spikes = LoadSpikes(cfg);

cfg = [];
cfg.convFact = ExpKeys.convFact;
pos = LoadPos(cfg);

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