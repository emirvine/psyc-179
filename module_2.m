cd('C:\Users\Emily\Desktop\R016-2012-10-08_promoted'); % same session as Module 1

cfg_csc = [];
cfg_csc.fc = {'R016-2012-10-08-CSC02d.ncs'};
csc = LoadCSC(cfg_csc);

cfg_evt = [];
evt = LoadEvents(cfg_evt);

plot(getd(evt,'1 pellet cue'), 0, 'k.')
%plot(csc.tvec,csc.data, 'b.')

spike = LoadSpikes([]);
 
% load video data (make sure the VT1.zip file is unzipped first and now present in MATLAB's working folder!)
[Timestamps, x, y, Angles, Targets, Points, Header] = Nlx2MatVT('VT1.nvt', [1 1 1 1 1 1], 1, 1, [] );

cfg_pos = [];
cfg_pos.fn = 'VT1.nvt';
cfg_pos.tsflag = 'sec';
pos_tsd = LoadPos(cfg_pos);

%%
cd('C:\Users\Emily\Desktop\R042-2013-08-18');
LoadMetadata;
 
pos = LoadPos([]);
 
left_pos = restrict(pos,metadata.taskvars.trial_iv_L); % left trials only
plot(getd(left_pos,'x'),getd(left_pos,'y'),'.'); % looks like right trials! camera reverses image, can fix with set(gca,'YDir','reverse')
 
% example 2: interplay between tsd and iv data
LoadExpKeys;
 
please = []; please.fc = ExpKeys.goodSWR(1); % local field potential with good "sharp wave-ripple" events
lfp = LoadCSC(please); % aacarey is Canadian and asks nicely; cfg name is just arbitrary
 
% detect possible artifacts
cfg = [];
cfg.method = 'zscore'; % first normalize the data
cfg.threshold = -8;
cfg.minlen = 0; % no minimum length on events to detect
cfg.dcn = '<'; % detect intervals with z-score lower than threshold
artifact_iv = TSDtoIV(cfg,lfp); % creates iv with start and end times of possible artifacts
 
% plot detected intervals
cfg = []; cfg.display = 'tsd'; % also try 'iv' mode!
PlotTSDfromIV(lfp,artifact_iv,lfp)

plot(getd(pos_tsd,'x'), getd(pos_tsd,'y'));
plot(pos_tsd.tvec, getd(pos_tsd,'x'));

clear all;
fname = 'R016-2012-10-08-CSC02c.ncs';
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);

fn = FindFile('*Events.nev');
[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);

cd('C:\Users\Emily\Desktop\R042-2013-08-18'); %Matt's fav session
LoadMetadata;
 
pos = LoadPos([]);
 
left_pos = restrict(pos,metadata.taskvars.trial_iv_L); % left trials only
plot(getd(left_pos,'x'),getd(left_pos,'y'),'.'); % looks like right trials! camera reverses image, can fix with set(gca,'YDir','reverse')
 
% example 2: interplay between tsd and iv data
LoadExpKeys;
 
please = []; please.fc = ExpKeys.goodSWR(1); % local field potential with good "sharp wave-ripple" events
lfp = LoadCSC(please); % aacarey is Canadian and asks nicely; cfg name is just arbitrary
 
% detect possible artifacts
cfg = [];
cfg.method = 'zscore'; % first normalize the data
cfg.threshold = -8;
cfg.minlen = 0; % no minimum length on events to detect
cfg.dcn = '<'; % detect intervals with z-score lower than threshold
artifact_iv = TSDtoIV(cfg,lfp); % creates iv with start and end times of possible artifacts
 
% plot detected intervals
cfg = []; cfg.display = 'tsd'; % also try 'iv' mode!
PlotTSDfromIV(lfp,artifact_iv,lfp)