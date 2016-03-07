cd('C:\Users\Emily\Desktop\R063-2015-03-23_recording');
LoadExpKeys;


%% Get SWR events
cfg = [];
cfg.fc = ExpKeys.good_lfp;
csc = LoadCSC(cfg);

csc = restrict(csc, ExpKeys.pauseA(1), ExpKeys.pauseA(2));

cfg = [];
cfg.load_questionable_cells = 1;
spikes = LoadSpikes(cfg);

load('R063-2015-03-23-emi-vt.mat');

pos_tsd.data(1,:) = pos_tsd.data(1,:)./ExpKeys.pxl_to_cm(1); 
pos_tsd.data(2,:) = pos_tsd.data(2,:)./ExpKeys.pxl_to_cm(2);

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

clearvars -except pos_tsd swr_evt ExpKeys spikes


%% Load track spikes .mat
% load('spike_pos_14_R063d3.mat');
load('spike_pos_R063d4.mat');
% load('saved_struct.mat');

trajectories = {'novel', 'shortcut', 'u'};
for side = 1:length(trajectories)
    spikes_track.(trajectories{side}) = ts;
    spikes_track.(trajectories{side}).usr = spikes.usr;
    spikes_track.(trajectories{side}).t = spike_pos.(trajectories{side});
end

%% Estimating place fields
linspeed = getLinSpd([], pos_tsd);

cfg = [];
cfg.method = 'raw';
cfg.dcn = '>';
cfg.threshold = 3.5;
iv_fast = TSDtoIV(cfg, linspeed);

pos_fast = restrict(pos_tsd, iv_fast); 

for side = 1:length(trajectories)
    spikes_fast.(trajectories{side}) = ts;
    spikes_fast.(trajectories{side}).usr = spikes.usr;
    spikes_fast.(trajectories{side}) = restrict(spikes_track.(trajectories{side}), iv_fast);
end


% plot(getd(pos_fast, 'x'), getd(pos_fast, 'y'), 'b')


%% Extracting trials

linear.novel = linear_shortcut(pos_fast, ExpKeys, 'novel1', 'novel2');
linear.shortcut = linear_shortcut(pos_fast, ExpKeys, 'shortcut1', 'shortcut2');
linear.u = linear_shortcut(pos_fast, ExpKeys, 'feeder1', 'feeder2');

for side = 1:length(trajectories)
    cfg = [];
    cfg.binSize = 1;
    tc.(trajectories{side}) = MakeTC(cfg, spikes_fast.(trajectories{side}), linear.(trajectories{side}));
    fields.(trajectories{side}) = tc.(trajectories{side}).field_template_idx;
end

%% check point
figure(1); clf; hold on;
tc_plot = tc.novel.tc;
for tuning = 1:size(tc_plot, 1)
    plot(tc_plot(tuning,:), 'b', 'LineWidth', 2);
end
set(gca,'FontSize',14); xlabel('Linearized position (cm)'); ylabel('Firing rate (Hz)');


%% Categorizing place cells
[~, remove.u, remove.shortcut] = intersect(fields.u, fields.shortcut);
u_shortcut = {'u', 'shortcut'};
for side = 1:length(u_shortcut)
    fields.(u_shortcut{side})(remove.(u_shortcut{side})) = [];
    fields.(u_shortcut{side}) = unique(fields.(u_shortcut{side}));
    spikes_side.(u_shortcut{side}) = SelectTS([], spikes, fields.(u_shortcut{side}));
end

[~, remove.u, remove.novel] = intersect(fields.u, fields.novel);
u_novel = {'u', 'novel'};
for side = 1:length(u_novel)
    fields.(u_novel{side})(remove.(u_novel{side})) = [];
    fields.(u_novel{side}) = unique(fields.(u_novel{side}));
    spikes_side.(u_novel{side}) = SelectTS([], spikes, fields.(u_novel{side}));
end

[~, remove.shortcut, remove.novel] = intersect(fields.shortcut, fields.novel);
shortcut_novel = {'shortcut', 'novel'};
for side = 1:length(shortcut_novel)
    fields.(shortcut_novel{side})(remove.(shortcut_novel{side})) = [];
    fields.(shortcut_novel{side}) = unique(fields.(shortcut_novel{side}));
    spikes_side.(shortcut_novel{side}) = SelectTS([], spikes, fields.(shortcut_novel{side}));
end


%% Make a Q matrix
for side = 1:length(trajectories)
    cfg = [];
    cfg.win = 0.1;
    q.(trajectories{side}) = MakeQfromS2(cfg, spikes_side.(trajectories{side}), swr_evt);
end


%% Co-activation probabilities
for side = 1:length(trajectories)
    cfg = [];
    cfg.nShuffle = 10000;
    cfg.useMask = 1;
    cfg.outputFormat = 'vectorU';
    cooccur.(trajectories{side}) = CoOccurQ2(cfg, q.(trajectories{side}));
end


%% Plotting
p_list = {'p0', 'p3', 'p4'};
titles = {'Activation probability (p0)', 'Observed coactivity (p3)', 'Coactivation above chance levels (p4)'};
ylabels = {{'Proportion of'; 'SWRs active'}, {'Cell pair'; 'joint probability'}, {'SWR coactivation'; 'z-scored'}};
trajectories = {'novel', 'shortcut', 'u'};
colors = flipud(linspecer(3));
location = [0.0, 0.75, 1.5];
xticklabel = {'Novel', 'Shortcut', 'U'};
 
figure;
for prob = length(p_list):-1:1
    p_data(prob, 1:length(trajectories)) = [nanmean(cooccur.novel.(p_list{prob})), nanmean(cooccur.shortcut.(p_list{prob})), nanmean(cooccur.u.(p_list{prob}))];
    fig(prob) = subplot(1, 3, prob);
    for side = 1:length(trajectories)
        bar(location(side), p_data(prob, side), 0.6, 'FaceColor', colors(side,:), 'Edgecolor', 'none');
        hold on;
        
        title(titles{prob});
        ylabel(ylabels{prob});
        xlabel('Field location');
    end
end

set(fig, 'XLim', [location(1)-1, location(3)+1], 'XTick', location, 'XTickLabel', xticklabel, 'FontSize', 12);
set(fig, 'PlotBoxAspectRatio', [1, 1, 1]);
maximize