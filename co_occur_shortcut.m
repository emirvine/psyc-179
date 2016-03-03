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


%% Load track spikes .mat
load('spike_pos_R063d3.mat');

for neuron = 1:length(spike_pos.u)
    spike_pos.u{neuron} = spike_pos.u{neuron}';
    spike_pos.shortcut{neuron} = spike_pos.shortcut{neuron}';
    spike_pos.novel{neuron} = spike_pos.novel{neuron}';
    spike_pos.other{neuron} = spike_pos.other{neuron}';
end

%% Estimating place fields
linspeed = getLinSpd([], pos);

cfg = [];
cfg.method = 'raw';
cfg.operation = '>';
cfg.threshold = 3.5;
iv_fast = TSDtoIV(cfg, linspeed);

pos_fast = restrict(pos, iv_fast); 
spikes_fast = restrict(spikes, iv_fast);

% plot(getd(pos_fast, 'x'), getd(pos_fast, 'y'), 'b')


%% Extracting trials
novels = iv;
novels.tstart = ExpKeys.novel_start;
novels.tend = ExpKeys.novel_stop;

shortcuts = iv;
shortcuts.tstart = ExpKeys.shortcut_start;
shortcuts.tend = ExpKeys.shortcut_stop;

us = iv;
us.tstart = ExpKeys.u_start;
us.tend = ExpKeys.u_stop;

trials = struct();
trials.novel = novels;
trials.shortcut = shortcuts;
trials.u = us;

trajectories = {'novel', 'shortcut', 'u'};
for side = 1:length(trajectories)
    pos_trial.(trajectories{side}) = restrict(pos_fast, trials.(trajectories{side}));
    spikes_trial.(trajectories{side}) = restrict(spikes_fast, trials.(trajectories{side}));
end

linear.novel = linear_shortcut(pos_fast, ExpKeys, 'novel1', 'novel2');
linear.shortcut = linear_shortcut(pos_fast, ExpKeys, 'shortcut1', 'shortcut2');
linear.u = linear_shortcut(pos_fast, ExpKeys, 'feeder1', 'feeder2');

for side = 1:length(trajectories)
    cfg = [];
    cfg.binSize = 1;
    tc.(trajectories{side}) = MakeTC(cfg, spikes_trial.(trajectories{side}), linear.(trajectories{side}));
    fields.(trajectories{side}) = tc.(trajectories{side}).field_template_idx;
end

%% check point
figure(1); clf; hold on;
for tuning = 1:size(tc.shortcut.tc, 1)
    plot(tc.shortcut.tc(tuning,:), 'b', 'LineWidth', 2);
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


%% Make a Q matrix
for side = 1:length(trajectories)
    cfg = [];
    cfg.win = 0.1;
    q.(trajectories{side}) = MakeQfromS2(cfg, spikes_side.(trajectories{side}), swr_evt);
end


%% Co-activation probabilities
for side = 1:length(trajectories)
    cfg = [];
    cfg.nShuffle = 1000;
    cfg.useMask = 1;
    cfg.outputFormat = 'vectorU';
    cooccur.(trajectories{side}) = CoOccurQ2(cfg, q.(trajectories{side}));
end


%% Plotting
colors = flipud(linspecer(3));
location = [1, 2, 3];
xticklabel = {'Novel', 'Shortcut', 'U'};

p0_data = [nanmean(cooccur.novel.p0), nanmean(cooccur.shortcut.p0), nanmean(cooccur.u.p0)];

figure(1);
for side = 1:length(p0_data)
    bar(location(side), p0_data(side), 'FaceColor', colors(side, :), 'EdgeColor', 'none');
    hold on;
end
set(gca, 'XLim', [location(1)-1, location(2)+1], 'XTick', location, 'XTickLabel', xticklabel);
xlabel('Field location');
ylabel({'Proportion of'; 'SWRs active'});
title('Activation probability (p0)');

p3_data = [nanmean(cooccur.novel.p3), nanmean(cooccur.shortcut.p3), nanmean(cooccur.u.p3)];

figure(2);
for side = 1:length(p3_data)
    bar(location(side), p3_data(side), 'FaceColor', colors(side, :), 'EdgeColor', 'none');
    hold on;
end
set(gca, 'XLim', [location(1)-1, location(2)+1], 'XTick', location, 'XTickLabel', xticklabel);
xlabel('Field location');
ylabel({'Cell pair'; 'joint probability'});
title('Observed coactivity (p3)');

p4_data = [nanmean(cooccur.novel.p4), nanmean(cooccur.shortcut.p4), nanmean(cooccur.u.p4)];

figure(3);
for side = 1:length(p4_data)
    bar(location(side), p4_data(side), 'FaceColor', colors(side, :), 'EdgeColor', 'none');
    hold on;
end
set(gca, 'XLim', [location(1)-1, location(2)+1], 'XTick', location, 'XTickLabel', xticklabel);
xlabel('Field location');
ylabel({'SWR coactivation'; 'z-scored'});
title('Coactivation above chance levels (p4)');

%% All plots together
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

set(fig, 'XLim', [location(1)-1, location(3)+1], 'XTick', location, 'XTickLabel', xticklabel);
set(fig, 'PlotBoxAspectRatio', [1, 1, 1]);
maximize