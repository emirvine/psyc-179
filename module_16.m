cd('C:\Users\Emily\Desktop\R064-2015-04-22');

%%
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
LoadMetadata;

[trials.left, trials.right] = GetMatchedTrials([], metadata, ExpKeys);

sides = {'left', 'right'};
shorthand = {'L', 'R'};
for side = 1:length(sides)
    pos_trial.(sides{side}) = restrict(pos_fast, trials.(sides{side}));
    spikes_trial.(sides{side}) = restrict(spikes_fast, trials.(sides{side}));
    
    cfg = [];
    cfg.binsize = 3;
    cfg.run_dist = ExpKeys.pathlength;
    coord.(sides{side}) = StandardizeCoord(cfg, metadata.coord.(['coord', shorthand{side},'_cm']));
    
    cfg = [];
    cfg.Coord = coord.(sides{side});
    linpos.(sides{side}) = LinearizePos(cfg, pos_trial.(sides{side}));
    
    cfg = [];
    cfg.binsize = 1;
    tc.(sides{side}) = MakeTC(cfg, spikes_trial.(sides{side}), linpos.(sides{side}));
end

%% Categorizing place cells
choice_tsd = tsd(0, metadata.coord.chp_cm, {'x', 'y'});

for side = 1:length(sides)
    cfg = [];
    cfg.Coord = coord.(sides{side});
    choice.(sides{side}) = LinearizePos(cfg, choice_tsd);
    fields.(sides{side}) = tc.(sides{side}).field_template_idx(tc.(sides{side}).field_loc > choice.(sides{side}).data(1));
end

[~, remove.left, remove.right] = intersect(fields.left, fields.right);
for side = 1:length(sides)
    fields.(sides{side})(remove.left) = [];
    fields.(sides{side}) = unique(fields.(sides{side}));
    spikes_side.(sides{side}) = SelectTS([], spikes, fields.(sides{side}));
end

%% Make a Q matrix
for side = 1:length(sides)
    cfg = [];
    cfg.win = 0.1;
    q.(sides{side}) = MakeQfromS2(cfg, spikes_side.(sides{side}), swr_evt);
end

%% Co-activation probabilities
for side = 1:length(sides)
    cfg = [];
    cfg.nShuffle = 1000;
    cfg.useMask = 1;
    cfg.outputFormat = 'vectorU';
    cooccur.(sides{side}) = CoOccurQ2(cfg, q.(sides{side}));
end

%% Plotting
p0_data = [nanmean(cooccur.left.p0), nanmean(cooccur.right.p0)];
colors = flipud(linspecer(2));
location = [1, 3];

figure(1);
for side = 1:length(p0_data)
    bar(location(side), p0_data(side), 'FaceColor', colors(side, :), 'EdgeColor', 'none');
    hold on;
end
set(gca, 'XLim', [location(1)-1, location(2)+1], 'XTick', location, 'XTickLabel', {'Left', 'Right'});
xlabel('Field location');
ylabel({'Proportion of'; 'SWRs active'});
title('Activation probability (p0)');


p3_data = [nanmean(cooccur.left.p3), nanmean(cooccur.right.p3)];
colors = flipud(linspecer(2));
location = [1, 3];

figure(2);
for side = 1:length(p3_data)
    bar(location(side), p3_data(side), 'FaceColor', colors(side, :), 'EdgeColor', 'none');
    hold on;
end
set(gca, 'XLim', [location(1)-1, location(2)+1], 'XTick', location, 'XTickLabel', {'Left', 'Right'});
xlabel('Field location');
ylabel({'Cell pair'; 'joint probability'});
title('Observed coactivity (p3)');


p4_data = [nanmean(cooccur.left.p4), nanmean(cooccur.right.p4)];
colors = flipud(linspecer(2));
location = [1, 3];

figure(2);
for side = 1:length(p4_data)
    bar(location(side), p4_data(side), 'FaceColor', colors(side, :), 'EdgeColor', 'none');
    hold on;
end
set(gca, 'XLim', [location(1)-1, location(2)+1], 'XTick', location, 'XTickLabel', {'Left', 'Right'});
xlabel('Field location');
ylabel({'SWR coactivation'; 'z-scored'});
title('Coactivation above chance levels (p4)');

%% All plots together
p_list = {'p0', 'p3', 'p4'};
titles = {'Activation probability (p0)', 'Observed coactivity (p3)', 'Coactivation above chance levels (p4)'};
ylabels = {{'Proportion of'; 'SWRs active'}, {'Cell pair'; 'joint probability'}, {'SWR coactivation'; 'z-scored'}};
sides = {'left', 'right'};
colors = flipud(linspecer(2));
location = [1, 2.5];
 
figure;
for prob = length(p_list):-1:1
    p_data(prob, 1:2) = [nanmean(cooccur.left.(p_list{prob})), nanmean(cooccur.right.(p_list{prob}))];
    fig(prob) = subplot(1, 3, prob);
    for side = 1:length(sides)
        bar(location(side), p_data(prob, side), 'FaceColor', colors(side,:), 'Edgecolor', 'none');
        hold on;
        
        title(titles{prob});
        ylabel(ylabels{prob});
        xlabel('Field location');
    end
end

set(fig, 'XLim', [location(1)-1, location(2)+1], 'XTick', location, 'XTickLabel', {'Left', 'Right'});
set(fig, 'PlotBoxAspectRatio', [1, 1, 1]);
maximize


%% Different trial numbers
fprintf('Number of left trials: %d\n', length(metadata.taskvars.trial_iv_L.tstart))
fprintf('Number of right trials: %d\n', length(metadata.taskvars.trial_iv_R.tstart))




