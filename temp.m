laptop_data = 'C:\Users\Emily\Desktop';
work_data = 'E:\data-shortcut\data-working\Shortcut-20150727';
laptop_code = 'C:\Users\Emily\Code\Shortcut';
work_code = 'E:\code\shortcut';
path_code = laptop_code;
path_data = laptop_data;

rat_id = 'R063_EI';
cd(fullfile(path_code, 'expkeys'))
expday = emi_expday(rat_id);
unique_folder = expday.three;
unique_id = unique_folder(1:15);

expkeys = emi_loadExpKeys(unique_folder);

% cd(fullfile(path_data, rat_id, unique_folder));
cd(fullfile(path_data, unique_folder));
start_phase = expkeys.phase2(1);
end_phase = expkeys.phase2(2);
boundary1 = 'feeder1';
boundary2 = 'feeder2';

% Get position
if exist([unique_id,'-emi-vt.mat'],'file');
    fprintf('*-emi-vt.mat file found, loading.\n');
    load([unique_id,'-emi-vt.mat']);
else
    pos_tsd = emi_position(unique_folder,expkeys);
end

% Linear position
[z, z_iv] = emi_linearize(pos_tsd, expkeys, start_phase, end_phase, boundary1, boundary2);

% Tuning doesn't use z_dist
z.u.data = z.u.data(1,:);

zlin = z.u;
zlin_iv = z_iv.u;

% Filtering neurons using emi_spikes_filtered.m
cfg = [];
cfg.load_questionable_cells = 1;
cfg.useClustersFile = 0;
spikes = LoadSpikes(cfg);
% spikes_filtered = spikes_filtered_shortcut(spikes, pos_tsd, zlin, zlin_iv, expkeys);


tstart = expkeys.phase2(1);
tend = expkeys.phase2(1)+30;

spikes_restricted = restrict(spikes, tstart, tend);
set_spky = 0.4;
figure; hold on;
for n = 1:length(spikes_restricted.t)
    if ~isempty(spikes_restricted.t{n})
        plot([spikes_restricted.t{n}, spikes_restricted.t{n}], [n-set_spky, n+set_spky],'k');
    end
end
ylabel('Neuron number');

% Making the tuning curves
z_min = min(zlin.data(1,:));
z_max = max(zlin.data(1,:));
binsize = 0.1;
cfg = [];
cfg.binEdges{1} = z_min:binsize:z_max;
cfg.smoothingKernel = gausskernel(15,3); % (wide by std dev)
clear tc;
tc = TuningCurves(cfg, spikes, zlin);

% check = 12;
figure(1); clf; hold on;
plot(tc.tc);
% imagescnan(squeeze(tc.tc(check,:)));
