laptop_data = 'C:\Users\Emily\Desktop';
work_data = 'E:\data-shortcut\data-working\Shortcut-20150727';
laptop_code = 'C:\Users\Emily\Code\Shortcut';
work_code = 'E:\code\shortcut';
path_code = laptop_code;
path_data = laptop_data;

% Loading the data. Here from R063-2015-03-20_recording.
rat_id = 'R063_EI';
cd(fullfile(path_code, 'expkeys'))
expday = emi_expday(rat_id);
unique_folder = expday.two;
unique_id = unique_folder(1:15);

expkeys = emi_loadExpKeys(unique_folder);

cd(fullfile(path_data, unique_folder));
start_phase = expkeys.phase2(1);
end_phase = expkeys.phase2(2);

% Loading spikes
cfg = [];
spikes = LoadSpikes(cfg);
spikes_phase = restrict(spikes, start_phase, end_phase);

% Get XY position
if exist([unique_id,'-emi-vt.mat'],'file');
    fprintf('*-emi-vt.mat file found, loading.\n');
    load([unique_id,'-emi-vt.mat']);
else
    pos_tsd = emi_position(unique_folder,expkeys);
end

% Linear position
boundary1 = 'feeder1';
boundary2 = 'feeder2';
linear_z = emi_linearize(pos_tsd, expkeys, start_phase, end_phase, boundary1, boundary2);

% Making the tuning curves without filters
z_min = min(linear_z.data(1,:));
z_max = max(linear_z.data(1,:));
binsize = 3;
cfg = [];
cfg.binEdges{1} = z_min:binsize:z_max;
% cfg.smoothingKernel = gausskernel(15,3); % (wide by std dev)
clear tc;
tc = TuningCurves(cfg, spikes_phase, linear_z);

figure(1); clf; hold on;
plot(tc.tc', 'LineWidth', 1.2);
ylabel('Firing rate (Hz)');
xlabel('Position (cm)');
title('Unfiltered tuning curves in Matlab');

% Making the tuning curves with a Gaussian filter
z_min = min(linear_z.data(1,:));
z_max = max(linear_z.data(1,:));
binsize = 3;
cfg = [];
cfg.binEdges{1} = z_min:binsize:z_max;
cfg.smoothingKernel = gausskernel(15, 3);
clear tc;
tc_gauss = TuningCurves(cfg, spikes_phase, linear_z);

figure(2); clf; hold on;
plot(tc_gauss.tc', 'LineWidth', 1.2);
ylabel('Firing rate (Hz)');
xlabel('Position (cm)');
title('Gaussian-filtered tuning curves in Matlab');
