%cd('C:\Users\Emily\Desktop\R042-2013-08-18');
cd('C:\Users\Emily\Desktop\R066-2014-12-01_recording');

cfg = [];
S = LoadSpikes(cfg);

cfg = [];
cfg.fc = {'R066-2014-12-01-CSC03a.ncs'};
csc = LoadCSC(cfg);

%tstart = 5900;
%tend = 6000;

tstart = 10900;
tend = 11000;

s_r = restrict(S, tstart, tend);
csc_r = restrict(csc, tstart, tend);

neuron_num = 17;
this_spk = s_r.t{neuron_num};
%fprintf('%4.3f\n',this_spk)

% h = plot([this_spk, this_spk], [0.5, 1.5], 'k');
% axis([xlim, 0, 5])
set_spky = 0.4;

figure; hold on;
for n = 1:length(s_r.t)
    if ~isempty(s_r.t{n})
        plot([s_r.t{n}, s_r.t{n}],[n-set_spky, n+set_spky],'k');
    end
end
ylabel('Neuron number');

set_cscy = [-5, 0];
plot(csc_r.tvec, rescale(csc_r.data, set_cscy(1), set_cscy(2)), 'b');
set(gca, 'YTickLabel', {});

ax1 = gca; ax1_pos = get(ax1, 'Position');
ax2 = axes('Position', ax1_pos);

cfg = []; cfg.tvec = csc.tvec; cfg.sigma = 0.1;
mua = getMUA(cfg,S);

xr = get(ax1, 'XLim');
mua_r = restrict(mua, xr(1), xr(2));

axes(ax2);
mua_hd1 = plot(mua_r.tvec, mua_r.data, 'Color', [0.7,0.7,0.7]);

set(gca, 'YAxisLocation','right', 'Box','off', 'XTick',[], 'Color','none', 'YColor',[0.5,0.5,0.5])
ylabel('multi-unit activity (spk/s)');
linkaxes([ax1,ax2], 'x');

zoom_start = 8965;
zoom_end = 8969;

%set(gca, 'XLim', [zoom_start, zoom_end]);

%%
%cd('C:\Users\Emily\Desktop\R042-2013-08-18');
cd('C:\Users\Emily\Desktop\R066-2014-12-01_recording');
S = LoadSpikes([]);
csc = []; csc.fc = {'R066-2014-12-01-CSC03a.ncs'};
csc = LoadCSC(csc);
 
%LoadMetadata;
unique_folder = 'R066-2014-12-01_recording';
key_id = strcat(unique_folder(1:4),'_',unique_folder(6:9),'_',...
    unique_folder(11:12),'_',unique_folder(14:15));

unique_key = strcat(key_id,'_ExpKeys.m');
run(unique_key);

cfg_plot = [];
cfg_plot.lfp = csc;
cfg_plot.evt = ExpKeys.pauseA;
%cfg_plot.evt = metadata.taskvars.trial_iv_L;

h = MultiRaster(cfg_plot,S);

%% movie version
cd('C:\Users\Emily\Desktop\R042-2013-08-18');
t = [5900 6000]; % start and end times (experiment time)
FPS = 30; % frame rate (per s)
twin = [-1 1]; % width of time window (in s)
 
tvec = t(1):1/FPS:t(2);
for iT = 1:length(tvec)
 
   set(gca,'XLim',twin+tvec(iT));
   drawnow; pause(1/FPS);
 
end

