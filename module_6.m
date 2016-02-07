%% Filtering
load count.dat;
x = count(:, 1);
a = 1;
b = [1/4, 1/4, 1/4, 1/4];
y = filter(b, a, x);

time = 1:length(x);
plot(time, x, 'r-.', time, y, 'b-'); grid on
legend('Original', 'Filtered', 2)

%% Butterworth
fs = 400;
wp = [50, 100] * 2 / fs; % passband - between 50 and 100 Hz
ws = [45, 105] * 2 / fs; % stopband
[n, wn] = buttord( wp, ws, 3, 20); % determine filter parameters
[b2, a2] = butter(n,wn); % builds filter
fvtool(b, a, b2, a2)

%% Chebyshev type I
fs = 400;
wp = [50, 100] * 2 / fs;
ws = [48, 102] * 2 / fs;
[n, wn] = cheb1ord(wp, ws, 3, 20);
[b_cheb1, a_cheb1] = cheby1(n, 0.5, wn);
fvtool(b2, a2, b_cheb1, a_cheb1)

%% Phase response
fs = 500; 
dt = 1./fs;
time_lim = [0, 10];
time = time_lim(1):dt:time_lim(2)-dt;
 
signal_1 = sin(2*pi*80*time+pi/6);
signal_2 = sin(2*pi*40*time);
signal = signal_1 + signal_2;
 
% filtered_signal = filter(b_cheb1, a_cheb1, signal);
%  
% plot(time,signal,'k',time,filtered_signal,'r--'); hold on;
% legend({'original','filtered'});
% xlim([0 0.2]);

filtered_signal = filtfilt(b_cheb1, a_cheb1, signal);
 
plot(time, signal, 'k', time, filtered_signal,'r--'); 
hold on;
legend({'original','filtered'});
xlim([0 0.2]);

%% compare freq responses
fs = 500; 
dt = 1./fs;
time_lim = [0 10];
time = time_lim(1):dt:time_lim(2)-dt;
 
x = rand(size(time)); % white noise input
[p, f] = pwelch(x, hanning(512), 256, 2^14, fs);
 
y1 = filter(b_cheb1, a_cheb1, x);
[p1, f1] = pwelch(y1, hanning(512), 256, 2^14, fs);
 
y2 = filtfilt(b_cheb1, a_cheb1, x);
[p1, f2] = pwelch(y2, hanning(512), 256, 2^14, fs);
 
plot(f, 10*log10(p), f, 10*log10(p1), f, 10*log10(p1));
legend({'original','filter','filtfilt'});

%% Remove 60Hz
fs = 500;
[b, a] = butter(10, [59, 61] * 2 / fs, 'stop');
fvtool(b, a);

[z, p, k] = butter(10, [59, 61] * 2 / fs, 'stop'); % note, we ask for 3 outputs instead of 2
[sos, g] = zp2sos(z, p, k); % convert to SOS format
h = dfilt.df2sos(sos, g); % create filter object
fvtool(h);


%% data application
cd('C:\Users\Emily\Desktop\R042-2013-08-18');
cfg = [];
cfg.fc = {'R042-2013-08-18-CSC08a.ncs'};
csc = LoadCSC(cfg);
 
cscR = restrict(csc, 3270, 3272);
plot(cscR.tvec, cscR.data)

fs = cscR.cfg.hdr{1}.SamplingFrequency;
wp = [180, 220] * 2 / fs;
ws = [178, 222] * 2 / fs;
[n, wn] = cheb1ord(wp, ws, 3, 20); % determine filter parameters
[b_cheb1, a_cheb1] = cheby1(n, 0.5, wn); % builds filter
 
% fvtool(b_cheb1, a_cheb1); % Check the filter
 
y = filtfilt(b_cheb1, a_cheb1, cscR.data);
plot(cscR.tvec, cscR.data, 'b', cscR.tvec, y, 'r');

chew_power = y.^2;
chew_power_filtered = medfilt1(chew_power, 101); % filter window is specified in samples, so this is ~50ms
[h1, h2] = plotyy(cscR.tvec, cscR.data, cscR.tvec, chew_power_filtered);

%% Filtering on real data

% cd('C:\Users\Emily\Desktop\R042-2013-08-18');
cd('C:\Users\Emily\Desktop\R063-2015-03-20_recording');
cfg_csc = [];
% cfg_csc.fc = {'R042-2013-08-18-CSC03a.ncs'};
cfg_csc.fc = {'R063-2015-03-20-CSC03a.ncs'};
csc = LoadCSC(cfg_csc);

cfg_swr = [];
cfg_swr.f = [140, 220];
cfg_swr.display_filter = 0;

filter_swr = FilterLFP(cfg_swr, csc);

% Obtain power and z-score
swr_power = LFPpower([], filter_swr);
swr_zpower = zscore_tsd(swr_power);

% Detect events
cfg_evt = [];
cfg_evt.method = 'raw';
cfg_evt.threshold = 4;
cft_evt.operation = '>';
cfg_evt.merge_thr = 0.05;
cfg_evt. minlen = 0.05;

swr_evt = TSDtoIV(cfg_evt, swr_zpower);

% Find max z-scored power for each event
cfg_maxp = [];
cfg_maxp.method = 'max';
cfg_maxp.label = 'maxSWRpower';

swr_evt = AddTSDtoIV(cfg_maxp, swr_evt, swr_zpower);

% Select events with >5 z-scored power
cfg_thres = [];
cfg_thres.operation = '>=';
cfg_thres.threshold = 5;

[swr_evt, ~] = SelectIV(cfg_thres, swr_evt, 'maxSWRpower');

% Plot swr events over all lfp
% PlotTSDfromIV([], swr_evt, csc);
% close all;

% Plot swr events alone
% cfg_center = [];
% cfg_center.display = 'iv';
% cfg_center.mode = 'center';
% cfg_center.fgcol = 'k';
% 
% PlotTSDfromIV(cfg_center, swr_evt, csc);

% Plot with highlighted swrs?
cfg_edge = [];
cfg_edge.display = 'iv';
cfg_edge.fgcol = 'r';
PlotTSDfromIV(cfg_edge, swr_evt, csc);





