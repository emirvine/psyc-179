freq_sample = 1000;
freq1 = 10;
window = [0, 1];

time1 = window(1):1/freq_sample:window(2);
signal1 = sin(2*pi*freq1*time1);

% plot(signal1);

freq2 = 12;
time2 = window(1):1/freq2:window(2);
signal2 = interp1(time1, signal1, time2, 'nearest');

plot(time1, signal1); hold on;
plot(time2, signal2, '.g', 'MarkerSize', 20);
plot(time1, -sin(2*pi*2*time1), 'r--', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('y');


%% Urg cell blocks are annoying!
freq_sample = 1200;
freq1 = 3;
freq2 = 10;
window = [0, 1];

time1 = window(1):1/freq_sample:window(2);
signal1a = sin(2*pi*freq1*time1);
signal1b = 0.5*sin(2*pi*freq2*time1);
signal1 = signal1a + signal1b;

dt = 100;
time2 = time1(1:dt:end);
signal2 = signal1(1:dt:end);

subplot(131)
plot(time1, signal1); hold on;
stem(time2, signal2, '.r', 'MarkerSize', 20);
title('W/out anti-aliasing filter');

time1d = decimate(time1, dt);
signal2d = decimate(signal1, dt);
subplot(132);
plot(time1, signal1a, 'b--'); hold on;
stem(time1d, signal2d, '.r', 'MarkerSize', 20);
xlabel('Time (s)'); ylabel('y');
title('With anti-aliasing filter');

subplot(133);
plot(time2, signal2-signal2d, 'g', 'LineWidth', 2);
title('Difference');


%% Stupid cells..
freq_sample = 2000;
time = 0:1/freq_sample:4;

freq1 = 100;
signal = sin(2*pi*freq1*time);

ax1 = subplot(211);
stem(time, signal);
title('Original'); hold on;

subsample_factor = 4;

time2 = time(1:subsample_factor:end);
signal2 = signal(1:subsample_factor:end);

ax2 = subplot(212);
stem(time2, signal2, 'r');
xlabel('Time (s)');
title('Subsampled');

xlim = [1, 1.04];
linkaxes([ax1, ax2], 'x');
set(ax1, 'XLim', xlim); 

hold on;

signal_interp = interp1(time2, signal2, time, 'linear');
plot1 = plot(time, signal_interp, 'b');

signal_interp2 = interp1(time2, signal2, time, 'spline');
plot2 = plot(time, signal_interp2, 'm');

legend([plot1, plot2], {'Linear', 'Spline'}, 'Location', 'Northeast');
legend boxoff


%% Annnoying cells!
cd('C:\Users\Emily\Desktop\R063-2015-03-20_recording');
filename = 'R063-2015-03-20-CSC01c.ncs';
[Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(filename, [1 1 1 1 1], 1, 1, []);

plot(diff(Timestamps))




