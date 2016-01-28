fs = 100;
time0 = 0; time1 = 1;
time = time0:1/fs:time1;

freq = 2;
%signal = sin(2*pi*freq*time);

%stem(time, signal);

phi = pi/2;

subplot(221);
signal = sin(2*pi*freq*time + phi);
stem(time, signal); hold on;
plot(time, cos(2*pi*freq*time), 'r--', 'LineWidth', 2);
legend('sin (phase-shifted)', 'cos');

amp = 2;
aubplot(222);
signal = a.*sin(2*pi*freq*time + phi);
stem(time, signal);

%% Stupid cell break
freq2 = 10;
change = 2;
 
subplot(311)
signal1 = sin(2*pi*freq*time);
plot(time, signal1); title('message');
 
subplot(312);
signal2 = sin(2*pi*freq2*time);
plot(time, signal2); title('carrier');
 
subplot(313);
signal3 = sin(2*pi*freq2*time + change.*sin(2*pi*freq*time - pi/2));
plot(time, signal3); title('FM signal');

%% Annoying cells
magnitudes = [0.1, 0, 1.3, 0.5];
phases = [-pi/6, 0, pi, 2*pi/3];
freq = 2;

signal_out = zeros(size(time));
for signal = 1:numel(magnitudes)
    current_signal = magnitudes(signal)*cos(2*pi*freq*signal*time + phases(signal));
    plot(time, current_signal, 'r:'); hold on;
    signal_out = signal_out + current_signal;
end
plot(time, signal_out, 'LineWidth', 2);

%% Annoying cell!
rng('default');

x = round(rand(1, 8)*10);
xlen = length(x);

x1 = fft(x);
x_mag = abs(x1);
x_phase = angle(x1);

n = 0:xlen-1;
time = 0:0.05:xlen-1;

for sig = xlen-1:-1:0
    val = sig+1;
    coarse(val,:) = x_mag(val)*cos(2*pi*n*sig/xlen + x_phase(val))/xlen;
    fine(val,:) = x_mag(val)*cos(2*pi*time*sig/xlen + x_phase(val))/xlen;
end

coarse_sum = sum(coarse);
fine_sum = sum(fine);

figure;
plot(n, x, 'mo', 'LineWidth', 2); hold on;
plot(time, fine_sum, 'b', 'LineWidth', 2);
plot(n, coarse_sum, 'r*', 'LineWidth', 2);
legend({'Original', 'Sum - All', 'Sum - Points only'});

%% Annoying cell!
fs = 20;
time0 = 0;
time1 = 1;
time = time0:1/fs:time1-(1/fs);

freq = 2;
signal = sin(2*pi*freq*time);

signal_fft = fft(signal, length(signal));
signal_fft_mag = abs(signal_fft);
signal_fft_phase = angle(signal_fft);
%stem(signal_fft_mag, 'LineWidth', 2);

npoints = length(signal);
fourier = [-npoints/2:npoints/2-1]./npoints;

signal_fft_mag = fftshift(signal_fft_mag);
stem(fourier, signal_fft_mag, 'LineWidth', 2);
xlabel('Frequency (Fs^{-1})');

%% Cell modeee
time = time0:1/fs:time1;
npoints = [length(time) 64 256 1024];

for point = 1:length(npoints)
    numpoint = npoints(point);
    subplot(2, 2, point);
    signal = sin(2*pi*freq*time);
    signal_fft = fft(signal, numpoint);
    signal_fft_mag = abs(signal_fft);
    signal_fft_phase = angle(signal_fft);
    
    fourier = [-numpoint/2:numpoint/2-1]./numpoint;
    signal_fft_mag = fftshift(signal_fft_mag);
    plot(fourier, signal_fft_mag, 'kx', fourier, signal_fft_mag, 'k');
    
    title(sprintf('%d point FFT', numpoint));
    xlabel('Frequency (Fs^{-1})');
end

%% Spectral leakage celll

time = time0:1/fs:time1-(1/fs);
num_repeat = [1, 2, 4, 8];

numpoints = 1024;

for point = 1:length(num_repeat)
    subplot(2, 2, point);
    
    signal = sin(2*pi*freq*time);
    signal = repmat(signal, [1, num_repeat(point)]);
 
    signal_fft = fft(signal, numpoint);
    signal_fft_mag = abs(signal_fft); 
    signal_fft_phase = angle(signal_fft);
 
    fourier = [-numpoint/2:numpoint/2-1]./numpoint;
    signal_fft_mag = fftshift(signal_fft_mag);
    plot(fourier, signal_fft_mag, 'kx', fourier, signal_fft_mag, 'k');
 
    title(sprintf('%d repeats',num_repeat(point)));
    xlabel('Frequency (Fs^{-1})');
end

%% Windowing
np = 25;
pfft = 1024;

windows = {'rectwin', 'triang', 'hamming', 'hanning', 'blackman'};
color = 'rgbcmyk';

for win = 1:length(windows)
    eval(cat(2, 'wn = ', windows{win}, '(np);'));
    wn = wn./sum(wn);
    
    subplot(211);
    plot(wn, color(win), 'LineWidth', 2); hold on;
    
    subplot(212);
    signal_fft = fft(win, pfft);
    signal_fft_mag = abs(signal_fft);
    signal_fft_phase = angle(signal_fft);
    
    fourier = [-pfft/2:pfft/2-1]./pfft;
    signal_fft_mag = fftshift(signal_fft_mag);
    
    meaningless(win) = plot(fourier, signal_fft_mag, color(win), 'LineWidth', 2);
    hold on;
end

xlabel('Frequency (Fs^{-1})');
legend(meaningless, windows);

%% Spectral estimation cell
[periodo, fourier] = periodogram(signal, [], np, fs);
plot(fourier, periodo);
xlabel('Frequency (Hz)'); hold on;
[periodo, fourier] = periodogram(signal, hanning(length(signal)), np, fs);
plot(fourier, periodo, 'r');

%% pwelch()
fs = 20;
time0 = 0; 
time1 = 1;
freq = 2;
num_repeat = 4;

time = time0:1/fs:time1-(1/fs);
np = 1024;
signal = sin(2*pi*freq*time);
signal = repmat(signal, [1, num_repeat]);

[periodo, fourier] = periodogram(signal, rectwin(length(signal)), np, fs);
plot(fourier, periodo); hold on;

window_size = 40;
[periodo, fourier] = pwelch(signal, rectwin(window_size), window_size/2, np, fs);
plot(fourier, periodo, 'r'); 
xlabel('Frequency (Hz)');

%% pitfalls to real world signals
fs = 20;
time0 = 0;
time1 = 1;
freq = 2;
np = 1024;
gaps = [5, 10, 15];

time = time0:1/fs:time1;
signal = sin(2*pi*freq*time);

subplot(211);
plot(time, signal, 'k*'); hold on;

signal_fft = fft(signal, np);
signal_fft_mag = abs(signal_fft);
signal_fft_phase = angle(signal_fft);

fourier = [-np/2:np/2-1]./np;
signal_fft_mag = fftshift(signal_fft_mag);

subplot(212);
plot(fourier, signal_fft_mag, 'm-x'); hold on;
xlabel('Frequency (Fs^{-1})');

signal = sin(2*pi*freq*time);
signal2 = signal;
signal2(gaps) = []; 
time(gaps) = [];

subplot(211);
plot(time, signal2, 'bo'); hold on;

signal_fft = fft(signal2, np);
signal_fft_mag = abs(signal_fft);
signal_fft_phase = angle(signal_fft);

fourier = [-np/2:np/2-1]./np;
signal_fft_mag = fftshift(signal_fft_mag);

subplot(212);
plot(fourier, signal_fft_mag, 'b-x');
legend('Signal', 'Signal with gaps');

%% Application to real data
cd('C:\Users\Emily\Desktop\R042-2013-08-18');
cfg = [];
cfg.fc = {'R042-2013-08-18-CSC11a.ncs'};
csc = LoadCSC(cfg);

csc_pre = restrict(csc,0,csc.cfg.ExpKeys.task(1)+10);

%plot(diff(csc_pre.tvec));
fs = 1./mean(diff(csc_pre.tvec));

factor = 4;
csc_pre.data = decimate(csc_pre.data, factor);
csc_pre.tvec = downsample(csc_pre.tvec, factor);
csc_pre.cfg.hdr{1}.SamplingFrequency = csc_pre.cfg.hdr{1}.SamplingFrequency./factor;

window = 1024;
[periodo, fourier] = periodogram(csc_pre.data, hamming(length(csc_pre.data)), length(csc_pre.data), fs);
plot(fourier, 10*log10(periodo), 'k');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
xlim([0,150]);


% Question: line 258/259. Why do we decimate the data but downsample the
% time?


    
    
    
    