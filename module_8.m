cd('C:\Users\Emily\Desktop\R042-2013-08-18');

spikes = LoadSpikes([]);

neuron = 47;
time = [5801, 5801.7];

restricted_spikes = restrict(spikes, time(1), time(2));

plot([restricted_spikes.t{neuron}, restricted_spikes.t{neuron}], [-1,-0.5], 'Color', [0,0,0]);

bin_size = 0.01;
time_edges = time(1):bin_size:time(2);
time_centers = time_edges(1:end-1) + bin_size/2;

spike_count = histc(restricted_spikes.t{neuron}, time_edges);
spike_count = spike_count(1:end-1);

hold on;

h = bar(time_centers, spike_count);
set(h, 'BarWidth', 1, 'EdgeColor', 'none', 'FaceColor', [1,0,0])

yticks = get(gca, 'YTick');
ylabels = get(gca, 'YTickLabel');
set(gca, 'YLim', [-1.5, 10], 'YTick', yticks(2:end), 'YTickLabel', ylabels(2:end,:));
xlabel('Time (s)');

% Convolving spike train with a gaussian
bin_size = 0.001;
time_edges = time(1):bin_size:time(2);
time_centers = time_edges(1:end-1) + bin_size/2;

spike_count = histc(restricted_spikes.t{neuron}, time_edges);
spike_count = spike_count(1:end-1);

rect_spikes = conv2(spike_count, rectwin(50), 'same');
plot(time_centers, rect_spikes, 'b'); hold on;

gauss_spikes = conv2(spike_count, gausswin(50), 'same');
plot(time_centers, gauss_spikes, 'm');
%%
plot([restricted_spikes.t{neuron}, restricted_spikes.t{neuron}], [-1,-0.5], 'Color', [0,0,0]); 
hold on;

bin_size = 0.001;
gauss_window = 1./bin_size;
gauss_std = 0.02./bin_size;
gauss_kernel = gausskernel(gauss_window, gauss_std);
gauss_kernel = gauss_kernel./bin_size;
gaussk_spikes = conv2(spike_count, gauss_kernel, 'same');
plot(time_centers, gaussk_spikes, 'g');

%% Peri-event of -stimulus time histograms
neuron = 47;
spike_time = spikes.t{neuron};
isi = diff(spike_time);

dt = 0.001;
isi_edges = 0:dt:0.25;
isi_centers = isi_edges(1:end-1) + dt./2;
isi_count = histc(isi, isi_edges);

bar(isi_centers, isi_count(1:end-1));
set(gca, 'FontSize', 16, 'XLim', [0,0.25]);
xlabel('ISI (s)');
ylabel('Count');
grid on;

%% Poisson point process
dt = 0.001;
interval = [0, 10];
time = interval(1):dt:interval(2);

prob_spike = 0.5;
rng default;
spike_poisson = rand(size(time));
spike_poisson_idx = find(spike_poisson < prob_spike);
spike_poisson_time = time(spike_poisson_idx)';

line([spike_poisson_time, spike_poisson_time],[-1,-0.5],'Color',[0,0,0]);
axis([0, 0.1, -1.5, 5]);
set(gca,'YTick',[]);
%% Not working to bin these poisson spikes
bin_size = 0.001;
time_edges = time(1):bin_size:time(2);
time_centers = time_edges(1:end-1) + bin_size/2;

spike_count = histc(spike_poisson_time, time_edges);
spike_count = spike_count(1:end-1);

hold on;

h = bar(time_centers, spike_count);
set(h, 'BarWidth', 1, 'EdgeColor', 'none', 'FaceColor', [1,0,0])

yticks = get(gca, 'YTick');
ylabels = get(gca, 'YTickLabel');
set(gca, 'YLim', [-1.5, 10], 'YTick', yticks(2:end), 'YTickLabel', ylabels(2:end,:));
xlabel('Time (s)');

%% Spike autocorrelation function
cd('C:\Users\Emily\Dropbox\Graduate courses\psyc-179');
bin_size = 0.01;
max_time = 1;
spike_time = spike_poisson_time;
[autocorrelation, auto_bins] = spike_autocorrelation(spike_time, bin_size, max_time);

plot(auto_bins, autocorrelation); % Doesn't give expected output...

%% Spike cross-correlation function
cd('C:\Users\Emily\Dropbox\Graduate courses\psyc-179');
bin_size = 0.01;
max_time = 1;

neuron1 = 5;
neuron2 = 42;

restricted_spikes = restrict(spikes, 3200, 5650);
spike_time1 = restricted_spikes.t{neuron1};
spike_time2 = restricted_spikes.t{neuron2};

[crosscorrelation, cross_bins] = spike_crosscorrelation(spike_time1, spike_time2, bin_size, max_time);

plot(cross_bins, crosscorrelation);
set(gca, 'FontSize', 20);
xlabel('Lag (s)');
ylabel('Crosscorrelation');
title(sprintf('%d-%d', neuron1, neuron2));

% generate tc from poisson spikes, then fake place field




