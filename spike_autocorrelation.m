function [autocorrelation, bins] = spike_autocorrelation(spike_time, bin_size, max_time)

bin_centers = -max_time-bin_size:bin_size:max_time+bin_size;
autocorrelation = zeros(size(bin_centers));

for spike = 1:length(spike_time)
    relative_spike_time = spike_time - spike_time(spike);
    autocorrelation = autocorrelation + hist(relative_spike_time, bin_centers);
end

bins = bin_centers(2:end-1);
autocorrelation = autocorrelation(2:end-1);

autocorrelation = autocorrelation./max(autocorrelation); % Normalized