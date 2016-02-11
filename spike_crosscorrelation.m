function [crosscorrelation, bins] = spike_crosscorrelation(spike_time1, spike_time2, bin_size, max_time)

bin_centers = -max_time-bin_size:bin_size:max_time+bin_size;
crosscorrelation = zeros(size(bin_centers));

for spike = 1:length(spike_time1)
    relative_spike_time = spike_time2 - spike_time1(spike);
    crosscorrelation = crosscorrelation + hist(relative_spike_time, bin_centers);
end

bins = bin_centers(2:end-1);
crosscorrelation = crosscorrelation(2:end-1);

crosscorrelation = crosscorrelation./(length(spike_time1)); % Normalized
end