cd('C:\Users\Emily\Desktop\R042-2013-08-18');
cfg_spikes = [];
cfg_spikes.load_questionable_cells = 1;
spikes = LoadSpikes(cfg_spikes);

position = LoadPos([]);

LoadExpKeys;
spikes = restrict(spikes, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);
position = restrict(position, ExpKeys.TimeOnTrack, ExpKeys.TimeOffTrack);

plot(getd(position,'x'), getd(position,'y'),'.', 'Color',[0.5,0.5,0.5],'MarkerSize',1);
axis off;
hold on;

neuron = 7;
spike_x = interp1(position.tvec, getd(position,'x'),spikes.t{neuron},'linear');
spike_y = interp1(position.tvec, getd(position,'y'),spikes.t{neuron},'linear');

h = plot(spike_x, spike_y, '.r');


% Estimating tuning curves
LoadMetadata;

encode_spikes = restrict(spikes, metadata.taskvars.trial_iv);
encode_position = restrict(position, metadata.taskvars.trial_iv);

% Remove empty spikes
spiking_idx = ~cellfun(@isempty, encode_spikes.t);
encode_spikes.t = encode_spikes.t(spiking_idx);
encode_spikes.label = encode_spikes.label(spiking_idx);

spikes.t = spikes.t(spiking_idx);
spikes.label = spikes.label(spiking_idx);

% tuning curves
clear pos;
pos(:,1) = getd(encode_position,'y'); % Why are the 1,2 - x,y switched here?
pos(:,2) = getd(encode_position,'x');

x_min = 80;
y_min = 0;
x_max = 660;
y_max = 520;
x_binsize = 10;
y_binsize = 10;

x_edges = x_min:x_binsize:x_max;
y_edges = y_min:y_binsize:y_max;

occupancy = histcn(pos, y_edges, x_edges);

no_occupancy_idx = find(occupancy == 0);
occupancy(no_occupancy_idx) = NaN;

occupancy = occupancy .* (1/30);

subplot(221);
pcolor(occupancy); 
shading flat; 
axis off; 
colorbar
title('Occupancy');


% Basic spike histogram
clear spk;
neuron = 7;
spike_x = interp1(encode_position.tvec, getd(encode_position,'x'), encode_spikes.t{neuron}, 'linear');
spike_y = interp1(encode_position.tvec, getd(encode_position,'y'), encode_spikes.t{neuron}, 'linear');

spk(:,2) = spike_x; 
spk(:,1) = spike_y;

spike_hist = histcn(spk, y_edges, x_edges);
spike_hist(no_occupancy_idx) = NaN;

subplot(222);
pcolor(spike_hist);
shading flat;
axis off;
colorbar
title('Spikes');

% Rate map
tc = spike_hist ./ occupancy;

subplot(223);
pcolor(tc);
shading flat;
axis off;
colorbar
title('Rate map');

%% Smoothing
kernel = gausskernel([4,4],2);

[occupancy, ~, ~, position_idx] = histcn(pos, y_edges, x_edges);
occupancy = conv2(occupancy, kernel, 'same');

occupancy(no_occupancy_idx) = NaN;
occupancy = occupancy .* (1/30);

subplot(221);
pcolor(occupancy); 
shading flat;
axis off;
colorbar
title('Occupancy')

spike_hist = histcn(spk, y_edges, x_edges);
spike_hist = conv2(spike_hist, kernel, 'same');
spike_hist(no_occupancy_idx) = NaN;

subplot(222);
pcolor(spike_hist);
shading flat;
axis off;
colorbar
title('Spikes');

tc = spike_hist ./ occupancy;

subplot(223);
pcolor(tc); 
shading flat;
axis off;
colorbar
title('Rate map');


%% Tuning curves for all neurons
clear tc all_tc
num_neurons = length(encode_spikes.t);
for neuron = 1:num_neurons
    spike_x = interp1(encode_position.tvec, getd(encode_position,'x'), encode_spikes.t{neuron}, 'linear');
    spike_y = interp1(encode_position.tvec, getd(encode_position,'y'), encode_spikes.t{neuron}, 'linear');
    
    clear spk;
    spk(:,2) = spike_x;
    spk(:,1) = spike_y;
    spike_hist = histcn(spk, y_edges, x_edges);
    spike_hist = conv2(spike_hist, kernel, 'same');
    
    spike_hist(no_occupancy_idx) = NaN;
    
    tc = spike_hist ./ occupancy;
    all_tc{neuron} = tc;
end

fig_plot = 25;
for neuron = 1:length(encode_spikes.t)
    num_fig = ceil(neuron/fig_plot);
    figure(num_fig);
    
    subtightplot(5, 5, neuron-(num_fig-1)*fig_plot);
    pcolor(all_tc{neuron});
    shading flat;
    axis off;
    caxis([0,10]);
end


%% Prep for decoding (making q matrix and stuff)
clear q_matrix;
binsize = 0.25;

time_edges = metadata.taskvars.trial_iv.tstart(1):binsize:metadata.taskvars.trial_iv.tend(end);
q_time_centers = time_edges(1:end-1) + binsize/2;

for neuron = length(encode_spikes.t):-1:1
    spike_t = encode_spikes.t{neuron};
    q_matrix(neuron,:) = histc(spike_t, time_edges);
    q_matrix(neuron, end-1) = q_matrix(neuron, end-1) + q_matrix(neuron,end);
end

q_matrix = q_matrix(:,1:end-1);

imagesc(q_time_centers,1:num_neurons,q_matrix)
set(gca,'FontSize',16);
xlabel('Time (s)');
ylabel('Neuron number');

q_tsd = tsd(q_time_centers, q_matrix);
q_tsd = restrict(q_tsd, metadata.taskvars.trial_iv);

clear tc;
num_bins = numel(occupancy);
num_neurons = length(spikes.t);
for neuron = num_neurons:-1:1
    tc(:,:,neuron) = all_tc{neuron};
end
tc = reshape(tc,[size(tc,1)*size(tc,2), size(tc,3)]);
occupancy_norm = repmat(1/num_bins,[num_bins,1]);


%% Actually decoding
q_time_centers = q_tsd.tvec;
q_matrix = q_tsd.data;
num_active_neurons = sum(q_matrix > 0);

len = length(q_time_centers);
product = nan(length(q_time_centers), num_bins);

for bin = 1:num_bins
    temp_product = nansum(log(repmat(tc(bin,:)',1,len) .^ q_matrix));
    temp_sum = exp(-binsize*nansum(tc(bin,:),2));
    product(:,bin) = exp(temp_product) * temp_sum * occupancy_norm(bin);
end

product = product ./ repmat(sum(product,2),1,num_bins);
product(num_active_neurons < 1,:) = 0;


%% Visual inspection
x_binned = interp1(encode_position.tvec, position_idx(:,1), q_time_centers);
y_binned = interp1(encode_position.tvec, position_idx(:,2), q_time_centers);

good_occupancy = find(occupancy > 0);
x_bins = length(x_edges)-1;
y_bins = length(y_edges)-1;

decode_error = nan(length(q_time_centers),1);

for bin = 1:length(q_time_centers)
    cla;
    bad_name = reshape(product(bin,:),[y_bins, x_bins]);
    to_plot = nan(y_bins, x_bins);
    to_plot(good_occupancy) = bad_name(good_occupancy);
    
    pcolor(to_plot); 
    axis xy;
    hold on;
    caxis([0,0.5]);
    shading flat;
    axis off;
    
    hold on;
    plot(y_binned(bin), x_binned(bin), 'ow','MarkerSize',15);
    
    % Get xy coordinates of "MAP"
    [~,idx] = max(to_plot(:));
    [x_map, y_map] = ind2sub(size(to_plot), idx);
    
    if num_active_neurons(bin) > 0
        decode_error(bin) = sqrt((y_binned(bin)-y_map) .^2 + (x_binned(bin)-x_map) .^2);
    end
    
    plot(y_map, x_map, 'y*', 'MarkerSize', 5);
    
    h = title(sprintf('t %.2f, num_neurons %d, dist %.2f', q_time_centers(bin), num_active_neurons(bin), decode_error(bin)));
    
    if num_active_neurons(bin) == 0
        set(h, 'Color', [1,0,0]);
    else
        set(h, 'Color', [0,0,0]);
    end
    drawnow; pause(0.1);
end
    

%% Quantifying decoder accuracy
decode_error = nan(length(q_time_centers),1);

x_bins = length(x_edges)-1;
y_bins = length(y_edges)-1;

x_binned = interp1(encode_position.tvec, position_idx(:,1), q_time_centers);
y_binned = interp1(encode_position.tvec, position_idx(:,2), q_time_centers);

for bin = 1:length(q_time_centers)
    placeholder = reshape(product(bin,:),[y_bins, x_bins]);
    to_plot = nan(y_bins, x_bins);
    to_plot(good_occupancy) = placeholder(good_occupancy);
    
    [~,idx] = max(to_plot(:));
    [x_map, y_map] = ind2sub(size(to_plot), idx);
    
    if num_active_neurons(bin) > 0
        decode_error(bin) = sqrt((y_binned(bin)-y_map) .^2 + (x_binned(bin)-x_map) .^2);
    end
end


%% Average by lap
trial_id = zeros(size(q_time_centers));
run_start = ExpKeys.TimeOnTrack;
trial_idx = nearest_idx3(run_start, q_time_centers);
trial_id(trial_idx) = 1;
trial_id = cumsum(trial_id);

figure;
set(gca,'FontSize',18);
boxplot(decode_error, trial_id);
xlabel('Trial');
ylabel('Decoding error (pixels)');

avg_error = nanmean(decode_error);
title(sprintf('Average error %.2f', avg_error));