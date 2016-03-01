function [u_times, shortcut_times, novel_times] = get_trial_times(pos_tsd, expkeys)

event_list = {'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0001).', 'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0002).', 'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0000).'};
events_tsd = shortcut_loadEvents(event_list);

[u_times, shortcut_times, novel_times, feeder1_times, feeder2_times, ...
    phase_start_times, pedestal_times, exploratory_times, trial_label] = ...
    shortcut_trial_idx(pos_tsd, expkeys);

figure(2); hold on;
if not(isempty(u_times))
    plot(u_times(1,:),700,'r.','MarkerSize',15);
    u_length = u_times(2,:) - u_times(1,:);
else
    u_length = 0;
end
if not(isempty(shortcut_times))
    plot(shortcut_times(1,:),700,'k.','MarkerSize',15);
    shortcut_length = shortcut_times(2,:) - shortcut_times(1,:);
else
    shortcut_length = 0;
end
if not(isempty(novel_times))
    plot(novel_times(1,:),700,'m.','MarkerSize',15);
    novel_length = novel_times(2,:) - novel_times(1,:);
else
    novel_length = 0;
end
if not(isempty(pedestal_times))
    plot(pedestal_times(1,:),700,'g.','MarkerSize',15);
    pedestal_length = pedestal_times(2,:) - pedestal_times(1,:);
else
    pedestal_length = 0;
end
if not(isempty(phase_start_times))
    plot(phase_start_times(1,:),700,'c.','MarkerSize',14);
end
if not(isempty(exploratory_times))
    plot(exploratory_times(1,:),700,'b.','MarkerSize',14);
    exploratory_length = exploratory_times(2,:) - exploratory_times(1,:);
else
    exploratory_length = 0;
end

figure(3); clf; hold on;
title(['Box plot of position trial lengths (s)']);
lengths = [u_length, shortcut_length, novel_length];
group = [zeros(1,length(u_length)), ones(1,length(shortcut_length)), ...
    ones(1,length(novel_length))*2];
boxplot(lengths, group);
set(gca, 'xtick', [1, 2, 3]);
n_u = strcat('n=', num2str(length(u_length)));
n_shortcut = strcat('n=', num2str(length(shortcut_length)));
n_novel = length(novel_length) + length(exploratory_length);
n_novel = strcat('n=', num2str(length(novel_length)));
set(gca,'xtickLabel', {['Full U, ', n_u], ['Shortcut, ', n_shortcut], ...
    ['Novel, ', n_novel]});

if not(isempty(shortcut_times))
    shortcut_start_idx = nearest_idx3(shortcut_times(1,:), pos_tsd.tvec);
    shortcut_end_idx = nearest_idx3(shortcut_times(2,:), pos_tsd.tvec);
end
u_start_idx = nearest_idx3(u_times(1,:), pos_tsd.tvec);
u_end_idx = nearest_idx3(u_times(2,:), pos_tsd.tvec);
if not(isempty(novel_times))
    novel_start_idx = nearest_idx3(novel_times(1,:), pos_tsd.tvec);
    novel_end_idx = nearest_idx3(novel_times(2,:), pos_tsd.tvec);
end
pedestal_start_idx = nearest_idx3(pedestal_times(1,:), pos_tsd.tvec);
pedestal_end_idx = nearest_idx3(pedestal_times(2,:), pos_tsd.tvec);
if not(isempty(exploratory_times))
    exploratory_start_idx = nearest_idx3(exploratory_times(1,:), pos_tsd.tvec);
    exploratory_end_idx = nearest_idx3(exploratory_times(2,:), pos_tsd.tvec);
end

figure(4); clf; hold on; axis off;
%title(['Position of ',unique_folder(1:4),', whole experiment']);
for seg = 1:length(shortcut_start_idx)
    plot(pos_tsd.data(1,shortcut_start_idx(seg):shortcut_end_idx(seg)),...
        pos_tsd.data(2,shortcut_start_idx(seg):shortcut_end_idx(seg)), 'r.', 'MarkerSize', 1);
end
for seg = 1:length(u_start_idx)
    plot(pos_tsd.data(1,u_start_idx(seg):u_end_idx(seg)),...
        pos_tsd.data(2,u_start_idx(seg):u_end_idx(seg)), 'b.', 'MarkerSize', 1);
end
for seg = 1:length(novel_start_idx)
    plot(pos_tsd.data(1,novel_start_idx(seg):novel_end_idx(seg)),...
        pos_tsd.data(2,novel_start_idx(seg):novel_end_idx(seg)), 'k.', 'MarkerSize', 1);
end
for seg = 1:length(pedestal_start_idx)
    plot(pos_tsd.data(1,pedestal_start_idx(seg):pedestal_end_idx(seg)),...
        pos_tsd.data(2,pedestal_start_idx(seg):pedestal_end_idx(seg)),'g.','MarkerSize', 1);
end

u_times = iv(u_times(1,:),u_times(2,:));
shortcut_times = iv(shortcut_times(1,:),shortcut_times(2,:));
novel_times = iv(novel_times(1,:),novel_times(2,:));
if not(isempty(exploratory_times))
    exploratory_times = iv(exploratory_times(1,:),exploratory_times(2,:));
else
    exploratory_times = iv(0,0);
end
