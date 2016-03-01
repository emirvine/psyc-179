function [ u_time, shortcut_time, novel_time, feeder1_time, feeder2_time, ...
    phase_start_time, pedestal_time, exploratory_time, trial_label] ...
    = shortcut_trial_idx(pos_tsd, expkeys)
% * Returns u_idx, shortcut_idx, novel_idx for trials extracted based on 
%   event and position data from shortcut experiment
[sorted_times, sorted_labels] = shortcut_trial_seg(pos_tsd, expkeys);

feeder1_idx = find(strcmp(sorted_labels,'feeder1'));
feeder2_idx = find(strcmp(sorted_labels,'feeder2'));
phase_start_idx = find(strcmp(sorted_labels,'phase_start'));

shortcut_trials = 0;
novel_trials = 0;
u_trials = 0;
pedestal_trials = 0;
exploratory_trials = 0;
non_trial = 0;

shortcut_idx = [];
novel_idx = [];
u_idx = [];
pedestal_idx = [];
exploratory_idx = [];
non_idx = [];

f1 = 1;
f2 = 1;
ped = 1;
trial_label = {};

while f1 <= length(feeder1_idx) && f2 <= length(feeder2_idx)
    start_idx = pick_min(feeder1_idx,feeder2_idx,phase_start_idx,f1,f2,ped);
    if feeder1_idx(f1) == start_idx
        f1 = f1 + 1;
    elseif feeder2_idx(f2) == start_idx
        f2 = f2 + 1;
    elseif phase_start_idx(ped) == start_idx
        ped = ped + 1;
    end
    end_idx = pick_min(feeder1_idx,feeder2_idx,phase_start_idx,f1,f2,ped);
    trial = sorted_labels(start_idx:end_idx);
    if strcmp(trial(end), 'phase_start')
        pedestal_trials = pedestal_trials + 1;
        trial_label = [trial_label, 'phase_start'];
        pedestal_idx(1,pedestal_trials) = start_idx;
        pedestal_idx(2,pedestal_trials) = end_idx;
%     elseif (any(strcmp(trial,'novel1')) || any(strcmp(trial,'novel2'))) && ...
%            (any(strcmp(trial,'shortcut1')) || any(strcmp(trial,'shortcut2')))
%         exploratory_trials = exploratory_trials + 1;
%         trial_label = [trial_label, 'exploratory'];
%         exploratory_idx(1,exploratory_trials) = start_idx;
%         exploratory_idx(2,exploratory_trials) = end_idx;
    elseif any(strcmp(trial,'novel1')) || any(strcmp(trial,'novel2'));
        novel_trials = novel_trials + 1;
        trial_label = [trial_label, 'novel'];
        novel_idx(1,novel_trials) = start_idx;
        novel_idx(2,novel_trials) = end_idx;
    elseif any(strcmp(trial,'shortcut1')) || any(strcmp(trial,'shortcut2'));
        shortcut_trials = shortcut_trials + 1;
        trial_label = [trial_label, 'shortcut'];
        shortcut_idx(1,shortcut_trials) = start_idx;
        shortcut_idx(2,shortcut_trials) = end_idx;
    elseif any(strcmp(trial,'u1')) || any(strcmp(trial,'u2'));
        u_trials = u_trials + 1;
        trial_label = [trial_label, 'u'];
        u_idx(1,u_trials) = start_idx;
        u_idx(2,u_trials) = end_idx; 
    else
        non_trial = non_trial + 1;
        non_idx(1,non_trial) = start_idx;
        non_idx(2,non_trial) = end_idx;
    end     
end

feeder1_time = find_near(pos_tsd.tvec, sorted_times(feeder1_idx), 1);
feeder2_time = find_near(pos_tsd.tvec, sorted_times(feeder2_idx), 1);
phase_start_time = find_near(pos_tsd.tvec, sorted_times(phase_start_idx), 1);
pedestal_time = find_near(pos_tsd.tvec, sorted_times(pedestal_idx), 2);
exploratory_time = find_near(pos_tsd.tvec, sorted_times(exploratory_idx), 2);
u_time = find_near(pos_tsd.tvec, sorted_times(u_idx), 2);
shortcut_time = find_near(pos_tsd.tvec, sorted_times(shortcut_idx), 2);
novel_time = find_near(pos_tsd.tvec, sorted_times(novel_idx), 2);
end

function [min_val] = pick_min(list1, list2, list3, idx1, idx2, idx3)
if idx1 > length(list1)
    val1 = Inf;
else
    val1 = list1(idx1);
end
if idx2 > length(list2)
    val2 = Inf;
else
    val2 = list2(idx2);
end
if idx3 > length(list3)
    val3 = Inf;
else
    val3 = list3(idx3);
end
min_val = min([val1,val2,val3]);
end

function [time] = find_near(data, idx, ax)
if size(idx) == [1,2]
    idx = reshape(idx,[2,1]);
end
time = zeros(size(idx));
if not(isempty(idx))
    for i = 1:ax
        time(i,:) = data(nearest_idx3(idx(i,:), data));
    end
end
end