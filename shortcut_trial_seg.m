function [ sorted_times, sorted_labels ] = shortcut_trial_seg(pos_tsd, expkeys)

ru1 = expkeys.u1;
ru2 = expkeys.u2;
rshortcut1 = expkeys.shortcut1;
rshortcut2 = expkeys.shortcut2;
rnovel1 = expkeys.novel1;
rnovel2 = expkeys.novel2;

labels = {};
times = [];

count = 1;

rect_u1 = pos_tsd.data(1,:) >= ru1(1) & pos_tsd.data(1,:) <= ru1(1)+ru1(3) & ...
    pos_tsd.data(2,:) >= ru1(2) & pos_tsd.data(2,:) <= ru1(2)+ru1(4);
rect_u1 = diff(rect_u1);
event_u1 = pos_tsd.tvec(rect_u1 == 1);
times = [times, event_u1];
labels(count:count+length(event_u1)-1) = {'u1'};
count = count + length(event_u1);

rect_u2 = pos_tsd.data(1,:) >= ru2(1) & pos_tsd.data(1,:) <= ru2(1)+ru2(3) & ...
    pos_tsd.data(2,:) >= ru2(2) & pos_tsd.data(2,:) <= ru2(2)+ru2(4);
rect_u2 = diff(rect_u2);
event_u2 = pos_tsd.tvec(rect_u2 == 1);
times = [times, event_u2];
labels(count:count+length(event_u2)-1) = {'u2'};
count = count + length(event_u2);

rect_shortcut1 = pos_tsd.data(1,:) >= rshortcut1(1) & ...
    pos_tsd.data(1,:) <= rshortcut1(1)+rshortcut1(3) & ...
    pos_tsd.data(2,:) >= rshortcut1(2) & ...
    pos_tsd.data(2,:) <= rshortcut1(2)+rshortcut1(4);
event_shortcut1 = diff(rect_shortcut1);
event_shortcut1 = pos_tsd.tvec(event_shortcut1 == 1);
times = [times, event_shortcut1];
labels(count:count+length(event_shortcut1)-1) = {'shortcut1'};
count = count + length(event_shortcut1);

rect_shortcut2 = pos_tsd.data(1,:) >= rshortcut2(1) & ...
    pos_tsd.data(1,:) <= rshortcut2(1)+rshortcut2(3) & ...
    pos_tsd.data(2,:) >= rshortcut2(2) & ...
    pos_tsd.data(2,:) <= rshortcut2(2)+rshortcut2(4);
event_shortcut2 = diff(rect_shortcut2);
event_shortcut2 = pos_tsd.tvec(event_shortcut2 == -1);
times = [times, event_shortcut2];
labels(count:count+length(event_shortcut2)-1) = {'shortcut2'};
count = count + length(event_shortcut2);

rect_novel1 = pos_tsd.data(1,:) >= rnovel1(1) & ...
    pos_tsd.data(1,:) <= rnovel1(1)+rnovel1(3) & ...
    pos_tsd.data(2,:) >= rnovel1(2) & ...
    pos_tsd.data(2,:) <= rnovel1(2)+rnovel1(4);
event_novel1 = diff(rect_novel1);
event_novel1 = pos_tsd.tvec(event_novel1 == 1);
times = [times, event_novel1];
labels(count:count+length(event_novel1)-1) = {'novel1'};
count = count + length(event_novel1);

rect_novel2 = pos_tsd.data(1,:) >= rnovel2(1) & ...
    pos_tsd.data(1,:) <= rnovel2(1)+rnovel2(3) & ...
    pos_tsd.data(2,:) >= rnovel2(2) & ...
    pos_tsd.data(2,:) <= rnovel2(2)+rnovel2(4);
event_novel2 = diff(rect_novel2);
event_novel2 = pos_tsd.tvec(event_novel2 == 1);
times = [times, event_novel2];
labels(count:count+length(event_novel2)-1) = {'novel2'};
count = count + length(event_novel2);

fn = FindFile('*Events.nev');
[EVTimeStamps, ~, ~, ~, EventStrings, ~] = Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);

time_convFactor = 10^-6; % converts nlx units to seconds
EVTimeStamps = EVTimeStamps * time_convFactor;

feeder1_idx = find(strcmp(expkeys.feeder1id, EventStrings));
feeder1_times = EVTimeStamps(feeder1_idx);
times = [times, feeder1_times];
labels(count:count+length(feeder1_times)-1) = {'feeder1'};
count = count + length(feeder1_times);

feeder2_idx = find(strcmp(expkeys.feeder2id, EventStrings));
feeder2_times = EVTimeStamps(feeder2_idx);
times = [times, feeder2_times];
labels(count:count+length(feeder2_times)-1) = {'feeder2'};
count = count + length(feeder2_times);


phase_start_times = [expkeys.phase1(1),expkeys.phase2(1),expkeys.phase3(1)];
times = [times, phase_start_times];
labels(count:count+length(phase_start_times)-1) = {'phase_start'};
count = count + length(phase_start_times);

[sorted_times, sort_idx] = sort(times);
sorted_labels = labels(sort_idx);
end