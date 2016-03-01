function [ events_ts ] = shortcut_loadEvents( event_list )
% requires ts.m, Nlx2MatEV
% returns events_ts
% Usage:
% event_list = {'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0001).', 'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0002).', 'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0000).'};
% events_tsd = shortcut_loadEvents(event_list);
% * Returns events_ts object (with tvec from 0 and data and header) 

time_convFactor = 10^-6; % converts nlx units to seconds

events_ts = ts;

fn = FindFile('*Events.nev');

[EVTimeStamps, EventIDs, TTLs, EVExtras, EventStrings, EVHeader] = ...
    Nlx2MatEV(fn,[1 1 1 1 1],1,1,[]);

for evt = 1:length(event_list)
    ev_string = event_list{evt};
    ev_id = strncmp(ev_string,EventStrings,length(ev_string));
    ev_t = EVTimeStamps(ev_id) * time_convFactor;
    events_ts.t{evt} = ev_t;
    events_ts.label{evt} = event_list;
end
% save([unique_folder(1:15),'-event-ts.mat'],'EVTimeStamps','EventStrings');
end