cd('C:\Users\Emily\Desktop\R042-2013-08-18');

%% input_csc
cfg_csc = [];
cfg_csc.fc = {'R042-2013-08-18-CSC11a.ncs'};
csc = LoadCSC(cfg_csc);

csc_type = csc.type;
csc_tvec = csc.tvec;
csc_data = csc.data;
csc_label = csc.label;

save('C:\Users\Emily\Dropbox\Graduate courses\psyc-179\inputs_csc', ...
    'csc_data', 'csc_tvec', 'csc_type', 'csc_label')

%% input_positions
cfg_pos = [];
position = LoadPos(cfg_pos);

pos_type = position.type;
pos_tvec = position.tvec;
pos_datax = position.data(1,:);
pos_datay = position.data(2,:);
pos_label = position.label;

save('C:\Users\Emily\Dropbox\Graduate courses\psyc-179\inputs_position', ...
    'pos_datax', 'pos_datay', 'pos_tvec', 'pos_type', 'pos_label')

%% input_events
cfg_evt = [];
cfg_evt.eventList = {'TTL Output on AcqSystem1_0 board 0 port 0 value (0x0004).','TTL Output on AcqSystem1_0 board 0 port 0 value (0x0040).'};
cfg_evt.eventLabel = {'FoodDelivery','WaterDelivery'};

evt = LoadEvents(cfg_evt);

evt_type = evt.type;
evt_food = evt.t{1};
evt_water = evt.t{2};
evt_label = evt.label;

save('C:\Users\Emily\Dropbox\Graduate courses\psyc-179\inputs_event', ...
    'evt_food', 'evt_water', 'evt_type', 'evt_label')

%% input_spikes
cfg_spk = [];

spikes = LoadSpikes(cfg_spk);

spikes_type = spikes.type;
spikes_times = spikes.t;
spikes_label = spikes.label;

save('C:\Users\Emily\Dropbox\Graduate courses\psyc-179\inputs_spike', ...
    'spikes_times', 'spikes_label', 'spikes_type')
