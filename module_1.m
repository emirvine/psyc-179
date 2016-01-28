cd('C:\Users\Emily\Desktop\R016-2012-10-08_promoted');

cfg = [];
cfg.fc = {'R016-2012-10-08-CSC02d.ncs'};
csc = LoadCSC(cfg);

plot(csc.tvec, csc.data);
xlim([1338.6,1339.2]);