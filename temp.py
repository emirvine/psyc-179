# from shapely.geometry import Point, LineString
# import numpy as np
# import matplotlib.pyplot as plt
# import scipy.io as sio
# import scipy
# from py_datatypes import *

# T-maze ideal points
# path_pts = dict()
# path_pts['start_box'] = (237.1, 243.4)
# path_pts['choice_point'] = (596.5, 243.4)
#
# path_pts['food_turn1'] = (595.6, 372.4)
# path_pts['food_turn2'] = (585.7, 408.9)
# path_pts['food_turn3'] = (553.2, 444.0)
# path_pts['food_turn4'] = (498.2, 472.2)
# path_pts['food_reward'] = (448.8, 473.6)
# path_pts['food_pedestal'] = (348.6, 376.8)
#
# path_pts['water_turn1'] = (579.1, 105.5)
# path_pts['water_turn2'] = (568.7, 83.9)
# path_pts['water_turn3'] = (517.0, 45.3)
# path_pts['water_turn4'] = (492.5, 31.3)
# path_pts['water_reward'] = (452.1, 31.3)
# path_pts['water_pedestal'] = (348.6, 137.5)
#
# food_line = LineString([path_pts['start_box'], path_pts['choice_point'], path_pts['food_turn1'],
#                         path_pts['food_turn2'], path_pts['food_turn3'], path_pts['food_turn4'],
#                         path_pts['food_reward']])
#
# water_line = LineString([path_pts['start_box'], path_pts['choice_point'], path_pts['water_turn1'],
#                          path_pts['water_turn2'], path_pts['water_turn3'], path_pts['water_turn4'],
#                          path_pts['water_reward']])
#
# # start/stop times for trials from Alyssa's metadata
# water_starts = [3240.5, 3591.8, 3744.1, 3891.7, 4145.1, 4966.5, 5085.7, 5214.4, 5330.3]
# water_stops = [3282.1, 3605.4, 3754.9, 3905.5, 4170.3, 4982.1, 5106.4, 5232.3, 5357.6]
#
# food_starts = [3433.5, 4015.4, 4267.6, 4404.5, 4540.3, 4703.8, 4822.6, 5749.6, 5583.6]
# food_stops = [3448.2, 4044.4, 4284.5, 4420.4, 4583.4, 4718.8, 4870.3, 5491.3, 5622.4]

# Linear paths for food and water trials.
# plt.plot(pos['x'], pos['y'], 'y')
# linear_food = linear_trajectory(pos, food_line, food_starts, food_stops)
# linear_water = linear_trajectory(pos, water_line, water_starts[0], water_stops[0])
# plt.plot(linear_food['x'], linear_food['y'], 'k', ms=40)
# plt.plot(linear_water['x'], linear_water['y'], 'g', ms=40)
# plt.show()

# print linear_water
#
# bin_size = 100
# dimension = len(linear_water)
#
# for edge in range(dimension):
#     bin_start = min(linear_water[edge, :])
#     bin_stop = max(linear_water[edge, :])
#     bin_edges.append(np.linspace(bin_start, bin_stop, num=bin_size))
#
# dt = 1/30.
# min_occupancy = 1
#
# pos_idx = np.searchsorted(linear_water, bin_edges)




# swr = scipy.signal.hilbert(csc['data'])
# print len(swr)

# Notes:
# - function to make tuning curve from spikes/times
# - tests? eg. time sliced spikes should have same neuron number as non-sliced

#################################

# Alyssa's Hilbert transform call in matlab
# case 'HT'
#         cfg.stepSize = [];
#         cfg.weightby = [];
#         cfg_temp = []; cfg_temp.verbose = cfg.verbose; cfg_temp.rippleband = [140 250]; cfg_temp.smooth = 1; cfg_temp.kernel = [];
#         SWR = OldWizard(cfg_temp,CSC);
#
# function SWR = OldWizard(cfg_in,CSC)
# %% Parse cfg parameters
# cfg_def.rippleband = [140 250]; % in Hz
# cfg_def.smooth = 1; % do you want to smooth the detector or not
# cfg_def.kernel = []; % which kernel (if empty, goes to default)
# cfg_def.verbose = 1; % talk to me or not
#
# mfun = mfilename;
# cfg = ProcessConfig(cfg_def,cfg_in,mfun);
#
# if cfg.verbose
#     tic
#     disp([mfun,': looking for sharp wave-ripple events...'])
# end
#
# % filter in the ripple band
#
# cfg_temp = [];
# cfg_temp.type = 'fdesign';
# cfg_temp.f = cfg.rippleband;
# cfg_temp.verbose = 0;
# CSCf = FilterLFP(cfg_temp,CSC);
#
# % ask hilbert what he thinks
# score = abs(hilbert(CSCf.data)).^2;
#
# % apply smoothing, if desired
#
# if cfg.smooth
#     if isempty(cfg.kernel)
#         kernel = gausskernel(60,20);
#     else
#         kernel = cfg.kernel;
#     end
#    score = conv(score,kernel,'same');
# end
#
# % make output
# SWR = tsd(CSC.tvec,score);
#
# water_starts = [3240.5, 3591.8, 3744.1, 3891.7, 4145.1, 4966.5, 5085.7, 5214.4, 5330.3]
# water_stops = [3282.1, 3605.4, 3754.9, 3905.5, 4170.3, 4982.1, 5106.4, 5232.3, 5357.6]
#
# trial_start = water_starts[1]
# trial_stop = water_stops[1]
#
# water_spikes = dict(time=[])
# for neuron in range(len(spikes['time'])):
#         water_spikes['time'].append(time_slice(spikes['time'][neuron], trial_start, trial_stop))
#         plt.plot(spikes['time'][neuron], np.zeros(len(spikes['time'][neuron]))+neuron+1, '|', color='r')
#
# print str(len(spikes['time'])) + ' should equal ' + str(len(water_spikes['time']))
#
# for neuron in range(len(water_spikes['time'])):
#     pass
#     #plt.plot(water_spikes['time'][neuron], np.zeros(len(water_spikes['time'][neuron]))+neuron+1, '|', color='k')
# #plt.xlim(trial_start, trial_stop)
# plt.show()
#
# # test = [30]
# # plt.plot(water_spikes['time'][test[0]], np.ones(len(water_spikes['time'][test[0]])), '|', color='k', ms=20)
# firing_rate = dict(time=[])
# bin_size = 5
# bins = []
#
# for neuron in range(len(water_spikes['time'])):
#     bin_start = trial_start
#     num_spikes = []
#     while bin_start < trial_stop:
#         bin_end = bin_start + bin_size
#         bin_stop = min(bin_end, trial_stop)
#         binned_neurons = time_slice(water_spikes['time'][neuron], bin_start, bin_stop)
#         bins.append(bin_stop)
#         bin_start = bin_stop + 0.0001
#         num_spikes.append(len(binned_neurons))
#     firing_rate['time'].append(num_spikes)
# plt.xlim(trial_start, trial_stop)
#
# # print len(bins)
# for neuron in range(len(water_spikes['time'])):
# # print firing_rate['time']
#     plt.plot(range(len(firing_rate['time'][neuron])), firing_rate['time'][neuron], color='k')
# plt.show()
#
#
#
