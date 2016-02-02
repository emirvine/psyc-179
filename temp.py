# from shapely.geometry import Point, LineString
# import numpy as np
# import matplotlib.pyplot as plt
# import scipy.io as sio
# import scipy
# from py_datatypes import *
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
