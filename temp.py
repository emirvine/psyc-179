from shapely.geometry import Point, LineString
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy
from py_datatypes import *

water_starts = [3240.5, 3591.8, 3744.1, 3891.7, 4145.1, 4966.5, 5085.7, 5214.4, 5330.3]
water_stops = [3282.1, 3605.4, 3754.9, 3905.5, 4170.3, 4982.1, 5106.4, 5232.3, 5357.6]

water_spikes = dict(time=[])
for neuron in range(len(spikes['time'])):
        water_spikes['time'].append(time_slice(spikes['time'][neuron], water_starts[0], water_stops[0]))

print str(len(spikes['time'])) + ' should equal ' + str(len(water_spikes['time']))

# for neuron in range(len(water_spikes['time'])):
#         plt.plot(water_spikes['time'][neuron], np.zeros(len(water_spikes['time'][neuron]))+neuron+1,
#                  '|', color='k')
# plt.show()

test = 30
plt.plot(water_spikes['time'][test], np.ones(len(water_spikes['time'][test])), '|', color='k', ms=20)


firing_rate = dict(time=[])
bin_size = 5
for neuron in range(len(water_spikes['time'])):
    neuron_firing = []
    bin_start = water_starts[0]

    while bin_start < water_stops[0]:
        bin_end = bin_start + bin_size
        bin_stop = min(bin_end, water_stops[0])
        binned_neurons = time_slice(water_spikes['time'][test], bin_start, bin_stop)
        plt.plot(binned_neurons, np.ones(len(binned_neurons)), '|', color='y', ms=10)

        bin_start = bin_stop + 0.0001
        print len(binned_neurons)
        neuron_firing.append(len(binned_neurons))
    firing_rate['time'].append(neuron_firing)

plt.show()

# for neuron in range(len(water_spikes['time'])):
print firing_rate['time'][test]
plt.plot(range(len(firing_rate['time'][test])), firing_rate['time'][test], color='k')
plt.show()




