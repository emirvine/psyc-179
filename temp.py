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
        water_spikes['time'].append(time_slice(spikes['time'][neuron], water_starts[0],
                                               water_stops[0]))

print str(len(spikes['time'])) + ' should equal ' + str(len(water_spikes['time']))

for neuron in range(len(water_spikes['time'])):
        plt.plot(water_spikes['time'][neuron], np.zeros(len(water_spikes['time'][neuron]))+neuron+1,
                 '|', color='k')
plt.show()






