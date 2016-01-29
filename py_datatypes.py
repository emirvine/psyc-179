import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio


def import_csc(matfile):
    load_csc = sio.loadmat(matfile)
    csc = dict(time=[])
    csc['data'] = load_csc['csc_data'][0]
    for val in range(len(load_csc['csc_tvec'])):
        csc['time'].append(load_csc['csc_tvec'][val][0])
    csc['type'] = load_csc['csc_type'][0]
    csc['label'] = load_csc['csc_label'][0][0][0]
    return csc


def import_position(matfile):
    load_pos = sio.loadmat(matfile)
    pos = dict()
    pos['x'] = load_pos['pos_datax'][0]
    pos['y'] = load_pos['pos_datay'][0]
    pos['time'] = load_pos['pos_tvec'][0]
    pos['type'] = load_pos['pos_type'][0]
    pos['label'] = load_pos['pos_label'][0][0][0]
    return pos


def import_events(matfile):
    load_events = sio.loadmat(matfile)
    events = dict()
    events['food'] = load_events['evt_food'][0]
    events['water'] = load_events['evt_water'][0]
    events['type'] = load_events['evt_type'][0]
    events['label'] = load_events['evt_label'][0][0][0]
    return events


def import_spikes(matfile):
    load_spikes = sio.loadmat(matfile)
    spikes = dict()
    spikes['time'] = load_spikes['spikes_times'][0]
    spikes['type'] = load_spikes['spikes_type'][0]
    spikes['label'] = load_spikes['spikes_label'][0][0][0]
    return spikes


csc = import_csc('inputs_csc.mat')
pos = import_position('inputs_position.mat')
events = import_events('inputs_event.mat')
spikes = import_spikes('inputs_spike.mat')

# Plotting lfp
# plt.plot(csc['time'], csc['data'], 'k')
# plt.xlim(3660, 3720)
# plt.show()

# Plotting position
# plt.plot(pos['x'], pos['y'], 'g.', ms=2)
# plt.axis('off')
# plt.show()

# Plotting event times
# plt.plot(events['food'], np.zeros(len(events['food'])), '|', color='g', ms=200)
# plt.plot(events['water'], np.zeros(len(events['water'])), '|', color='b', ms=200)
# plt.ylim(-0.1, 0.2)
# plt.show()

# Plotting spike times
# for neuron in range(len(spikes['time'])):
#     plt.plot(spikes['time'][neuron], np.ones(len(spikes['time'][neuron]))+neuron+1, '|', color='k')
# plt.xlim(3660, 3720)
# plt.ylabel('Neuron number')
# plt.xlabel('Time (ms?)')
# plt.title('Check it out! I can plot spikes!')
# plt.show()



