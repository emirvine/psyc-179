import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy
from shapely.geometry import Point, LineString


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
    events['led1'] = load_events['evt_led1id'][0]
    events['led2'] = load_events['evt_led2id'][0]
    events['ledoff'] = load_events['evt_ledoff'][0]
    events['pb1'] = load_events['evt_pb1id'][0]
    events['pb2'] = load_events['evt_pb2id'][0]
    events['pboff'] = load_events['evt_pboff'][0]
    events['feeder1'] = load_events['evt_feeder1id'][0]
    events['feeder2'] = load_events['evt_feeder2id'][0]
    events['feederoff'] = load_events['evt_feederoff'][0]
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


def find_nearest_idx(array, val):
    return (np.abs(array-val)).argmin()


def time_slice(spikes, t_start, t_stop):
    """
    Creates a new spike train sliced to the time interval given by start
    and end times (inclusive). Start or end times can also be None to use
    infinite endpoints for the time interval.
    :param spikes: spike train to be modified
    :param t_start: start interval
    :param t_stop: end interval
    :return: sliced_spikes: based on interval times
    """
    if t_start is None:
        t_start = -np.inf
    if t_stop is None:
        t_stop = np.inf
    indices = (spikes >= t_start) & (spikes <= t_stop)
    sliced_spikes = spikes[indices]
    return sliced_spikes


def linear_trajectory(pos, ideal_path, trial_start, trial_stop):
    for trial in range(len(trial_start)):
        t_start_idx = find_nearest_idx(np.array(pos['time']), trial_start[trial])
        t_end_idx = find_nearest_idx(np.array(pos['time']), trial_stop[trial])

        pos_trial = dict()
        pos_trial['x'] = pos['x'][t_start_idx:t_end_idx]
        pos_trial['y'] = pos['y'][t_start_idx:t_end_idx]
        pos_trial['time'] = pos['time'][t_start_idx:t_end_idx]

        linear_pos = dict(x=[], y=[])
        linear_pos['time'] = pos_trial['time']
        for point in range(len(pos_trial['x'])-1):
            position = Point(pos_trial['x'][point], pos_trial['y'][point])
            linearized_point = ideal_path.interpolate(ideal_path.project(position))
            linear_pos['x'].append(linearized_point.xy[0])
            linear_pos['y'].append(linearized_point.xy[1])
    return linear_pos

csc = import_csc('emi_inputs_csc.mat')
# pos = import_position('emi_inputs_position.mat')
# events = import_events('emi_inputs_event.mat')
# spikes = import_spikes('emi_inputs_spike.mat')


# Plotting lfp
plt.plot(csc['time'], csc['data'], 'k')
plt.show()

# Plotting event times
# plt.plot(events['food'], np.zeros(len(events['food'])), '|', color='g', ms=200)
# plt.plot(events['water'], np.zeros(len(events['water'])), '|', color='b', ms=200)
# plt.ylim(-0.1, 0.2)
# plt.show()

# event = events['food'][0]
# t_start = event - 1
# t_stop = event + 1

# sliced_spikes = dict(time=[])
# for neuron in range(len(spikes['time'])):
#     sliced_spikes['time'].append(time_slice(spikes['time'][neuron], t_start, t_stop))
#
# print str(len(spikes['time'])) + ' should equal ' + str(len(sliced_spikes['time']))
#
# for neuron in range(len(spikes['time'])):
#     plt.plot(sliced_spikes['time'][neuron], np.ones(len(sliced_spikes['time'][neuron]))+neuron+1,
#              '|', color='k')
# plt.plot(events['food'][0], np.zeros(np.array(events['food'][0]).size)+(len(sliced_spikes['time'])/2.), '|', color='r', ms=300)
# plt.ylabel('Neuron number')
# plt.xlabel('Time (ms?)')
# plt.title('Check it out!!! I can slice spikes!')
# plt.show()

# t_start_idx = find_nearest_idx(np.array(csc['time']), t_start)
# t_end_idx = find_nearest_idx(np.array(csc['time']), t_stop)
#
# sliced_csc = dict()
# sliced_csc['data'] = csc['data'][t_start_idx:t_end_idx]
# sliced_csc['time'] = csc['time'][t_start_idx:t_end_idx]
#
# plt.plot(sliced_csc['time'], sliced_csc['data'], 'y')
# plt.plot(events['food'][0], np.zeros(1), '|', color='k', ms=300)
# plt.show()


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