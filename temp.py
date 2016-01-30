from shapely.geometry import Point, LineString
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import scipy
from py_datatypes import *

csc = import_csc('inputs_csc.mat')
pos = import_position('inputs_position.mat')
events = import_events('inputs_event.mat')
spikes = import_spikes('inputs_spike.mat')

path_pts = dict()
path_pts['start_box'] = (237.1, 243.4)
path_pts['choice_point'] = (596.5, 243.4)

path_pts['left_turn1'] = (595.6, 372.4)
path_pts['left_turn2'] = (585.7, 408.9)
path_pts['left_turn3'] = (553.2, 444.0)
path_pts['left_turn4'] = (498.2, 472.2)
path_pts['left_reward'] = (448.8, 473.6)
path_pts['left_pedestal'] = (348.6, 376.8)

path_pts['right_turn1'] = (579.1, 105.5)
path_pts['right_turn2'] = (568.7, 83.9)
path_pts['right_turn3'] = (517.0, 45.3)
path_pts['right_turn4'] = (492.5, 31.3)
path_pts['right_reward'] = (452.1, 31.3)
path_pts['right_pedestal'] = (348.6, 137.5)

# Plotting position
plt.plot(pos['x'], pos['y'], 'g.', ms=2)
plt.plot(path_pts['start_box'][0], path_pts['start_box'][1], 'b.', ms=12)
plt.plot(path_pts['choice_point'][0], path_pts['choice_point'][1], 'b.', ms=12)
plt.plot(path_pts['left_pedestal'][0], path_pts['left_pedestal'][1], 'b.', ms=12)
plt.plot(path_pts['right_pedestal'][0], path_pts['right_pedestal'][1], 'b.', ms=12)
plt.plot(path_pts['left_reward'][0], path_pts['left_reward'][1], 'b.', ms=12)
plt.plot(path_pts['right_reward'][0], path_pts['right_reward'][1], 'b.', ms=12)
plt.plot(path_pts['left_turn1'][0], path_pts['left_turn1'][1], 'r.', ms=12)
plt.plot(path_pts['left_turn2'][0], path_pts['left_turn2'][1], 'r.', ms=12)
plt.plot(path_pts['left_turn3'][0], path_pts['left_turn3'][1], 'r.', ms=12)
plt.plot(path_pts['left_turn4'][0], path_pts['left_turn4'][1], 'r.', ms=12)
plt.plot(path_pts['right_turn1'][0], path_pts['right_turn1'][1], 'r.', ms=12)
plt.plot(path_pts['right_turn2'][0], path_pts['right_turn2'][1], 'r.', ms=12)
plt.plot(path_pts['right_turn3'][0], path_pts['right_turn3'][1], 'r.', ms=12)
plt.plot(path_pts['right_turn4'][0], path_pts['right_turn4'][1], 'r.', ms=12)
plt.axis('off')


left_line = LineString([path_pts['start_box'], path_pts['choice_point'], path_pts['left_turn1'],
                        path_pts['left_turn2'], path_pts['left_turn3'], path_pts['left_turn4'],
                        path_pts['left_reward']])

right_line = LineString([path_pts['start_box'], path_pts['choice_point'], path_pts['right_turn1'],
                         path_pts['right_turn2'], path_pts['right_turn3'], path_pts['right_turn4'],
                         path_pts['right_reward']])

plt.plot(left_line.xy[0], left_line.xy[1], 'c', ms=10)
plt.plot(right_line.xy[0], right_line.xy[1], 'y', ms=10)

left_interpolate = left_line.interpolate(1.5, normalized=True)

plt.plot(left_interpolate.xy[0], left_interpolate.xy[1], 'k', ms=10)
plt.show()

