from bmtk.analyzer import nodes_table
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pdb
from mpl_toolkits.mplot3d import Axes3D
import h5py
from bmtk.analyzer.cell_vars import _get_cell_report, plot_report
import matplotlib.pyplot as plt
from scipy.signal import welch

# Load data
config_file = "simulation_config.json"
lfp_file = "./output/ecp.h5"
mem_pot_file = './output/v_report.h5'
raster_file = './output/spikes.h5'


# load 
f = h5py.File(mem_pot_file,'r')


mem_potential = f['report']['ff_model']['data']
f = h5py.File(raster_file,'r')
gids = f['spikes']['ff_model']['node_ids']
timestamps = f['spikes']['ff_model']['timestamps']
plt.figure(1)
plt.plot(mem_potential)
plt.figure(2)
plt.plot(timestamps,gids,'.')


plt.show()
