import os
import numpy as np

from bmtk.builder.networks import NetworkBuilder


# Step 1: Create a v1 mock network of 14 cells (nodes) with across 7 different cell "types"
net = NetworkBuilder("v1")
net.add_nodes(N=1,  # specifiy the number of cells belong to said group.
              pop_name='', location='L2', ei='e',  # pop_name, location, and ei are optional parameters that help's identifies properties of the cells. The modeler can choose whatever key-value pairs as they deem appropiate.
              positions=[(2.753, -464.868, -161.705)],  # The following properties we are passing in lists
              tuning_angle=[0.0],                 #  values to each individual cell
	      rotation_angle_xaxis=-0.010516,
              rotation_angle_yaxis=0,
              rotation_angle_zaxis=-3.011895, # Note that the y-axis rotation is differnt for each cell (ie. given a list of size N), but with z-axis rotation all cells have the same value
              model_type='biophysical',  # The type of cell we are using
              model_template='ctdb:Biophys1.hoc',  # Tells the simulator that when building cells models use a hoc_template specially created for parsing Allen Cell-types file models. Value would be different if we were using NeuronML or different model files
              model_processing='aibs_allactive_ani_directed',  # further instructions for how to processes a cell model. In this case aibs_perisomatic is a built-in directive to cut the axon in a specific way
              dynamics_params='optim_param_571654895.json',  # Name of file (downloaded from Allen Cell-Types) used to set model parameters and channels
              morphology='571654895.swc'),  # Name of morphology file downloaded


# Step 2: We want to connect our network. Just like how we have node-types concept we group our connections into
# "edge-types" that share rules and properties
net.add_edges(source={'ei': 'e'}, target={'pop_name': ''},
              connection_rule=5,
              syn_weight=4e-05,
              weight_function='gaussianLL',
              weight_sigma=50.0,
              distance_range=[30.0, 150.0],
              target_sections=['basal', 'apical'],
              delay=2.0,
              dynamics_params='AMPA_ExcToExc.json',
              model_template='exp2syn')
net.build()
net.save(output_dir='network')


def generate_positions(N, x0=0.0, x1=300.0, y0=0.0, y1=100.0):
    X = np.random.uniform(x0, x1, N)
    Y = np.random.uniform(y0, y1, N)
    return np.column_stack((X, Y))


def select_source_cells(src_cells, trg_cell, n_syns):
    synapses = [n_syns]*len(src_cells)
    
    #if 'tuning_angle' in trg_cell:
    #    synapses = [n_syns if src['pop_name'] == 'tON' or src['pop_name'] == 'tOFF' else 0 for src in src_cells]
    #else:
    #    synapses = [n_syns if src['pop_name'] == 'tONOFF' else 0 for src in src_cells]
    
    return synapses


lgn = NetworkBuilder("lgn")
lgn.add_nodes(N=15, pop_name='tON', ei='e', location='LGN',
              positions=generate_positions(15),
              model_type='virtual')

lgn.add_nodes(N=15, pop_name='tOFF', ei='e', location='LGN',
              positions=generate_positions(15),
              model_type='virtual')

lgn.add_nodes(N=15, pop_name='tONOFF', ei='e', location='LGN',
              positions=generate_positions(15),
              model_type='virtual')


lgn.add_edges(source=lgn.nodes(), target=net.nodes(pop_name=''),
              iterator='all_to_one',
              connection_rule=select_source_cells,
              connection_params={'n_syns': 10},
              syn_weight=4e-05,
              weight_function='wmax',
              distance_range=[0.0, 150.0],
              target_sections=['basal', 'apical'],
              delay=2.0,
              dynamics_params='AMPA_ExcToExc.json',
              model_template='exp2syn')
lgn.build()
lgn.save(output_dir='network')


tw = NetworkBuilder("tw")
tw.add_nodes(N=15, pop_name='TW', ei='e', location='TW', model_type='virtual')

tw.add_edges(source=tw.nodes(), target=net.nodes(pop_name='H16'),
             connection_rule=5,
             syn_weight=0.0000,
             weight_function='wmax',
             distance_range=[30.0, 150.0],
             target_sections=['basal', 'apical'],
             delay=2.0,
             dynamics_params='AMPA_ExcToExc.json',
             model_template='exp2syn')

tw.build()
tw.save(output_dir='network')
