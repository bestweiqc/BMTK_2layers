from bmtk.builder import NetworkBuilder
import numpy as np
from bmtk.builder.auxi.node_params import positions_cuboid, positions_list
import math
import random
np.random.seed(91)

# Initialize our network
netff = NetworkBuilder("ff_model")

# Create the possible x,y,z coordinates
xside_length = 600; yside_length = 600; height = 60; min_dist = 25;
x_grid = np.arange(0,xside_length+min_dist,min_dist)
y_grid = np.arange(0,yside_length+min_dist,min_dist)
z_grid = np.arange(0,height+min_dist,min_dist)
xx, yy, zz = np.meshgrid(x_grid, y_grid, z_grid)
pos_list = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T

# Number of cells in each population

numBask = 20
numL5PNC = 40
numL5PNA = 40

inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numL5PNA,replace=False)
pos = pos_list[inds,:]

netff.add_nodes(N=numL5PNA, pop_name='numL5PNA',
                positions=positions_list(positions=pos),
                potental='exc',
                model_type='biophysical',
                model_template='hoc:Cell_CP',
                morphology=None
                )


inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numL5PNC,replace=False)
pos = pos_list[inds,:]

netff.add_nodes(N=numL5PNC, pop_name='L5PNC',
                positions=positions_list(positions=pos),
                potental='exc',
                model_type='biophysical',
                model_template='hoc:Cell_CS',
                morphology=None
                )

############## ###################################################################
###########################Fast - spiking PV ints################################

# Get rid of coordinates already used
pos_list = np.delete(pos_list,inds,0)

# Pick new coordinates
inds = np.random.choice(np.arange(0,np.size(pos_list,0)),numBask,replace=False)
pos = pos_list[inds,:]

# Add a population of numBask nodes
netff.add_nodes(N=numBask, pop_name='Cell_Bask',
                positions=positions_list(positions=pos),
                potental='i',
                model_type='biophysical',
                model_template='hoc:Cell_PN',
                morphology=None
                )


##############################################################################
############################## CONNECT CELLS #################################

def dist_conn_perc(source, target, prob=0.1, min_dist=0.0, max_dist=300.0, min_syns=1, max_syns=2):
    sid = source.node_id
    tid = target.node_id
    # No autapses
    if sid==tid:
        return None
    else:
        src_pos = source['positions']
        trg_pos = target['positions']
    dist =np.sqrt((src_pos[0]-trg_pos[0])**2+(src_pos[1]-trg_pos[1])**2+(src_pos[2]-trg_pos[2])**2)
        #print("src_pos: {} trg_pos: {} dist: {}".format(src_pos,trg_pos,dist))        

    if dist <= max_dist and np.random.uniform() < prob:
        tmp_nsyn = np.random.randint(min_syns, max_syns)
        print("{}to{}done".format(sid,tid))
        #print("creating {} synapse(s) between cell {} and {}".format(tmp_nsyn,sid,tid))
    else:
        tmp_nsyn = 0

    return tmp_nsyn

def dist_conn_perc1(source, target,prob=0 ,min_dist=0, max_dist=600, min_syns=1, max_syns=1):
	
        x_ind,y_ind,z_ind = 0,1,2
        dx = target['positions'][x_ind] - source['positions'][x_ind]
        dy = target['positions'][y_ind] - source['positions'][y_ind]
        dz = target['positions'][z_ind] - source['positions'][z_ind]
        #distxyz = math.sqrt(dx**2 + dy**2)
        distxyz = math.sqrt(dx**2 + dy**2 + dz**2)
	#print("conn one")
        prob = 1

        if distxyz <= 50:
            prob = 0.03
            
        elif distxyz <= 100:
            prob = 0.02
        elif distxyz <= 200:
            prob = 0.01

        else:
            prob = 0.005
        
        if random.random() < prob:
            #Since there will be recurrect connections we need to keep track externally to BMTK
            #BMTK will call build_edges twice if we use net.edges() before net.build()
            #Resulting in double the edge count
            #syn_list.append({'source_gid':source['node_id'],'target_gid':target['node_id']})
            tmp_nsyn= 1 ## random.randint(min_syns,max_syns)
            print("{}to{}done".format(source.node_id,target.node_id))
        else:
            tmp_nsyn=0
        #print('PyrA, PyrC connection=',tmp_nsyn)
        return tmp_nsyn


# Create connections between Pyr --> Bask cells


netff.add_edges(source={'pop_name': ['L5PNA','L5PNC']}, target={'pop_name': ['Cell_Bask']},
                connection_rule=dist_conn_perc,
                connection_params={'prob':0.12,'min_dist':0.0,'max_dist':300.0,'min_syns':1,'max_syns':2},
                syn_weight=20,
                #weight_function = 'Lognormal',
                #weight_sigma=2,
                dynamics_params='AMPA_ExcToInh.json',
                model_template='exp2syn',
                distance_range=[0.0, 300.0],
                target_sections=['basal', 'apical'],
                delay=2.0)


netff.add_edges(source={'pop_name': ['L5PNA','L5PNC']}, target={'pop_name': ['L5PNA','L5PNC']},
                connection_rule=dist_conn_perc,
                connection_params={'prob':0.02,'min_dist':0.0,'max_dist':300.0,'min_syns':1,'max_syns':2},
                syn_weight=40,
                #weight_function = 'Lognormal',
                #weight_sigma=5,
                dynamics_params='AMPA_ExcToInh.json',
                model_template='exp2syn',
                distance_range=[0.0, 300.0],
                target_sections=['basal', 'apical'],
                delay=2.0)

# Create connections between Bask --> Bask cellss
netff.add_edges(source={'pop_name': 'Cell_Bask'}, target={'pop_name': ['Cell_Bask']},
                connection_rule=dist_conn_perc,
                connection_params={'prob':0.26,'min_dist':0.0,'max_dist':600.0,'min_syns':1,'max_syns':2},
                syn_weight=5,
                #weight_function = 'Lognormal',
                #weight_sigma=2,
                dynamics_params='GABA_InhToExc.json',
                model_template='exp2syn',
                distance_range=[0.0, 600.0],
                target_sections=['basal', 'apical'],
                delay=2.0)
# Create connections between Bask --> Pyr cellss
netff.add_edges(source={'pop_name': 'Cell_Bask'}, target={'pop_name': ['L5PNA','L5PNC']},
                connection_rule=dist_conn_perc,
                connection_params={'prob':0.34,'min_dist':0.0,'max_dist':300.0,'min_syns':1,'max_syns':2},
                syn_weight=10,
                #weight_function = 'Lognormal',
                #weight_sigma=5,
                dynamics_params='GABA_InhToExc.json',
                model_template='exp2syn',
                distance_range=[0.0, 300.0],
                target_sections=['basal', 'apical'],
                delay=2.0)

netff.build()
netff.save_nodes(output_dir='network')
netff.save_edges(output_dir='network')

print("Internal nodes and edges built")

# Create connections between "thalamus" and Pyramidals
# First define the connection rule
def one_to_one(source, target):
    #print("one to one")
    sid = source.node_id
    tid = target.node_id
    
    if tid > 89 :
        print("working on Bask")
    if sid == tid:
        print("connecting cell {} to {}".format(sid,tid))
        tmp_nsyn = 1
    else:
        tmp_nsyn = 0
    return tmp_nsyn
################################################################################
############################# BACKGROUND INPUTS ################################

# External inputs
thalamus = NetworkBuilder('mthalamus')
thalamus.add_nodes(N=numL5PNA+numL5PNC+numBask,
                   pop_name='tON',
                   potential='exc',
                   model_type='virtual')
print(thalamus.nodes)

thalamus.add_edges(source=thalamus.nodes(), target=netff.nodes(pop_name=['L5PNA','L5PNC']),
                   connection_rule=one_to_one,
                   syn_weight=70,
                   delay=2.0,
                   weight_function=None,
                   target_sections=['basal', 'apical'],
                   distance_range=[0.0, 300.0],
                   dynamics_params='AMPA_ExcToExc.json',
                   model_template='exp2syn')

thalamus.add_edges(source=thalamus.nodes(), target=netff.nodes(pop_name='Cell_Bask'),
                   connection_rule=one_to_one,
                   syn_weight=5,
                   target_sections=['somatic'],
                   delay=2.0,
                   distance_range=[0.0, 600.0],
                   dynamics_params='AMPA_ExcToInh.json',
                   model_template='exp2syn')
thalamus.build()
thalamus.save_nodes(output_dir='network')
thalamus.save_edges(output_dir='network')

print("External nodes and edges built")

# Build and save our network

from bmtk.utils.sim_setup import build_env_bionet

build_env_bionet(base_dir='./',
		network_dir='./network',
		tstop=1000.0, dt = 0.1,
		spikes_inputs=[('mthalamus',   # Name of population which spikes will be generated for
                        'mthalamus_spikes.h5')],
		report_vars=['v'],     # Record membrane potential and calcium (default soma)
		components_dir='biophys_components',
		compile_mechanisms=True)

from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator

psg = PoissonSpikeGenerator(population='mthalamus')
psg.add(node_ids=range(numL5PNA+numL5PNC+numBask),  # Have nodes to match mthalamus
        firing_rate=15.0,    # 15 Hz, we can also pass in a nonhomoegenous function/array
        times=(0.0, 3.0))    # Firing starts at 0 s up to 3 s
psg.to_sonata('mthalamus_spikes.h5')
