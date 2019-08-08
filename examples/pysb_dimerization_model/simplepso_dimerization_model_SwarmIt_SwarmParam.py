"""
PSO run script for dimerization_model.py using swarm_it
classes.
"""
from dimerization_model_swarmit import model, swarm_param
import numpy as np

from swarm_it import SwarmIt

if __name__ == '__main__':
    # Initialize PySB solver object for running simulations.
    # Simulation timespan should match experimental data.
    tspan = np.linspace(0,1, num=51)

    # USER must add commands to import/load any experimental
    # data for use in the likelihood function!
    experiments_avg = np.load('dimerization_model_dimer_data.npy')
    experiments_sd = np.load('dimerization_model_dimer_sd.npy')

    # Setup the PSO run
    observable_data = dict()
    time_idxs = list(range(len(tspan)))
    observable_data['A_dimer'] = tuple((experiments_avg, experiments_sd, time_idxs))
    # Initialize the GAlibrateIt instance with the model details.
    swarmit = SwarmIt(model, observable_data, tspan, swarm_param=swarm_param)
    # Now build the GAO object. -- All inputs are
    # optional keyword arguments.
    pso = swarmit(pso_kwargs=dict({'save_sampled':False, 'verbose':True}),
                  cost_type='norm_logpdf')

    # run it
    num_particles = 20
    num_iterations = 50
    pso.run(num_particles, num_iterations)
    print("best_theta: ",pso.best)
