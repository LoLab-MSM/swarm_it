# swarm_it

![Python version badge](https://img.shields.io/badge/python-3.6-blue.svg)
[![license](https://img.shields.io/github/license/LoLab-VU/swarm_it.svg)](LICENSE)
![version](https://img.shields.io/badge/version-0.2.0-orange.svg)
[![release](https://img.shields.io/github/release-pre/LoLab-VU/swarm_it.svg)](https://github.com/LoLab-VU/swarm_it/releases/tag/v0.2.0)

**swarm_it** is a python utility designed to abstract away some of the effort in setting up Particle Swarm Optimization (PSO)-based calibrations of biological models in the [PySB](http://pysb.org/) format using [simplePSO](https://github.com/LoLab-VU/ParticleSwarmOptimization).

------

# Install

| **! Warning** |
| :--- |
|  swarm_it is still under heavy development and may rapidly change. |

**swarm_it** installs as the `swarm_it` package. It is compatible (i.e., tested) with Python 3.6.

Note that `swarm_it` has the following core dependencies:
   * [NumPy](http://www.numpy.org/)
   * [SciPy](https://www.scipy.org/)   
   * [simplePSO](https://github.com/LoLab-VU/ParticleSwarmOptimization) version 1.0
   * [PySB](http://pysb.org/)

### pip install
You can install the latest release of the `swarm_it` package using `pip` sourced from the GitHub repo:
```
pip install -e git+https://github.com/LoLab-VU/swarm_it@v0.1.0#egg=swarm_it
```
However, this will not automatically install the core dependencies. You will have to do that separately:
```
pip install numpy scipy simplepso
```
See the [http://pysb.org/download](http://pysb.org/download) page for instructions on installing PySB.

------

# License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

------

# Documentation and Usage

### Quick Overview
swarm_it is a utility that helps generate a simplePSO run script or simplepso.pso.PSO objects for a PySB model.

#### Commmand line use
swarm_it can be used as a command line utility to generate a template Particle Swarm Optimization run script for a PySB model. swarm_it reads the model file, imports and pulls out all the kinetic parameters, and then writes out a simplepso run script for that model. Users can then edit edit the generated template script to load any data and modify the cost function.

Run swarm_it from the command line with following format:
```
python -m swarm_it model.py output_path
```      
where output_path is the directory/folder location where you want the generated script to be saved as "simplepso_model_name.py".

The command line version of swarm_it also has support for setting the calibrated model parameters via an instance of [SwarmParam](https://github.com/LoLab-VU/swarm_it/tree/master#swarmparam) defined along with the model.   


#### Progammatic use via the SwarmIt class
The swarm_it utility can be used progammatically via the SwarmIt
class. It's importable from swarm_it:
```python
from swarm_it import SwarmIt
```
The SwarmIt class can build an instance of a simplepso.pso.PSO object.  
 Here's a faux minimal example:
```python
from my_pysb_model import model as my_model
from swarm_it import SwarmIt
import numpy as np

timespan = np.linspace(0., 10., 10)
data = np.load('my_data.npy')
data_sd = np.load('my_data_sd.npy')
observable_data = dict()
time_idxs = list(range(len(timespan)))
# my_observable is the name of an Observable in the model.
observable_data['my_observable'] = (data, data_sd, time_idxs)
# Initialize the SwarmIt instance with the model details.
swarmit = SwarmIt(my_model, observable_data, timespan)
# Now build the PSO object. -- All inputs are
# optional keyword arguments.
pso = swarmit(pso_kwargs=dict({'save_sampled':False, 'verbose':True}),
              cost_type='norm_logpdf')
# Then you can run it.
num_particles = 20
num_iterations = 50
pso.run(num_particles, num_iterations)
print("best_theta: ",pso.best)
```

#### SwarmParam
By default SwarmIt constructs the PSO object to sample all of a model's kinetic rate parameters; it assumes that the priors are uniform with size 4 orders of magnitude and centered on the values defined in the model. However, the swarm_it module has a built-in helper class, SwarmParam, which can be used in conjunction of with SwarmIt class to set which model parameters are sampled and the corresponding sampling ranges. SwarmParam can be used at the level of PySB model definition to log which parameters to include in a Particle Swarm Optimization run, or it can be used after model definition at as long as it is defined before defining a SwarmIt object. It can be imported from swarm_it:
```python
from swarm_it import SwarmParam
```
It is passed at instantiation of the SwarmIt class via the `swarm_param` input parameter,
```python
swarmit = SwarmIt(model, observable_data, tspan, swarm_param=swarm_param)
```
which uses it to build the set of parameters to sample, their staring PSO position and bounds, as well as the parameter mask for use with built-in cost functions.

Note that if you flag a parameter for sampling without setting sampling bounds, SwarmParam will by default assign the parameter bounds centered on the parameter's value (as defined in the model) with a width of 4 orders of magnitude.

For a more in depth example usage of SwarmParam, see the [simplepso_dimerization_model_SwarmIt_SwarmParam example](./examples/pysb_dimerization_model/simplepso_dimerization_model_SwarmIt_SwarmParam.py) and corresponding model [](./examples/pysb_dimerization_model/dimerization_model_swarm_it.py).

##### Built-in cost functions
SwarmIt currently has three pre-defined cost functions with different estimators, specified with the keyword parameter cost_type:
```python
# Now build the PSO object.
pso = swarmit(cost_type='mse')
```
The options are
  * 'norm_logpdf'=> (default) Compute the cost as the negative sum of normal distribution logpdf's over each timecourse data points and each given model observable.
  * 'mse'=>Compute the cost using the mean squared error estimator over the timecourse data and each given model observable.
  * 'sse'=>Compute the cost using the sum of squared errors estimator over the timecourse data and each given observable.

Each of these functions computes the cost estimate using the timecourse output of a model simulation for each observable defined in the `observable_data` dictionary.

##### Custom cost function
If you want to use a different or more complicated cost function than the built-ins, the `SwarmIt` class object can also accept a custom cost function defined by the user with the parameters `cost_type='custom'` and `custom_cost`:
```python
# Now build the PSO object with a custom cost function.
def user_cost(model, simulation_output):
    # Do something and compute the cost for this model and its evaluation.
    return cost
pso = swarmit(cost_type='custom', custom_cost=user_cost)
```
Note that the custom cost function defined by the user, `user_cost`, should accept two parameters: the model itself and the simulation output (or trajectory). The user can then define how to compute the cost using these input objects or using other data/objects defined outside `user_cost`.

### Examples
Additional example scripts that show how to setup and launch PSO runs using **swarm_it** can be found under [examples](./examples).

------

# Contact

To report problems or bugs please open a
[GitHub Issue](https://github.com/LoLab-VU/swarm_it/issues). Additionally, any
comments, suggestions, or feature requests for **swarm_it** can also be submitted as a
[GitHub Issue](https://github.com/LoLab-VU/swarm_it/issues).

------

# Citing

If you use the **swarm_it** software in your research, please cite the GitHub repo.

Also, please cite the following references as appropriate for software used with/via **swarm_it**:

#### Packages from the SciPy ecosystem

This includes NumPy and SciPy for which references can be obtained from:
https://www.scipy.org/citing.html

#### simplePSO
The simplePSO reference can be exported from its Zenodo DOI entry:
[10.5281/zenodo.2612912](https://doi.org/10.5281/zenodo.2612912)

#### PySB
  1. Lopez, C. F., Muhlich, J. L., Bachman, J. A. & Sorger, P. K. Programming biological models in Python using PySB. Mol Syst Biol 9, (2013). doi:[10.1038/msb.2013.1](dx.doi.org/10.1038/msb.2013.1)
