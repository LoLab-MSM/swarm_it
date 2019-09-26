"""Tools to create PSO objects in order to optimize PySB model parameters.

This module defines the SwarmParam and SwarmIt classes. SwarmParam is used to flag/log
PySB model parameters for sampling during the optimization. With a SwarmParam instance
users can optionally set the loc and width of each parameter. SwarmIt is used
generate an instance of the PSO object for a PySB. It automatically pulls out
kinetic parameters, or can get the parameters from a SwarmParam instance, and can
automatically define a cost function for the model and a given timeseries
dataset.

"""

import importlib
import os.path
try:
    import pysb
    from pysb.simulator import ScipyOdeSimulator
except ImportError:
    pass
import numpy as np
from scipy.stats import norm
from simplepso.pso import PSO

def is_numbers(inputString):
    return all(char.isdigit() for char in inputString)

def parse_directive(directive, priors, no_sample, Keq_sample):
    words = directive.split()
    if words[1] == 'no-sample':
        if is_numbers(words[2]):
            par_idx = int(words[2])
            par = model.parameters[par_idx].name
        else:
            par = words[2]
        no_sample.append(par)
    elif words[1] == 'sample_Keq':
        if is_numbers(words[2]):
            par_idx = int(words[2])
            par_f = model.parameters[par_idx].name
        if is_numbers(words[3]):
            par_idx = int(words[3])
            par_r = model.parameters[par_idx].name
        if is_numbers(words[4]):
            par_idx = int(words[4])
            par_rep = model.parameters[par_idx].name
        no_sample.append(par_rep)
        Keq_sample.append((par_f, par_r, par_rep))
    return

def prune_no_samples(parameters, no_sample):
    pruned_pars = [parameter for parameter in parameters if parameter[0].name not in no_sample]

    return pruned_pars

def update_with_Keq_samples(parameters, Keq_sample):
    k_reps = [sample[2] for sample in Keq_sample]
    pruned_pars = [parameter for parameter in parameters if parameter[0].name not in k_reps]
    return pruned_pars

def write_uniform_param(p_name, p_val):
    line = "sp_{} = SampledParameter(\'{}\', loc=np.log10({})-1.0, width=2.0)\n".format(p_name, p_name, p_val)
    return line

class SwarmParam(object):
    """Container to store parameters and priors for sampling.
    This object is geared towards flagging PySB model parameters at the level of
    model definition for genetic algorithm optimization. However, it could be used outside
    model definition.

    Attributes:
        parms (dict of :obj:): A dictionary keyed to the parameter names. The
            values are the parameter search space meta-parameters as given in
            the call function.

    """

    def __init__(self):
        self.parms = dict()
        return

    def __call__(self, parameter, loc=None, width=None):
        """Add a parameter to the list.

        Args:
            parameter (:obj:pysb.Parameter): The parameter to be registered for
                genetic algorithm optimization.
            loc (float): The lower boundary for the parameter search space.
                Should be defined using log10 scale. Default: None
                If None, loc will be set to numpy.log10(parameter.value)-2.0.
            width (float): The width of the parameter search space. Should be
                defined using log10 scale. Default: None
                if None, width will be set to 4.0 (i.e., 4 orders of magnitude).
        """
        if loc is None:
            if width is None:
                loc = np.log10(parameter.value)-2.
            else:
                loc = np.log10(parameter.value) - (width/2.)
        if width is None:
            width = 4.
        self.parms[parameter.name] = tuple((loc, width))
        return parameter

    def __getitem__(self, key):
        return self.parms[key]

    def __setitem__(self, key, loc_width):
        self.parms[key] = loc_width

    def __delitem__(self, key):
        del self.parms[key]

    def __contains__(self, key):
        return (key in self.names)

    def __iadd__(self, parm):
        self.__call__(parm)
        return self

    def __isub__(self, parm):
        try:
            name = parm.name
            self.__delitem__(name)
        except TypeError:
            self.__delitem__(parm)
        return self

    def names(self):
        return list(self.parms.keys())

    def keys(self):
        return self.parms.keys()

    def mask(self, model_parameters):
        names = self.names()
        return [(parm.name in names) for parm in model_parameters]

    def locs(self):
        return np.array([self.parms[name][0] for name in self.parms.keys()])

    def widths(self):
        return np.array([self.parms[name][1] for name in self.parms.keys()])

    def centers(self):
        return self.locs()+self.widths()/2.

    def lower(self):
        return self.locs()

    def upper(self):
        return self.locs() + self.widths()

    def add_all_kinetic_params(self, pysb_model):
        for rule in pysb_model.rules:
            if rule.rate_forward:
                 self.__call__(rule.rate_forward)
            if rule.rate_reverse:
                 self.__call__(rule.rate_reverse)
        return

    def add_all_nonkinetic_params(self, pysb_model):
        kinetic_params = list()
        for rule in pysb_model.rules:
            if rule.rate_forward:
                 kinetic_params.append(rule.rate_forward)
            if rule.rate_reverse:
                 kinetic_params.append(rule.rate_reverse)
        for param in pysb_model.paramters:
            if param not in kinetic_params:
                self.__call__(param)
        return

    def add_by_name(self, pysb_model, name_or_list):
        if isinstance(name_or_list, (list, tuple)):
            for param in pysb_model.parameters:
                if param.name in name_or_list:
                    self.__call__(param)
        else:
            for param in pysb_model.parameters:
                if param.name == name_or_list:
                    self.__call__(param)
        return

class SwarmIt(object):
    """Create instances of PSO objects for PySB models.

    Args:
        model (pysb.Model): The instance of the PySB model that you want to run
            genetic algorithm optimization on.
        observable_data (dict of tuple): Defines the observable data to
            use when computing the cost function. It is a dictionary
            keyed to the model Observables (or species names) that the
            data corresponds to. Each element is a 3 item tuple of format:
            (:numpy.array:data, None or :numpy.array:data_standard_deviations,
            None or :list like:time_idxs or :list like:time_mask).
        timespan (numpy.array): The timespan for model simulations.
        solver (:obj:): The ODE solver to use when running model simulations.
            Defaults to pysb.simulator.ScipyOdeSimulator.
        solver_kwargs (dict): Dictionary of optional keyword arguments to
            pass to the solver when it is initialized. Defaults to dict().
        swarm_param (:obj:swarm_it.SwarmParam): An user-defined
            instance of the SwarmParam class with the data about the parameters to
            be sampled. If None, the default parameters to be sampled
            are the kinetic rate parameters with ranges of four orders of
            magnitude. Default: None

    Attributes:
        model
        observable_data
        timespan
        solver
        solver_kwargs

    """
    def __init__(self, model, observable_data, timespan,
                 solver=ScipyOdeSimulator,
                 solver_kwargs=None, swarm_param=None):
        """Inits the SwarmIt."""
        if solver_kwargs is None:
            solver_kwargs = dict()
        self.model = model
        self.observable_data = observable_data
        self.timespan = timespan
        self.solver = solver
        self.solver_kwargs = solver_kwargs
        # self.ns_version = None
        self._pso_kwargs = None
        self._like_data = dict()
        self._data = dict()
        self._data_mask = dict()
        for observable_key in observable_data.keys():
            self._like_data[observable_key] = norm(loc=observable_data[observable_key][0],
                                               scale=observable_data[observable_key][1])
            self._data[observable_key] = observable_data[observable_key][0]
            self._data_mask[observable_key] = observable_data[observable_key][2]
            # print(observable_data[observable_key][2])
            if observable_data[observable_key][2] is None:
                self._data_mask[observable_key] = range(len(self.timespan))
        self._model_solver = solver(self.model, tspan=self.timespan, **solver_kwargs)
        if swarm_param is not None:
            parm_mask = swarm_param.mask(model.parameters)
            self._starting_position = swarm_param.centers()
            self._lower = swarm_param.lower()
            self._upper = swarm_param.upper()
            self._rate_mask = parm_mask
        else:
            swarm_param = SwarmParam()
            for rule in model.rules:
                if rule.rate_forward:
                     swarm_param(rule.rate_forward)
                if rule.rate_reverse:
                     swarm_param(rule.rate_reverse)
            parm_mask = swarm_param.mask(model.parameters)
#            self._sampled_parameters = [SampledParameter(parm.name, *swarm_param[parm.name]) for i,parm in enumerate(model.parameters) if parm_mask[i]]
            self._rate_mask = parm_mask
            self._starting_position = swarm_param.centers()
            self._lower = swarm_param.lower()
            self._upper = swarm_param.upper()
        self._param_values = np.array([param.value for param in model.parameters])
        return

    def norm_logpdf_cost(self, position):
        """Compute the cost using the normal distribution estimator.

        Args:
            position (numpy.array): The parameter vector the compute cost
                of.

        Returns:
            float: The natural logarithm of the likelihood estimate.

        """
        Y = np.copy(position)
        params = self._param_values.copy()
        params[self._rate_mask] = 10**Y
        sim = self._model_solver.run(param_values=[params]).all
        logl = 0.
        for observable in self._like_data.keys():
            sim_vals = sim[observable][self._data_mask[observable]]
            logl += np.sum(self._like_data[observable].logpdf(sim_vals))
        if np.isnan(logl):
            return np.inf
        return -logl,

    def mse_cost(self, position):
        """Compute the cost using the negative mean squared error estimator.

        Args:
            position (numpy.array): The parameter vector the compute cost of.

        Returns:
            float: The natural logarithm of the likelihood estimate.

        """
        Y = np.copy(position)
        params = self._param_values.copy()
        params[self._rate_mask] = 10**Y
        sim = self._model_solver.run(param_values=[params]).all
        logl = 0.0
        for observable in self._like_data.keys():
            sim_vals = sim[observable][self._data_mask[observable]]
            logl += np.mean((self._data[observable]-sim_vals)**2)
        if np.isnan(logl):
            return np.inf
        return logl,

    def sse_cost(self, position):
        """Compute the cost using the negative sum of squared errors estimator.

        Args:
            position (numpy.array): The parameter vector the compute cost
                of.

        Returns:
            float: The natural logarithm of the likelihood estimate.

        """
        Y = np.copy(position)
        params = self._param_values.copy()
        params[self._rate_mask] = 10**Y
        sim = self._model_solver.run(param_values=[params]).all
        logl = 0.0
        for observable in self._like_data.keys():
            sim_vals = sim[observable][self._data_mask[observable]]
            logl += np.sum((self._data[observable]-sim_vals)**2)
        if np.isnan(logl):
            return np.inf
        return logl,

    def __call__(self, pso_kwargs=None,
                 cost_type='norm_logpdf'):
        """Call the SwarmIt instance to construct to instance of the NestedSampling object.

        Args:
                pso_kwargs (dict): Dictionary of any additional optional keyword
                    arguments to pass to the PSO object constructor.
                    Defaults to dict().
                cost_type (str): Define the type of cost estimator
                    to use. Options are 'norm_logpdf'=>Compute the cost using
                    the normal distribution estimator, 'mse'=>Compute the
                    cost using the negative mean squared error estimator,
                    'sse'=>Compute the cost using the negative sum of
                     squared errors estimator. Defaults to 'norm_logpdf'.

        Returns:
            type: Description of returned object.

        """
        if pso_kwargs is None:
            pso_kwargs = dict()
        # self.ns_version = ns_version
        self._pso_kwargs = pso_kwargs
        #population_size = pso_population_size
        if cost_type == 'mse':
            cost = self.mse_cost
        elif cost_type == 'sse':
            cost = self.sse_cost
        else:
            cost = self.norm_logpdf_cost

        # Construct the PSO
        if 'save_sampled' not in pso_kwargs.keys():
            pso_kwargs['save_sampled'] = False
        if 'verbose' not in pso_kwargs.keys():
            pso_kwargs['verbose'] = False
        pso = PSO(**pso_kwargs)
        pso.set_start_position(self._starting_position)
        pso.set_cost_function(cost)
        pso.set_bounds(lower=self._lower, upper=self._upper)
        pso.set_speed(-.25, .25)
        return pso


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('model_file', metavar='model_file', type=str, help='The model.py file that you want run PSO on.')
    parser.add_argument('output_path', metavar='output_path', type=str, help='The file location where you want to output the PSO run script.')
    args = parser.parse_args()
    # get the input script from command line imputs
    model_file = os.path.abspath(args.model_file)
    # get the location to dump the output
    output_path = os.path.abspath(args.output_path)

    print("Using model from file: {}".format(model_file))
    base=os.path.basename(model_file)
    model_module_name = os.path.splitext(base)[0]
    print("With name: {}".format(model_module_name))
    #default_prior_shape = 'norm'
    #print("The default prior shape is: {}".format(default_prior_shape))
    #print(model_file)

    #print(model_module_name)
    model_module = importlib.import_module(model_module_name)
    model = getattr(model_module, 'model')

    try:
        swarm_param = getattr(model_module, 'swarm_param')
    except:
        swarm_param = None
    #print(model)
    priors = dict()
    no_sample = list()
    Keq_sample = list()
    #Read the file and parse any #SWARM_IT directives
    print("Parsing the model for any #SWARM_IT directives...")
    with open(model_file, 'r') as file_obj:
        for line in file_obj:
            words = line.split()
            if len(words) > 1:
                if words[0] == '#SWARM_IT':
                    parse_directive(line, priors, no_sample, Keq_sample)

    #now we need to extract a list of kinetic parameters
    parameters = list()
    print("Inspecting the model and pulling out kinetic parameters...")
    for rule in model.rules:
        print(rule.rate_forward, rule.rate_reverse)
        #print(rule_keys)
        if rule.rate_forward:
            param = rule.rate_forward
            #print(param)
            parameters.append([param,'f'])
        if rule.rate_reverse:
            param = rule.rate_reverse
            #print(param)
            parameters.append([param, 'r'])
    #print(no_sample)
    parameters = prune_no_samples(parameters, no_sample)
    parameters = update_with_Keq_samples(parameters, Keq_sample)
    #print(parameters)
    print("Found the following kinetic parameters:")
    print("{}".format(parameters))
    #default the priors to norm - i.e. normal distributions
    #for parameter in parameters:
    #    name = parameter[0].name
    #    if name not in priors.keys():
    #        priors[name] = default_prior_shape

    # Obtain mask of sampled parameters to run simulation in the likelihood function
    parameters_idxs = [model.parameters.index(parameter[0]) for parameter in parameters]
    if swarm_param is None:
        calibrate_mask = [i in parameters_idxs for i in range(len(model.parameters))]
    else:
        calibrate_mask = swarm_param.mask(model.parameters)
    param_values = [p.value for p in model.parameters]

    out_file = open("simplepso_"+base, 'w')

    print("Writing to simplePSO run script: simplepso_{}".format(base))
    out_file.write("\'\'\'\nGenerated by swarm_it\n")
    out_file.write("simplePSO run script for {} \n".format(base))
    out_file.write("\'\'\'")
    out_file.write("\n")
    out_file.write("from pysb.simulator import ScipyOdeSimulator\n")
    out_file.write("import numpy as np\n")
    out_file.write("from scipy.stats import norm\n")
    out_file.write("from simplepso.pso import PSO\n")

    #out_file.write("import inspect\n")
    #out_file.write("import os.path\n")
    if swarm_param is not None:
        out_file.write("from "+model_module_name+" import model,swarm_param\n")
    else:
        out_file.write("from "+model_module_name+" import model\n")
    out_file.write("\n")

    out_file.write("# Initialize PySB solver object for running simulations.\n")
    out_file.write("# USER-Adjust: simulation timespan should match experimental data.\n")
    out_file.write("tspan = np.linspace(0,10, num=100)\n")
    out_file.write("solver = ScipyOdeSimulator(model, tspan=tspan)\n")
    out_file.write("parameters_idxs = " + str(parameters_idxs)+"\n")
    out_file.write("calibrate_mask = " + str(calibrate_mask)+"\n" )
    out_file.write("param_values = np.array([p.value for p in model.parameters])\n" )
    out_file.write("\n")
    out_file.write("# USER-Set: must add commands to import/load any experimental\n")
    out_file.write("# data for use in the likelihood function!\n")
    out_file.write("experiments_avg = np.load()\n")
    out_file.write("experiments_sd = np.load()\n")
    out_file.write("like_data = norm(loc=experiments_avg, scale=experiments_sd)\n")
    out_file.write("# USER-Set: must appropriately update cost function!\n")
    out_file.write("def cost(position):\n")
    out_file.write("    Y=np.copy(position)\n")
    out_file.write("    param_values[calibrate_mask] = 10 ** Y\n")
    out_file.write("    sim = solver.run(param_values=param_values).all\n")
    out_file.write("    logp_data = np.sum(like_data.logpdf(sim['observable']))\n")
    out_file.write("    if np.isnan(logp_data):\n")
    out_file.write("        logp_data = np.inf\n")
    out_file.write("    return -logp_data,\n")
    out_file.write("\n")
    #write the sampled params lines
    out_file.write("# Setup the particle swarm optimization run\n \n")
    #out_file.write("# Setup the parameters that being sampled. \n")
    #out_file.write("sampled_parameters = list()\n")
    #for parameter in parameters:
    #    name = parameter[0].name
    #    prior_shape = 'uniform'
    #    print("Will sample parameter {} with starting position {}".format(name, value))
        #if prior_shape == 'uniform':
        #    line = write_uniform_param(name, value)
        #    ps_name = line.split()[0]
        #    out_file.write(line)
        #    out_file.write("sampled_parameters.append({})\n".format(ps_name))

    #out_file.write("n_params = len(sampled_parameters)\n")
    out_file.write("# Set the number of particles in the swarm. \n")
    out_file.write("num_particles = 25\n")
    out_file.write("# Set the number of iterations for PSO run. \n")
    out_file.write("num_iterations = 50\n")
    out_file.write("# Construct the optimizer \n")
    out_file.write("pso = PSO(save_sampled=False,\n")
    out_file.write("          verbose=True,\n")
    out_file.write("          num_procs=1) \n")
    out_file.write("pso.set_cost_function(cost)\n")
    if swarm_param is None:
        out_file.write("starting_position = param_values[calibrate_mask]\n")
        out_file.write("pso.set_start_position(starting_position)\n")
        out_file.write("# allows particles to move +/- 2 orders of magnitude\n")
        out_file.write("pso.set_bounds(2)\n")
    else:
        out_file.write("starting_position = swarm_param.centers()\n")
        out_file.write("pso.set_start_position(starting_position)\n")
        out_file.write("pso.set_bounds(lower=swarm_param.lower(), upper=swarm_param.upper())\n")

    out_file.write("# sets maximum speed that a particle can travel\n")
    out_file.write("pso.set_speed(-.25, .25)\n")
    out_file.write("# run it\n")
    out_file.write("pso.run(num_particles,\n")
    out_file.write("        num_iterations,\n")
    out_file.write("        stop_threshold=1e-5)\n")
    out_file.write("print(\"Best parameters: \",pso.best)\n")
    #out_file.write("print(\"Best parameters cost: \", cost_val)\n")


    out_file.close()
    print("swarm_it is complete!")
    print("END OF LINE.")
