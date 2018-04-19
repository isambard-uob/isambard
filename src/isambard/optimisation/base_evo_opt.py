"""Base class for bio-inspired optimizers."""

from concurrent import futures
import datetime
import enum
import os
import sys

import matplotlib.pylab as plt
import numpy
from deap import base, creator, tools

from ..modelling import pack_side_chains_scwrl


def default_build(spec_seq_params):
    specification, sequences, params = spec_seq_params
    model = specification(*params)
    model = pack_side_chains_scwrl(model, sequences)
    return model


class BaseOptimizer:
    """Abstract base class for the evolutionary optimizers.

    Notes
    -----
    Not intended to be used directly, see the evolutionary
    optimizers for full documentation.
    """

    def __init__(self, specification, sequences, parameters,
                 build_fn=None, eval_fn=None, **kwargs):
        self.specification = specification
        self.sequences = sequences
        self.parameters = parameters
        self._make_parameters()
        if sys.platform == 'win32':
            self.mp_disabled = True
            print('Multiprocessing for this module is currently unavailable'
                  'on Windows, only a single process will be used.')
        else:
            if 'mp_disabled' in kwargs:
                self.mp_disabled = kwargs['mp_disabled']
            else:
                self.mp_disabled = False
        self.build_fn = build_fn
        self.eval_fn = eval_fn
        self.population = None
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        self.toolbox = base.Toolbox()
        self.parameter_log = []
        self._cores = 1
        self._evals = 0
        self._model_count = 0
        self._store_params = True

    def parse_individual(self, individual):
        """Converts a deap individual into a full list of parameters.

        Parameters
        ----------
        individual: deap individual from optimization
            Details vary according to type of optimization, but
            parameters within deap individual are always between -1
            and 1. This function converts them into the values used to
            actually build the model

        Returns
        -------
        fullpars: list
            Full parameter list for model building.
        """
        scaled_ind = []
        for i in range(len(self.value_means)):
            scaled_ind.append(self.value_means[i] + (
                individual[i] * self.value_ranges[i]))
        fullpars = list(self.arrangement)
        for k in range(len(self.variable_parameters)):
            for j in range(len(fullpars)):
                if fullpars[j] == self.variable_parameters[k]:
                    fullpars[j] = scaled_ind[k]
        return fullpars

    def run_opt(self, pop_size, generations, cores=1, plot=False, log=False,
                log_path=None, run_id=None, store_params=True, **kwargs):
        """Runs the optimizer.

        Parameters
        ----------
        pop_size: int
            Size of the population each generation.
        generation: int
            Number of generations in optimisation.
        cores: int, optional
            Number of CPU cores used to run the optimisation.
            If the 'mp_disabled' keyword is passed to the
            optimizer, this will be ignored and one core will
            be used.
        plot: bool, optional
            If true, matplotlib will be used to plot information
            about the minimisation.
        log: bool, optional
            If true, a log file describing the optimisation will
            be created. By default it will be written to the
            current directory and named according to the time the
            minimisation finished. This can be manually specified
            by passing the 'output_path' and 'run_id' keyword
            arguments.
        log_path : str
            Path to write output file.
        run_id : str
            An identifier used as the name of your log file.
        store_params: bool, optional
            If true, the parameters for each model created during
            the optimisation will be stored. This can be used to
            create funnel data later on.
        """
        self._cores = cores
        self._store_params = store_params
        self.parameter_log = []
        self._model_count = 0
        self.halloffame = tools.HallOfFame(1)
        self.stats = tools.Statistics(lambda thing: thing.fitness.values)
        self.stats.register("avg", numpy.mean)
        self.stats.register("std", numpy.std)
        self.stats.register("min", numpy.min)
        self.stats.register("max", numpy.max)
        self.logbook = tools.Logbook()
        self.logbook.header = ["gen", "evals"] + self.stats.fields
        start_time = datetime.datetime.now()
        self._initialize_pop(pop_size)
        for g in range(generations):
            self._update_pop(pop_size)
            self.halloffame.update(self.population)
            self.logbook.record(gen=g, evals=self._evals,
                                **self.stats.compile(self.population))
            print(self.logbook.stream)
        end_time = datetime.datetime.now()
        time_taken = end_time - start_time
        print("Evaluated {} models in total in {}".format(
            self._model_count, time_taken))
        print("Best fitness is {0}".format(self.halloffame[0].fitness))
        print("Best parameters are {0}".format(self.parse_individual(
            self.halloffame[0])))
        for i, entry in enumerate(self.halloffame[0]):
            if entry > 0.95:
                print(
                    "Warning! Parameter {0} is at or near maximum allowed "
                    "value\n".format(i + 1))
            elif entry < -0.95:
                print(
                    "Warning! Parameter {0} is at or near minimum allowed "
                    "value\n".format(i + 1))
        if log:
            self.log_results(output_path=output_path, run_id=run_id)
        if plot:
            print('----Minimisation plot:')
            plt.figure(figsize=(5, 5))
            plt.plot(range(len(self.logbook.select('min'))),
                     self.logbook.select('min'))
            plt.xlabel('Iteration', fontsize=20)
            plt.ylabel('Score', fontsize=20)
        return

    def _make_parameters(self):
        """Converts a list of Parameters into DEAP format."""
        self.value_means = []
        self.value_ranges = []
        self.arrangement = []
        self.variable_parameters = []
        current_var = 0
        for parameter in self.parameters:
            if parameter.type == ParameterType.DYNAMIC:
                self.value_means.append(parameter.value[0])
                if parameter.value[1] < 0:
                    raise AttributeError(
                        '"{}" parameter has an invalid range. Range values '
                        'must be greater than zero'.format(parameter.label))
                self.value_ranges.append(parameter.value[1])
                var_label = 'var{}'.format(current_var)
                self.arrangement.append(var_label)
                self.variable_parameters.append(var_label)
                current_var += 1
            elif parameter.type == ParameterType.STATIC:
                self.arrangement.append(parameter.value)
            else:
                raise AttributeError(
                    'Unknown parameter type ({}). Parameters can be STATIC or'
                    ' DYNAMIC.'.format(parameter.type))
        return

    def assign_fitnesses(self, targets):
        """Assigns fitnesses to parameters.

        Notes
        -----
        Uses `self.eval_fn` to evaluate each member of target.
        Parameters
        ---------

        targets
            Parameter values for each member of the population.
        """
        self._evals = len(targets)
        px_parameters = zip([self.specification] * len(targets),
                            [self.sequences] * len(targets),
                            [self.parse_individual(x) for x in targets])
        if (self._cores == 1) or (self.mp_disabled):
            models = map(self.build_fn, px_parameters)
            fitnesses = map(self.eval_fn, models)
        else:
            with futures.ProcessPoolExecutor(
                    max_workers=self._cores) as executor:
                models = executor.map(self.build_fn, px_parameters)
                fitnesses = executor.map(self.eval_fn, models)
        tars_fits = list(zip(targets, fitnesses))
        if self._store_params:
            self.parameter_log.append(
                [(self.parse_individual(x[0]), x[1]) for x in tars_fits])
        for ind, fit in tars_fits:
            ind.fitness.values = (fit,)
        return

    def _generate(self):
        """Generates a particle using the creator function."""
        raise NotImplementedError("Will depend on optimizer type")

    def _initialize_pop(self):
        """Assigns indices to individuals in population."""
        raise NotImplementedError("Will depend on optimizer type")

    def _update_pop(self):
        """Updates population according to crossover and fitness criteria."""
        raise NotImplementedError("Will depend on optimizer type")

    def log_results(self, output_path=None, run_id=None):
        """Saves files for the minimization.

        Notes
        -----
        Currently saves a logfile with best individual and a pdb of
        the best model.
        """
        best_ind = self.halloffame[0]
        model_params = self.parse_individual(
            best_ind)  # need to change name of 'params'
        if output_path is None:
            output_path = os.getcwd()
        if run_id is None:
            run_id = '{:%Y%m%d-%H%M%S}'.format(
                datetime.datetime.now())
        with open('{0}/{1}_opt_log.txt'.format(
                output_path, run_id), 'w') as log_file:
            log_file.write('\nEvaluated {0} models in total\n'.format(
                self._model_count))
            log_file.write('Run ID is {0}\n'.format(run_id))
            log_file.write('Best fitness is {0}\n'.format(
                self.halloffame[0].fitness))
            log_file.write(
                'Parameters of best model are {0}\n'.format(model_params))
            log_file.write(
                'Best individual is {0}\n'.format(self.halloffame[0]))
            for i, entry in enumerate(self.halloffame[0]):
                if entry > 0.95:
                    log_file.write(
                        "Warning! Parameter {0} is at or near maximum allowed "
                        "value\n".format(i + 1))
                elif entry < -0.95:
                    log_file.write(
                        "Warning! Parameter {0} is at or near minimum allowed "
                        "value\n".format(i + 1))
            log_file.write('Minimization history: \n{0}'.format(self.logbook))
        with open('{0}/{1}_opt_best_model.pdb'.format(
                output_path, run_id), 'w') as output_file:
            output_file.write(self.best_model.pdb)
        return

    @property
    def best_model(self):
        """Rebuilds the top scoring model from an optimisation.

        Returns
        -------
        model: AMPAL
            Returns an AMPAL model of the top scoring parameters.

        Raises
        ------
        AttributeError
            Raises a name error if the optimiser has not been run.
        """
        if not hasattr(self, 'halloffame'):
            raise AttributeError(
                'No best model found, have you ran the optimiser?')
        model = self.build_fn(
            (self.specification,
             self.sequences,
             self.parse_individual(self.halloffame[0])
             ))
        return model

    def make_energy_funnel_data(self, cores=1):
        """Compares models created during the minimisation to the best model.

        Parameters
        ----------
        cores : int
            Number of CPU cores to rebuild models and measure RMSD.

        Returns
        -------
        energy_rmsd_gen: [(float, float, int)]
            A list of triples containing the BUFF score, RMSD to the
            top model and generation of a model generated during the
            minimisation.
        """
        if not self.parameter_log:
            raise AttributeError(
                'No parameter log data to make funnel, have you ran the '
                'optimiser?')
        model_cls = self.specification
        gen_tagged = []
        for gen, models in enumerate(self.parameter_log):
            for model in models:
                gen_tagged.append((model[0], model[1], gen))
        sorted_pps = sorted(gen_tagged, key=lambda x: x[1])
        top_result = sorted_pps[0]
        top_result_model = model_cls(*top_result[0])
        if (cores == 1) or (sys.platform == 'win32'):
            energy_rmsd_gen = map(
                self.funnel_rebuild,
                [(x, top_result_model, self.specification)
                 for x in gen_tagged])
        else:
            with futures.ProcessPoolExecutor(max_workers=cores) as executor:
                energy_rmsd_gen = executor.map(
                    self.funnel_rebuild,
                    [(x, top_result_model, self.specification)
                     for x in gen_tagged])
        return list(energy_rmsd_gen)

    @staticmethod
    def funnel_rebuild(psg_trm_spec):
        """Rebuilds a model and compares it to a reference model.

        Parameters
        ----------
        psg_trm: (([float], float, int), AMPAL, specification)
            A tuple containing the parameters, score and generation for a
            model as well as a model of the best scoring parameters.

        Returns
        -------
        energy_rmsd_gen: (float, float, int)
            A triple containing the BUFF score, RMSD to the top model
            and generation of a model generated during the minimisation.
        """
        param_score_gen, top_result_model, specification = psg_trm_spec
        params, score, gen = param_score_gen
        model = specification(*params)
        rmsd = top_result_model.rmsd(model)
        return rmsd, score, gen


class ParameterType(enum.IntEnum):
    """Type of parameter used in evolutionary optimizers."""
    STATIC = 1
    DYNAMIC = 2


class Parameter:
    """Defines a parameter used in the evolutionary optimizers."""

    def __init__(self, label, parameter_type, value):
        self.label = label
        self.type = parameter_type
        self.value = value

    def __repr__(self):
        return "<Parameter: {}, {}>".format(self.label, str(self.type))

    @classmethod
    def static(cls, label, value):
        """Creates a static parameter.

        Parameters
        ----------
        label : str
            A human-readable label for the parameter.
        value
            The static value to be used.
        """
        return cls(label, ParameterType.STATIC, value)

    @classmethod
    def dynamic(cls, label, val_mean, val_range):
        """Creates a static parameter.

        Parameters
        ----------
        label : str
            A human-readable label for the parameter.
        val_mean : float
            The mean value of the parameter.
        val_range : float
            The minimum and maximum variance from the mean allowed for
            parameter.
        """
        return cls(label, ParameterType.DYNAMIC, (val_mean, val_range))

    @property
    def default_value(self):
        if self.type == ParameterType.STATIC:
            return self.value
        elif self.type == ParameterType.DYNAMIC:
            return self.value[0]
        else:
            raise AttributeError('"{}" is an unknown parameter type.')


__author__ = 'Andrew R. Thomson, Christopher W. Wood, Gail J. Bartlett'
