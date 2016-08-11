from deap import algorithms, base, creator, tools
import random
import numpy
import matplotlib.pylab as plt
import datetime
from external_programs.profit import run_profit
from concurrent import futures

# TODO: way of getting best model in different subclasses


class ParentDE:
    """Parent class for DE optimizers. Not meant to be instantiated.

    Parameters
    ----------
    specification: Isambard specification
        Tells the optimizer what kind of structure it is building
    output_path: str
        Specifies the location to save files if logging is enabled
    mode: 'buff', 'rmsd' or 'comparator'
        Specifies whether the score should be the average per chain or total cumulative score.
    """
    def __init__(self):
        self._de_params = locals()
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMin)
        self.toolbox = base.Toolbox()

    def parameters(self, sequence, value_means, value_ranges, arrangement):
        """Relates the individual to be evolved to the full parameter string for building the specification object

        Parameters
        ----------
        sequence: str
            Full amino acid sequence for specification object to be optimized. Must be equal to the number of residues in the
            model.
        value_means: list
            List containing mean values for parameters to be optimized.
        value_ranges: list
            List containing ranges for parameters to be optimized. Values must be positive.
        arrangement: list
            Full list of fixed and variable parameters for model building. Fixed values are the appropriate value.
            Values to be varied should be listed as 'var0', 'var1' etc, and must be in ascending numerical order.
            Variables can be repeated if required.
        """
        self._de_params.update(locals())
        if any(x <= 0 for x in self._de_params['value_ranges']):
            raise ValueError("range values must be greater than zero")
        self._de_params['variable_parameters'] = []
        for i in range(len(self._de_params['value_means'])):
            self._de_params['variable_parameters'].append("".join(['var', str(i)]))
        if len(set(arrangement).intersection(self._de_params['variable_parameters'])) !=\
                len(self._de_params['value_means']):
            raise ValueError("argument mismatch!")
        if len(self._de_params['value_ranges']) != len(self._de_params['value_means']):
            raise ValueError("argument mismatch!")
        self.toolbox.register("individual", self.generate)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)

    def generate(self):
        """Generates a particle using the creator function. Position and speed are uniformly randomly seeded within
        allowed bounds. The particle also has speed limit settings taken from global values.

        Returns
        -------
        particle object
        """
        ind = creator.Individual([random.uniform(-1, 1) for _ in range(len(self._de_params['value_means']))])
        ind.index = None
        ind.neighbours = None
        return ind

    def parse_individual(self, individual):
        """Converts an individual from the PSO into a full list of parameters for building the specification object.
        Parameters
        ----------
        individual: position elements from swarm particle

        Returns
        -------
        fullpars: list
            Full parameter list for model building.
        """
        scaled_ind = []
        for i in range(len(self._de_params['value_means'])):
            scaled_ind.append(self._de_params['value_means'][i] + (individual[i] * self._de_params['value_ranges'][i]))
        fullpars = list(self._de_params['arrangement'])
        for k in range(len(self._de_params['variable_parameters'])):
            for j in range(len(fullpars)):
                if fullpars[j] == self._de_params['variable_parameters'][k]:
                    fullpars[j] = scaled_ind[k]
        return fullpars

    def crossover(self, ind):
        """Used by the evolution process to generate a new individual

        Parameters
        ----------
        pop_index: int
            The position of the individual within the population. A candidate solution is generated using two other
            randomly selected individuals.
        Returns
        -------
        y: deap individual
            An individual representing a candidate solution, to be assigned a fitness.
        """
        if self._de_params['neighbours']:
            a, b, c = random.sample([self.pop[i] for i in ind.neighbours], 3)
        else:
            a, b, c = random.sample(self.pop, 3)
        y = self.toolbox.clone(a)
        y.index = ind.index
        y.neighbours = ind.neighbours
        del y.fitness.values
        # y should now be a copy of ind with the vector elements from a
        index = random.randrange(len(self._de_params['value_means']))
        for i, value in enumerate(y):
            if i == index or random.random() < self._de_params['cxpb']:
                entry = a[i] + random.lognormvariate(-1.2, 0.5)*self._de_params['diff_weight']*(b[i]-c[i])
                tries = 0
                while abs(entry) > 1.0:
                    tries += 1
                    entry = a[i] + random.lognormvariate(-1.2, 0.5)*self._de_params['diff_weight']*(b[i]-c[i])
                    if tries > 10000:
                        entry = a[i]
                y[i] = entry
        return y

    def run_de(self, gensize, numgen, procs, diff_weight=1, cxpb=0.75, neighbours=False, plot=False, log=False,
               run_id='DE_model'):
        """

        Parameters
        ----------
        gensize: int
            The number of individuals to evolve. Values of 20-50 are typical.
        numgen: int
            The number of generations over which the evolution will be run. A value of 50 is typical.
        procs: int
            The number of processors to use
        diff_weight: float
            A parameter determining how much of the crossover product is added to the candidate solution. Low values
            favour minimization, high values favour exploration. Values in the range 0.5 - 2.0 are typical.
        cxpb: float
            The crossover probability. Determines how likely a given element of an individual is to undergo mutation.
        plot: bool
            Whether or not to plot the minimization graph
        log: bool
            Whether or not to save files for the best model and minimization history.
        run_id: str
            Identifier for saving files if logging is enabled
        """
        self._de_params.update(locals())
        if self._de_params['log'] is True:
            if self._de_params['output_path'] is None:
                raise ValueError("Output path must be specified to enable logging")
        if self._de_params['cxpb'] > 1 or self._de_params['cxpb'] < 0:
            raise ValueError("Crossover probability must be between 0 and 1")

        start_time = datetime.datetime.now()
        print('Starting minimisation ({:%Y-%m-%d %H:%M:%S})'.format(start_time))
        self.pop = self.toolbox.population(n=self._de_params['gensize'])
        if self._de_params['neighbours']:
            for i in range(len(self.pop)):
                self.pop[i].index = i
                self.pop[i].neighbours = list(set([(i-x) % len(self.pop)
                                                   for x in range(1, self._de_params['neighbours']+1)] +
                                                  [(i+x) % len(self.pop)
                                                   for x in range(1, self._de_params['neighbours']+1)]))
        self.halloffame = tools.HallOfFame(1)
        random.seed()
        self.stats = tools.Statistics(lambda thing: thing.fitness.values)
        self.stats.register("avg", numpy.mean)
        self.stats.register("std", numpy.std)
        self.stats.register("min", numpy.min)
        self.stats.register("max", numpy.max)
        self.logbook = tools.Logbook()
        self.logbook.header = ["gen", "evals"] + self.stats.fields
        self.assign_fitness(self.pop)
        self.halloffame.update(self.pop)
        modelcount = len(self.pop)

        for g in range(self._de_params['numgen']):
            candidates = []
            for ind in self.pop:
                candidates.append(self.crossover(ind))
            modelcount += len(candidates)
            self.assign_fitness(candidates)
            for i in range(len(self.pop)):
                if candidates[i].fitness.values[0] < self.pop[i].fitness.values[0]:
                    self.pop[i] = candidates[i]
            self.halloffame.update(self.pop)
            self.logbook.record(gen=g, evals=len(candidates), **self.stats.compile(self.pop))
            print(self.logbook.stream)
        end_time = datetime.datetime.now()
        time_taken = end_time - start_time
        self._de_params['time_taken'] = time_taken
        self._de_params['model_count'] = modelcount
        print("End of minimisation ({:%Y-%m-%d %H:%M:%S})".format(end_time))
        print("Run ID is {0}".format(self._de_params['run_id']))
        print('Minimization time = {0}'.format(time_taken))
        print("Evaluated {0} models in total".format(self._de_params['model_count']))
        print("Best score is {0}".format(self.halloffame[0].fitness.values[0]))
        print("Best parameters are {0}".format(self.parse_individual(self.halloffame[0])))
        for i, entry in enumerate(self.halloffame[0]):
            if entry > 0.95:
                print("Warning! Parameter {0} is at or near maximum allowed value\n".format(i+1))
            elif entry < -0.95:
                print("Warning! Parameter {0} is at or near minimum allowed value".format(i+1))
        if self._de_params['log']:
            self.log_results()
        if self._de_params['plot']:
            print('----Minimisation plot:')
            plt.figure(figsize=(5, 5))
            plt.plot(range(len(self.logbook.select('min'))), self.logbook.select('min'))
            plt.xlabel('Iteration', fontsize=20)
            plt.ylabel('Score', fontsize=20)

    def assign_fitness(self, targets):
        pass

    def log_results(self):
        """Saves files for the minimization. Currently saves a logfile with best individual and a pdb of the best model
        """
        best_ind = self.halloffame[0]
        params = self.parse_individual(best_ind)
        with open('{0}{1}_log.txt'.format(self._de_params['output_path'], self._de_params['run_id']), 'a+') as log_file:
            log_file.write('\nEvaluated {0} models in total\n'.format(self._de_params['model_count']))
            log_file.write('Run ID is {0}\n'.format(self._de_params['run_id']))
            log_file.write('Best score is {0}\n'.format(self.halloffame[0].fitness.values[0]))
            log_file.write('Parameters of best model are {0}'.format(params))
            log_file.write('Best individual is {0}\n'.format(self.halloffame[0]))
            for i, entry in enumerate(self.halloffame[0]):
                if entry > 0.95:
                    log_file.write("Warning! Parameter {0} is at or near maximum allowed value\n".format(i+1))
                elif entry < -0.95:
                    log_file.write("Warning! Parameter {0} is at or near minimum allowed value\n".format(i+1))
            log_file.write('Minimization history: \n{0}'.format(self.logbook))
        with open('{0}{1}_bestmodel.pdb'.format(self._de_params['output_path'], self._de_params['run_id']),
                  'w') as output_file:
            model = self._de_params['specification'](*params)
            model.build()
            model.pack_new_sequences(self._de_params['sequence'])
            output_file.write(model.pdb)


class OptDE(ParentDE):
    def __init__(self, specification, output_path=None, mode='buff'):
        super(OptDE, self).__init__()
        self._de_params = locals()

    def assign_fitness(self, targets):
        px_parameters = zip([self._de_params['specification']] * len(targets),
                            [self._de_params['sequence']] * len(targets),
                            [self.parse_individual(x) for x in targets])
        with futures.ProcessPoolExecutor(max_workers=self._de_params['procs']) as executor:
            fitnesses = executor.map(buff_eval, px_parameters)
        for ind, fit in zip(targets, fitnesses):
            ind.fitness.values = (fit,)


class OptDE_rmsd(ParentDE):
    def __init__(self, specification, ref_pdb, output_path=None, mode='buff', cmd_file_path=None):
        super(OptDE_rmsd, self).__init__()
        self._de_params = locals()

    def assign_fitness(self, targets):
        px_parameters = zip([self._de_params['specification']] * len(targets),
                            [self._de_params['sequence']] * len(targets),
                            [self.parse_individual(x) for x in targets],
                            [self._de_params['ref_pdb']] * len(targets),
                            [self._de_params['cmd_file_path']] * len(targets))
        with futures.ProcessPoolExecutor(max_workers=self._de_params['procs']) as executor:
            fitnesses = executor.map(rmsd_eval, px_parameters)
        for ind, fit in zip(targets, fitnesses):
            ind.fitness.values = (fit,)


class OptDE_comparator(ParentDE):
    def __init__(self, top1, top2, params1, params2, seq1, seq2, output_path=None, mode='buff'):
        super(OptDE_comparator, self).__init__()
        self._de_params = locals()
        obj1 = top1(*params1)
        obj1.pack_new_sequences(seq1)
        obj2 = top2(*params2)
        obj2.pack_new_sequences(seq2)
        self._de_params['ref1'] = obj1.buff_interaction_energy.total_energy
        self._de_params['ref2'] = obj2.buff_interaction_energy.total_energy

    def assign_fitness(self, targets):
        px_parameters = zip([self._de_params['top1']] * len(targets),
                            [self._de_params['top2']] * len(targets),
                            [self._de_params['params1']] * len(targets),
                            [self._de_params['params2']] * len(targets),
                            [self._de_params['seq1']] * len(targets),
                            [self._de_params['seq2']] * len(targets),
                            [self.parse_individual(x) for x in targets])
        with futures.ProcessPoolExecutor(max_workers=self._de_params['procs']) as executor:
            fitnesses = executor.map(comparator_eval, px_parameters)
        for ind, fit in zip(targets, fitnesses):
            ind.fitness.values = (fit - (self._de_params['ref1'] + self._de_params['ref2']),)

    def parameters(self, value_means, value_ranges, arrangement):
        """Relates the individual to be evolved to the full parameter string for building the specification object

        Parameters
        ----------
        sequence: str
            Full amino acid sequence for specification object to be optimized. Must be equal to the number of residues in the
            model.
        value_means: list
            List containing mean values for parameters to be optimized.
        value_ranges: list
            List containing ranges for parameters to be optimized. Values must be positive.
        arrangement: list
            Full list of fixed and variable parameters for model building. Fixed values are the appropriate value.
            Values to be varied should be listed as 'var0', 'var1' etc, and must be in ascending numerical order.
            Variables can be repeated if required.
        """
        self._de_params.update(locals())
        if any(x <= 0 for x in self._de_params['value_ranges']):
            raise ValueError("range values must be greater than zero")
        self._de_params['variable_parameters'] = []
        for i in range(len(self._de_params['value_means'])):
            self._de_params['variable_parameters'].append("".join(['var', str(i)]))
        if len(set(arrangement).intersection(self._de_params['variable_parameters'])) !=\
                len(self._de_params['value_means']):
            raise ValueError("argument mismatch!")
        if len(self._de_params['value_ranges']) != len(self._de_params['value_means']):
            raise ValueError("argument mismatch!")
        self.toolbox.register("individual", self.generate)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)


def buff_eval(params):
    """Special version of the eval function that makes no reference to .self and so can be farmed out to other
    namespaces.

    Notes
    -----
    This is set up for parallel execution so makes no reference to local namespace

    Parameters
    ----------
    params: list
        Tuple containing the specification to be built, the sequence, and the parameters for model building.

    Returns
    -------
    model.bude_score: float
        BUFF score for model to be assigned to particle fitness value.
    """
    specification, sequence, parsed_ind = params
    model = specification(*parsed_ind)
    model.build()
    model.pack_new_sequences(sequence)
    return model.buff_interaction_energy.total_energy


def rmsd_eval(rmsd_params):
    """
    Builds a model based on an individual from the optimizer and runs profit against a reference model.

    Parameters
    ----------
    rmsd_params

    Returns
    -------
    rmsd: float
        rmsd against reference model as calculated by profit.
    """
    specification, sequence, parsed_ind, reference_pdb, cmd_file_path = rmsd_params
    model = specification(*parsed_ind)
    model.pack_new_sequences(sequence)
    try:
        ca, bb, aa = run_profit(reference_pdb, model.pdb, path1=False, path2=False, path_to_cmd_file=cmd_file_path)
    except ValueError:
        return 100.0
    return bb


def comparator_eval(comparator_params):
    """Gets BUFF score for interaction between two AMPAL objects
    """
    top1, top2, params1, params2, seq1, seq2, movements = comparator_params
    xrot, yrot, zrot, xtrans, ytrans, ztrans = movements
    obj1 = top1(*params1)
    obj2 = top2(*params2)
    obj2.rotate(xrot, [1, 0, 0])
    obj2.rotate(yrot, [0, 1, 0])
    obj2.rotate(zrot, [0, 0, 1])
    obj2.translate([xtrans, ytrans, ztrans])
    model = obj1 + obj2
    model.relabel_all()
    model.pack_new_sequences(seq1 + seq2)
    return model.buff_interaction_energy.total_energy


__author__ = 'Andrew R. Thomson'
__status__ = 'Development'
