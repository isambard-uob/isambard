from deap import algorithms, base, creator, tools
import random
import numpy
import matplotlib.pylab as plt
import datetime
from concurrent import futures

# TODO: Other crossover types? Might be interesting to test e.g. cxSimulatedBinaryBounded
# TODO: Overhaul- eliminate non parallel version (copy subclass to main class)


class OptGA:
    """Class for running a genetic algorithm using the DEAP framework

    Notes
    -----
    This system makes use of an internal toolbox of functions that all originate from DEAP.
    This includes selection, mutation, mating, etc etc, as well as decorators.
    A population of individuals is generated, and each individual is evaluated. Individuals are a list of parameters,
    with an addition ind.fitness.values property that is written when they are evaluated.
    This has been set up such that each individual can represent an evolving sequence of parameters that are parsed
    to give a full set of parameters for a specification object, allowing for invariant parameters in the specification object, or
    e.g. one parameter from the evolving individual to be used to specify two (identical but variable) parameters in the
    specification object.
    On instantiation it needs a specification to evaluate and an output path.

    Parameters
    ----------
    specification: class
        the specific specification to optimise
    output path: str
        location to save output file
    """

    def __init__(self, specification, output_path, bude_mode='average'):
        self.output_path = output_path
        self.specification = specification
        self.bude_mode = bude_mode
        self.sequence = None
        self.value_ranges = None
        self.value_means = None
        self.arrangement = None
        self.keys = None
        self.log = None
        self.modelcount = None
        self.pop = None
        self.procs = None
        self.plot = None
        self.threshold = None
        self.gensize = None
        self.numgen = None
        self.cxpb = None
        self.mutpb = None
        self.log = None
        self.genbest = None
        self.run_id = None
        self.start_time = None
        self.end_time = None
        self.halloffame = None
        self.stats = None
        self.logbook = None
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
        self.value_ranges = value_ranges
        if any(x <= 0 for x in self.value_ranges):
            raise ValueError("range values must be greater than zero")
        self.value_means = value_means
        self.sequence = sequence
        self.arrangement = arrangement
        self.keys = []
        for i in range(len(self.value_means)):
            self.keys.append("".join(['var', str(i)]))
        if len(set(arrangement).intersection(self.keys)) != len(self.value_means):
            raise ValueError("argument mismatch!")
        if len(self.value_ranges) != len(self.value_means):
            raise ValueError("argument mismatch!")

        self.toolbox.register("params", self.makepars)
        self.toolbox.register("individual", tools.initIterate, creator.Individual, self.toolbox.params)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)

        self.toolbox.register("mate", tools.cxBlend, alpha=0.2)
        self.toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=0.2, indpb=0.4)
        # self.toolbox.decorate("mutate", self.checkbounds())
        # self.toolbox.decorate("mate", self.checkbounds())
        self.toolbox.register("select", tools.selTournament, tournsize=3)
        self.toolbox.register("selectbest", tools.selBest)
        # self.toolbox.register("evaluate", self.evalmodel)

    def makepars(self):
        """Generates random list of parameters uniformly between defined bounds

        Notes
        -----
        This could be handled in a more sophisticated way, e.g. with a Gaussian distribution. However, a uniform random
        distribution is likely to be less biased.

        Returns
        -------
        pars: list
            A list of parameters within the specified ranges.
        """
        pars = [random.uniform(-1.0, 1.0) for _ in range(len(self.value_means))]
        return pars

    def parse_individual(self, individual):
        """converts an individual from the GA into a full list of parameters for building the specification object.
        Parameters
        ----------
        individual: DEAP individual

        Returns
        -------
        fullpars: list
            Full parameter list to define the specification object.
        """
        scaled_ind = []
        for i in range(len(self.value_means)):
            scaled_ind.append(self.value_means[i] + (individual[i] * self.value_ranges[i]))
        fullpars = list(self.arrangement)
        for k in range(len(self.keys)):
            for j in range(len(fullpars)):
                if fullpars[j] == self.keys[k]:
                    fullpars[j] = scaled_ind[k]
        return fullpars

    # def evalmodel(self, individual):
    #     """For a given individual builds the specification object and assesses the bude score
    #
    #     Parameters
    #     ----------
    #     individual: DEAP individual
    #         Individual representing a set of parameters to be evaluated
    #
    #     Returns
    #     -------
    #     model.bude_score: float
    #         The BUDE score for the evaluated model.
    #     """
    #     params = self.parse_individual(individual)
    #     model = self.specification(*params)
    #     model.build()
    #     model.pack_new_sequences(self.sequence)
    #     return model.bude_score

    def run_ga(self, gensize, numgen, procs, threshold=0.0, cxpb=0.5, mutpb=0.2, plot=False, log=False, run_id='model'):
        """Sets the GA running.

        Notes
        -----
        This is a modified version of example DEAP routine on github

        Parameters
        ----------
        gensize: int
            The number of individuals to be included in each generation. Twice this number is evaluated in the first
            instance. The actual number evaluated at each generation will be smaller by ca. 50% depending on the
            crossover and mutation probabilities.
        numgen: int
            The number of generations that the GA should run for.
        threshold: float
            A threshold value for terminating the minimization. If a value is supplied then the GA will terminate if
            five successive models are within the threshold value. Use with caution!
        cxpb: float
            A number between 1 and 0 specifying the probability that a pair of individuals will undergo a crossover
            event.
        mutpb: float
            A number between 1 and 0 specifying the probability that a given individual will be mutated. Note that once
            selected the likelihood of any particular parameter being mutated is given by a separate value owned by the
            mutation function itself and specified when the function is registered with the toolbox
        plot: Bool
            Logical flag specifying whether to plot the minimization profile.
        run_id: str
            Provides file name for log files and pdb model of best model
        """
        self.plot = plot
        self.log = log
        self.threshold = threshold
        self.gensize = gensize
        self.numgen = numgen
        self.procs = procs
        self.cxpb = cxpb
        if self.cxpb > 1 or self.cxpb < 0:
            raise ValueError("Crossover probability must be between 0 and 1")
        self.mutpb = mutpb
        if self.mutpb > 1 or self.mutpb < 0:
            raise ValueError("Crossover probability must be between 0 and 1")
        self.genbest = []
        self.run_id = run_id
        self.start_time = datetime.datetime.now()
        random.seed()
        self.halloffame = tools.HallOfFame(maxsize=1)
        self.stats = tools.Statistics(lambda ind: ind.fitness.values)
        self.stats.register("avg", numpy.mean)
        self.stats.register("std", numpy.std)
        self.stats.register("min", numpy.min)
        self.stats.register("max", numpy.max)
        self.logbook = tools.Logbook()
        self.logbook.header = ["gen", "evals"] + self.stats.fields
        self.pop = self.toolbox.population(n=self.gensize * 2)  # initial population is double normal size

        print ('Starting minimisation ({:%Y-%m-%d %H:%M:%S})'.format(self.start_time))
        px_parameters = zip([self.specification] * len(self.pop), [self.sequence] * len(self.pop),
                            [self.parse_individual(x) for x in self.pop], [self.bude_mode]*len(self.pop))
        with futures.ProcessPoolExecutor(max_workers=procs) as executor:
            fitnesses = executor.map(px_eval, px_parameters)
        for ind, fit in zip(self.pop, fitnesses):
            ind.fitness.values = (fit,)
        # initialmodels = []
        # for ind in self.pop:
        #     entry = [self.parse_individual(ind), ind.fitness.values[0]]
        #     initialmodels.append(entry)
        # initialmodels.sort(key=lambda thing: thing[1])
        self.halloffame.update(self.pop)
        self.modelcount = len(self.pop)

        for g in range(self.numgen):
            # Begin the evolution

            # Select the next generation individuals
            offspring = self.toolbox.select(self.pop, self.gensize)
            # Clone the selected individuals
            offspring = list(map(self.toolbox.clone, offspring))

            # Apply crossover and mutation on the offspring
            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < self.cxpb:
                    self.toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            for mutant in offspring:
                if random.random() < self.mutpb:
                    self.toolbox.mutate(mutant)
                    del mutant.fitness.values

            # Evaluate the individuals which have no fitness value
            without_fitness = [ind for ind in offspring if not ind.fitness.valid]

            valid_individuals = []
            invalid_individuals = []
            for ind in without_fitness:
                if any(x > 1 or x < -1 for x in ind):
                    invalid_individuals.append(ind)
                else:
                    valid_individuals.append(ind)

            # Assign dummy fitness to out of bounds individuals
            for ind in invalid_individuals:
                excess = []
                for i in ind:
                    if abs(i) > 1:
                        excess.append(i)
                ind.fitness.values = (10*sum(excess)**2,)

            # Evaluate the valid individuals
            self.modelcount += len(valid_individuals)
            topologies = [self.specification]*len(valid_individuals)
            sequences = [self.sequence]*len(valid_individuals)
            parameters = [self.parse_individual(x) for x in valid_individuals]
            mode = [self.bude_mode]*len(valid_individuals)
            px_parameters = zip(topologies, sequences, parameters, mode)
            with futures.ProcessPoolExecutor(max_workers=procs) as executor:
                fitnesses = executor.map(px_eval, px_parameters)
            for ind, fit in zip(valid_individuals, fitnesses):
                ind.fitness.values = (fit,)

            # The population is entirely replaced by the offspring
            self.pop[:] = valid_individuals + invalid_individuals
            self.halloffame.update(self.pop)
            self.logbook.record(gen=g, evals=len(valid_individuals), **self.stats.compile(self.pop))
            print(self.logbook.stream)

            if self.threshold and (g > 10):
                if self.threshold >= numpy.ptp(self.genbest[-5:]) > 0.0:
                    print('Difference =  {0}'.format(numpy.ptp(self.genbest[-5:])))
                    print("-- Threshold conditions met, stopping minimization")
                    print("-- Evaluated {0} models in total --".format(self.modelcount))
                    break
        self.end_time = datetime.datetime.now()
        time_taken = self.end_time - self.start_time
        print("End of minimisation ({:%Y-%m-%d %H:%M:%S})".format(self.end_time))
        print("Run ID is {0}".format(self.run_id))
        print('Minimization time = {0}'.format(time_taken))
        print("Evaluated {0} models in total".format(self.modelcount))
        print("Best score is {0}".format(self.halloffame[0].fitness.values[0]))
        print("Best parameters are {0}".format(self.parse_individual(self.halloffame[0])))
        for i, entry in enumerate(self.halloffame[0]):
            if entry > 0.95:
                print("Warning! Parameter {0} is at or near maximum allowed value\n".format(i+1))
            elif entry < -0.95:
                print("Warning! Parameter {0} is at or near minimum allowed value".format(i+1))
        if self.log:
            self.log_results()
        if self.plot:
            print('----Minimisation plot:')
            plt.figure(figsize=(5, 5))
            plt.plot(range(len(self.logbook.select('min'))), self.logbook.select('min'))
            plt.xlabel('Iteration', fontsize=20)
            plt.ylabel('Score', fontsize=20)

    def log_results(self):
        """Saves files for the minimization. Currently saves a logfile with best individual and a pdb of the best model
        """
        best_ind = self.halloffame[0]
        params = self.parse_individual(best_ind)
        with open('{0}{1}_log.txt'.format(self.output_path, self.run_id), 'a+') as log_file:
            log_file.write('\nEvaluated {0} models in total\n'.format(self.modelcount))
            log_file.write('Run ID is {0}\n'.format(self.run_id))
            log_file.write('Best score is {0}\n'.format(self.halloffame[0].fitness.values[0]))
            log_file.write('Parameters of best model are {0}'.format(params))
            log_file.write('Best individual is {0}\n'.format(self.halloffame[0]))
            for i, entry in enumerate(self.halloffame[0]):
                if entry > 0.95:
                    log_file.write("Warning! Parameter {0} is at or near maximum allowed value\n".format(i+1))
                elif entry < -0.95:
                    log_file.write("Warning! Parameter {0} is at or near minimum allowed value\n".format(i+1))
            log_file.write('Minimization history: \n{0}'.format(self.logbook))
        with open('{0}{1}_bestmodel.pdb'.format(self.output_path, self.run_id), 'w') as output_file:
            model = self.specification(*params)
            model.build()
            model.pack_new_sequences(self.sequence)
            output_file.write(model.pdb)


def px_eval(model_specifications):
    """Special version of the eval function that makes no reference to .self and so can be farmed out to other
    namespaces.

    Notes
    -----
    This is set up for parallel execution so makes no reference to local namespace

    Parameters
    ----------
    model_specifications: list
        Tuple containing the specification to be built, the sequence, and the parameters for model building.

    Returns
    -------
    model.bude_score: float
        BUDE score for model to be assigned to particle fitness value.
    """
    specification, sequence, parsed_ind, mode = model_specifications
    model = specification(*parsed_ind)
    model.build()
    model.pack_new_sequences(sequence)
    if mode == 'average':
        return model.average_bude_score
    elif mode == 'additive':
        return model.bude_score
    elif mode == 'internal':
        return model.bude_internal_energy

__author__ = 'Andrew R. Thomson'
__status__ = 'Development'
