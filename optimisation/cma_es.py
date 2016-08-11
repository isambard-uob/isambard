from deap import algorithms, base, creator, tools
# import random
import numpy
import matplotlib.pylab as plt
import datetime
from concurrent import futures

# TODO: Tidy and naming... Make initial population random uniform within bounds?


class OptCMAES:
    """Initializes the CMA-ES class to run a covariance matrix adaptation evolutionary strategy minimization

    Parameters
    ----------
    specification: Isambard specification
        Tells the optimizer what kind of structure it is building
    output_path: str
        Specifies the location to save files if logging is enabled
    bude_mode: 'average' or 'additive'
        Specifies whether the score should be the average per chain or total cumulative score.
    """
    def __init__(self, specification, output_path, bude_mode='average'):
        self.output_path = output_path
        self.bude_mode = bude_mode
        self.sequence = None
        self.minvalues = None
        self.maxvalues = None
        self.value_ranges = None
        self.value_means = None
        self.arrangement = None
        self.specification = specification
        self.keys = None
        self.log = None
        self.modelcount = None
        self.population = None
        self.plot = None
        self.threshold = None
        self.gensize = None
        self.numgen = None
        self.log = None
        self.genbest = None
        self.run_id = None
        self.start_time = None
        self.procs = None
        self.logbook = None

        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMin)
        self.toolbox = base.Toolbox()
        self.toolbox.register("evaluate", px_eval)

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
        self.minvalues = [a - b for a, b in zip(self.value_means, self.value_ranges)]
        self.maxvalues = [a + b for a, b in zip(self.value_means, self.value_ranges)]
        self.arrangement = arrangement
        self.keys = []
        for i in range(len(self.maxvalues)):
            self.keys.append("".join(['var', str(i)]))
        if len(set(arrangement).intersection(self.keys)) != len(self.minvalues):
            raise ValueError("argument mismatch!")
        if len(self.maxvalues) != len(self.minvalues):
            raise ValueError("argument mismatch!")

    def make_individual(self, paramlist):
        part = creator.Individual(paramlist)
        part.index = None
        return part

    def run_cmaes(self, lambda_, numgen, procs, sigma=0.3, weights='superlinear', run_id='CMA_ES_model', plot=False,
                  log=False):
        """Starts the CMA-ES running.

        Parameters
        ----------
        numgen: int
            Number of iterations of the CMA to run. A value of 50 is typical.
        lambda_: int
            Size of population to use. Values of 20-50 are common depending on the number of paramters to optimize
        procs: int
            The number of processors on which the strategy will run
        sigma: float
            Determines the initial distribution of the population. Larger values scatter them more widely.
        weights: 'superlinear', 'linear' or 'equal'
            Determine the rate at which the population contracts
        run_id: str
            Provides file name for logging purposes
        plot: bool
            Whether or not to plot the minimization graph
        log: bool
            Whether or not to save files for the best model and minimization history.
        """
        # The cma module uses the numpy random number generator
        numpy.random.seed()
        self.plot = plot
        self.log = log
        N = len(self.value_means)
        self.initialize_cma_es(sigma=sigma, weights=weights, lambda_=lambda_, centroid=[0]*N)  # add some kwargs from args here?
        # self.toolbox.register("individual", creator.Individual)
        self.toolbox.register("individual", self.make_individual)
        self.toolbox.register("generate", self.generate, self.toolbox.individual)
        self.halloffame = tools.HallOfFame(maxsize=1)
        self.toolbox.register("update", self.update)
        self.stats = tools.Statistics(lambda ind: ind.fitness.values)
        self.stats.register("avg", numpy.mean)
        self.stats.register("std", numpy.std)
        self.stats.register("min", numpy.min)
        self.stats.register("max", numpy.max)
        self.logbook = tools.Logbook()
        self.logbook.header = ["gen", "evals"] + self.stats.fields
        self.modelcount = 0
        self.numgen = numgen
        self.procs = procs
        self.run_id = run_id
        self.start_time = datetime.datetime.now()
        print ('Starting minimisation ({:%Y-%m-%d %H:%M:%S})'.format(self.start_time))
        for gen in range(self.numgen):
            # Generate a new population
            self.toolbox.generate()
            for i in range(len(self.population)):
                self.population[i].index = i
            # separate out constraint violating individuals
            valid_individuals = []
            invalid_individuals = []
            for ind in self.population:
                if any(x > 1 or x < -1 for x in ind):
                    invalid_individuals.append(ind)
                else:
                    valid_individuals.append(ind)

            # Assign dummy fitness to invalid individuals
            for ind in invalid_individuals:
                excess = []
                for i in ind:
                    if abs(i) > 1:
                        excess.append(i)
                ind.fitness.values = (10+10*sum(excess)**2,)

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
            self.population[:] = valid_individuals + invalid_individuals
            self.population.sort(key=lambda _: _.index)
            self.halloffame.update(self.population)

            # Update the strategy with the evaluated individuals
            self.toolbox.update(self.population)
            record = self.stats.compile(self.population)
            self.logbook.record(gen=gen, evals=len(valid_individuals), **record)
            print(self.logbook.stream)
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

    def parse_individual(self, individual):
        """Converts an individual from the CMA ES into a full list of parameters for building the specification object.
        Parameters
        ----------
        individual: position elements from swarm particle

        Returns
        -------
        fullpars: list
            Full parameter list for model building.
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

    def evalmodel(self, individual):
        """For a given individual builds the specification object and assesses the bude score
        """
        params = self.parse_individual(individual)
        model = self.specification(*params)
        model.evaluate(self.sequence)
        return model.bude_score

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

    def initialize_cma_es(self, **kwargs):
        """
        A strategy that will keep track of the basic parameters of the CMA-ES
        algorithm.
        :param centroid: An iterable object that indicates where to start the
                         evolution.
        :param sigma: The initial standard deviation of the distribution.
        :param parameter: One or more parameter to pass to the strategy as
                          described in the following table, optional.
        +----------------+---------------------------+----------------------------+
        | Parameter      | Default                   | Details                    |
        +================+===========================+============================+
        | ``lambda_``    | ``int(4 + 3 * log(N))``   | Number of children to      |
        |                |                           | produce at each generation,|
        |                |                           | ``N`` is the individual's  |
        |                |                           | size (integer).            |
        +----------------+---------------------------+----------------------------+
        | ``mu``         | ``int(lambda_ / 2)``      | The number of parents to   |
        |                |                           | keep from the              |
        |                |                           | lambda children (integer). |
        +----------------+---------------------------+----------------------------+
        | ``cmatrix``    | ``identity(N)``           | The initial covariance     |
        |                |                           | matrix of the distribution |
        |                |                           | that will be sampled.      |
        +----------------+---------------------------+----------------------------+
        | ``weights``    | ``"superlinear"``         | Decrease speed, can be     |
        |                |                           | ``"superlinear"``,         |
        |                |                           | ``"linear"`` or            |
        |                |                           | ``"equal"``.               |
        +----------------+---------------------------+----------------------------+
        | ``cs``         | ``(mueff + 2) /           | Cumulation constant for    |
        |                | (N + mueff + 3)``         | step-size.                 |
        +----------------+---------------------------+----------------------------+
        | ``damps``      | ``1 + 2 * max(0, sqrt((   | Damping for step-size.     |
        |                | mueff - 1) / (N + 1)) - 1)|                            |
        |                | + cs``                    |                            |
        +----------------+---------------------------+----------------------------+
        | ``ccum``       | ``4 / (N + 4)``           | Cumulation constant for    |
        |                |                           | covariance matrix.         |
        +----------------+---------------------------+----------------------------+
        | ``ccov1``      | ``2 / ((N + 1.3)^2 +      | Learning rate for rank-one |
        |                | mueff)``                  | update.                    |
        +----------------+---------------------------+----------------------------+
        | ``ccovmu``     | ``2 * (mueff - 2 + 1 /    | Learning rate for rank-mu  |
        |                | mueff) / ((N + 2)^2 +     | update.                    |
        |                | mueff)``                  |                            |
        +----------------+---------------------------+----------------------------+
        """
        self.params = kwargs

        # Create a centroid as a numpy array
        self.centroid = numpy.array([0]*len(self.value_means))

        self.dim = len(self.centroid)
        self.sigma = self.params.get("sigma", 0.5)
        self.pc = numpy.zeros(self.dim)
        self.ps = numpy.zeros(self.dim)
        self.chiN = numpy.sqrt(self.dim) * (1 - 1. / (4. * self.dim) +
                                      1. / (21. * self.dim ** 2))

        self.C = self.params.get("cmatrix", numpy.identity(self.dim))
        self.diagD, self.B = numpy.linalg.eigh(self.C)

        indx = numpy.argsort(self.diagD)
        self.diagD = self.diagD[indx] ** 0.5
        self.B = self.B[:, indx]
        self.BD = self.B * self.diagD

        self.cond = self.diagD[indx[-1]] / self.diagD[indx[0]]

        self.lambda_ = self.params.get("lambda_", int(4 + 3 * numpy.log(self.dim)))
        self.update_count = 0
        self.computeParams(self.params)

    def generate(self, func):
        """Generate a population of :math:`\lambda` individuals of type
        *ind_init* from the current strategy.
        :param ind_init: A function object that is able to initialize an
                         individual from a list.
        :returns: A list of individuals.
        """
        arz = numpy.random.standard_normal((self.lambda_, self.dim))
        arz = self.centroid + self.sigma * numpy.dot(arz, self.BD.T)
        self.population = list(map(func, arz))

    def update(self, population):
        """Update the current covariance matrix strategy from the
        *population*.
        :param population: A list of individuals from which to update the
                           parameters.
        """
        population.sort(key=lambda ind: ind.fitness, reverse=True)

        old_centroid = self.centroid
        self.centroid = numpy.dot(self.weights, population[0:self.mu])

        c_diff = self.centroid - old_centroid

        # Cumulation : update evolution path
        self.ps = (1 - self.cs) * self.ps \
            + numpy.sqrt(self.cs * (2 - self.cs) * self.mueff) / self.sigma \
            * numpy.dot(self.B, (1. / self.diagD)
                        * numpy.dot(self.B.T, c_diff))

        hsig = float((numpy.linalg.norm(self.ps) /
                      numpy.sqrt(1. - (1. - self.cs) ** (2. * (self.update_count + 1.))) / self.chiN
                      < (1.4 + 2. / (self.dim + 1.))))

        self.update_count += 1

        self.pc = (1 - self.cc) * self.pc + hsig \
            * numpy.sqrt(self.cc * (2 - self.cc) * self.mueff) / self.sigma \
            * c_diff

        # Update covariance matrix
        artmp = population[0:self.mu] - old_centroid
        self.C = (1 - self.ccov1 - self.ccovmu + (1 - hsig)
                  * self.ccov1 * self.cc * (2 - self.cc)) * self.C \
            + self.ccov1 * numpy.outer(self.pc, self.pc) \
            + self.ccovmu * numpy.dot((self.weights * artmp.T), artmp) \
            / self.sigma ** 2

        self.sigma *= numpy.exp((numpy.linalg.norm(self.ps) / self.chiN - 1.)
                                * self.cs / self.damps)

        self.diagD, self.B = numpy.linalg.eigh(self.C)
        indx = numpy.argsort(self.diagD)

        self.cond = self.diagD[indx[-1]] / self.diagD[indx[0]]

        self.diagD = self.diagD[indx] ** 0.5
        self.B = self.B[:, indx]
        self.BD = self.B * self.diagD

    def computeParams(self, params):
        """Computes the parameters depending on :math:`\lambda`. It needs to
        be called again if :math:`\lambda` changes during evolution.
        :param params: A dictionary of the manually set parameters.
        """
        self.mu = params.get("mu", int(self.lambda_ / 2))
        rweights = params.get("weights", "superlinear")
        if rweights == "superlinear":
            self.weights = numpy.log(self.mu + 0.5) - \
                numpy.log(numpy.arange(1, self.mu + 1))
        elif rweights == "linear":
            self.weights = self.mu + 0.5 - numpy.arange(1, self.mu + 1)
        elif rweights == "equal":
            self.weights = numpy.ones(self.mu)
        else:
            raise RuntimeError("Unknown weights : %s" % rweights)

        self.weights /= sum(self.weights)
        self.mueff = 1. / sum(self.weights ** 2)

        self.cc = params.get("ccum", 4. / (self.dim + 4.))
        self.cs = params.get("cs", (self.mueff + 2.) /
                                   (self.dim + self.mueff + 3.))
        self.ccov1 = params.get("ccov1", 2. / ((self.dim + 1.3) ** 2 +
                                               self.mueff))
        self.ccovmu = params.get("ccovmu", 2. * (self.mueff - 2. +
                                                 1. / self.mueff) /
                                 ((self.dim + 2.) ** 2 + self.mueff))
        self.ccovmu = min(1 - self.ccov1, self.ccovmu)
        self.damps = 1. + 2. * max(0, numpy.sqrt((self.mueff - 1.) /
                                           (self.dim + 1.)) - 1.) + self.cs
        self.damps = params.get("damps", self.damps)


def px_eval(model_specifications):
    """Special version of the eval function that makes no reference to .self and so can be farmed out to other
    namespaces.

    Parameters
    ----------
    model_specifications: list
        Provided by the run_cma_es_px method, contains reference to the specification, the sequence, and the parameters for
        the model, and the bude mode to specify average versus cumulative score.
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
