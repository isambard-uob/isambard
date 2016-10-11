from concurrent import futures
import datetime
import operator
import random
import sys

from deap import base, creator, tools
import numpy
import matplotlib.pylab as plt

from external_programs.profit import run_profit


class BaseOptimizer:

    def __init__(self, **kwargs):
        self._params = {}
        self._params.update(**kwargs)
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        self.toolbox = base.Toolbox()

    def parse_individual(self, individual):
        """Converts a deap individual into a full list of parameters for building the specification object.
        Parameters
        ----------
        individual: deap individual from optimization. Details vary according to type of optimization, but parameters
        within deap individual are always between -1 and 1. This function converts them into the values used to
        actually build the model

        Returns
        -------
        fullpars: list
            Full parameter list for model building.
        """
        scaled_ind = []
        for i in range(len(self._params['value_means'])):
            scaled_ind.append(self._params['value_means'][i] + (individual[i] * self._params['value_ranges'][i]))
        fullpars = list(self._params['arrangement'])
        for k in range(len(self._params['variable_parameters'])):
            for j in range(len(fullpars)):
                if fullpars[j] == self._params['variable_parameters'][k]:
                    fullpars[j] = scaled_ind[k]
        return fullpars

    def run_opt(self, popsize, numgen, processors, plot=False, log=False, **kwargs):
        """
        Runs the optimizer.
        :param popsize:
        :param numgen:
        :param processors:
        :param plot:
        :param log:
        :param kwargs:
        :return:
        """
        self._params['popsize'] = popsize
        self._params['numgen'] = numgen
        self._params['processors'] = processors
        self._params['plot'] = plot
        self._params['log'] = log
        self._params.update(**kwargs)  # allows us to pass in additional arguments e.g. neighbours
        self.halloffame = tools.HallOfFame(1)
        self.stats = tools.Statistics(lambda thing: thing.fitness.values)
        self.stats.register("avg", numpy.mean)
        self.stats.register("std", numpy.std)
        self.stats.register("min", numpy.min)
        self.stats.register("max", numpy.max)
        self.logbook = tools.Logbook()
        self.logbook.header = ["gen", "evals"] + self.stats.fields
        self._params['model_count'] = 0
        start_time = datetime.datetime.now()
        self.initialize_pop()
        for g in range(self._params['numgen']):
            self.update_pop()
            self.halloffame.update(self.population)
            self.logbook.record(gen=g, evals=self._params['evals'], **self.stats.compile(self.population))
            print(self.logbook.stream)
        end_time = datetime.datetime.now()
        time_taken = end_time - start_time
        self._params['time_taken'] = time_taken
        print("Evaluated {0} models in total".format(self._params['model_count']))
        print("Best fitness is {0}".format(self.halloffame[0].fitness))
        print("Best parameters are {0}".format(self.parse_individual(self.halloffame[0])))
        for i, entry in enumerate(self.halloffame[0]):
            if entry > 0.95:
                print("Warning! Parameter {0} is at or near maximum allowed value\n".format(i+1))
            elif entry < -0.95:
                print("Warning! Parameter {0} is at or near minimum allowed value\n".format(i+1))
        if self._params['log']:
            self.log_results()
        if self._params['plot']:
            print('----Minimisation plot:')
            plt.figure(figsize=(5, 5))
            plt.plot(range(len(self.logbook.select('min'))), self.logbook.select('min'))
            plt.xlabel('Iteration', fontsize=20)
            plt.ylabel('Score', fontsize=20)

    def parameters(self, sequence, value_means, value_ranges, arrangement):
        """Relates the individual to be evolved to the full parameter string for building the topology object

        Parameters
        ----------
        sequence: str
            Full amino acid sequence for topology object to be optimized. Must be equal to the number of residues in the
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
        self._params['sequence'] = sequence
        self._params['value_means'] = value_means
        self._params['value_ranges'] = value_ranges
        self._params['arrangement'] = arrangement
        if any(x <= 0 for x in self._params['value_ranges']):
            raise ValueError("range values must be greater than zero")
        self._params['variable_parameters'] = []
        for i in range(len(self._params['value_means'])):
            self._params['variable_parameters'].append("".join(['var', str(i)]))
        if len(set(arrangement).intersection(self._params['variable_parameters'])) != \
                len(self._params['value_means']):
            raise ValueError("argument mismatch!")
        if len(self._params['value_ranges']) != len(self._params['value_means']):
            raise ValueError("argument mismatch!")

    def assign_fitnesses(self):
        raise NotImplementedError("Will depend on evaluation subclass")
        # must always operate on whole population- bounds checking etc to be done internally

    def generate(self):
        raise NotImplementedError("Will depend on optimizer type")

    def initialize_pop(self):
        raise NotImplementedError("Will depend on optimizer type")

    def update_pop(self):
        raise NotImplementedError("Will depend on optimizer type")

    def log_results(self):
        """Saves files for the minimization. Currently saves a logfile with best individual and a pdb of the best model
        """
        best_ind = self.halloffame[0]
        model_params = self.parse_individual(best_ind)  # need to change name of 'params'
        with open('{0}{1}_log.txt'.format(self._params['output_path'], self._params['run_id']), 'a+') as log_file:
            log_file.write('\nEvaluated {0} models in total\n'.format(self._params['model_count']))
            log_file.write('Run ID is {0}\n'.format(self._params['run_id']))
            log_file.write('Best fitness is {0}\n'.format(self.halloffame[0].fitness))
            log_file.write('Parameters of best model are {0}\n'.format(model_params))
            log_file.write('Best individual is {0}\n'.format(self.halloffame[0]))
            for i, entry in enumerate(self.halloffame[0]):
                if entry > 0.95:
                    log_file.write("Warning! Parameter {0} is at or near maximum allowed value\n".format(i+1))
                elif entry < -0.95:
                    log_file.write("Warning! Parameter {0} is at or near minimum allowed value\n".format(i+1))
            log_file.write('Minimization history: \n{0}'.format(self.logbook))
        with open('{0}{1}_bestmodel.pdb'.format(self._params['output_path'], self._params['run_id']),
                  'w') as output_file:
            model = self._params['topology'](*model_params)
            model.build()
            model.pack_new_sequences(self._params['sequence'])
            output_file.write(model.pdb)


class BaseScore(BaseOptimizer):
    """
    Assigns BUFF score as fitness to individuals in optimization
    """
    def __init__(self):
        super().__init__()

    def assign_fitnesses(self, targets):
        self._params['evals'] = len(targets)
        px_parameters = zip([self._params['topology']] * len(targets),
                            [self._params['sequence']] * len(targets),
                            [self.parse_individual(x) for x in targets])
        if (self._params['processors'] == 1) or (sys.platform == 'win32'):
            fitnesses = map(buff_eval, px_parameters)
        else:
            with futures.ProcessPoolExecutor(max_workers=self._params['processors']) as executor:
                fitnesses = executor.map(buff_eval, px_parameters)
        for ind, fit in zip(targets, fitnesses):
            ind.fitness.values = (fit,)


class BaseRMSD(BaseOptimizer):
    """
    Assigns RMSD as fitness to individuals in optimization. Allows optimization of parameters to best fit a target model
    """
    def __init__(self):
        super().__init__()

    def assign_fitnesses(self, targets):
        self._params['evals'] = len(targets)
        px_parameters = zip([self._params['topology']] * len(targets),
                            [self._params['sequence']] * len(targets),
                            [self.parse_individual(x) for x in targets],
                            [self._params['ref_pdb']] * len(targets))
        if (self._params['processors'] == 1) or (sys.platform == 'win32'):
            fitnesses = map(rmsd_eval, px_parameters)
        else:
            with futures.ProcessPoolExecutor(max_workers=self._params['processors']) as executor:
                fitnesses = executor.map(rmsd_eval, px_parameters)
        for ind, fit in zip(targets, fitnesses):
            ind.fitness.values = (fit,)


class BaseComparator(BaseOptimizer):
    """
    Assigns individual fitness to be change in BUFF score on positioning two predefined models according to parameters
     from individual. Allows basic rigid body docking between two AMPAL objects with side chain repacking in order to
     estimate interactions.
    """
    def __init__(self):
        super().__init__()

    def assign_fitnesses(self, targets):
        self._params['evals'] = len(targets)
        px_parameters = zip([self._params['top1']] * len(targets),
                            [self._params['top2']] * len(targets),
                            [self._params['params1']] * len(targets),
                            [self._params['params2']] * len(targets),
                            [self._params['seq1']] * len(targets),
                            [self._params['seq2']] * len(targets),
                            [self.parse_individual(x) for x in targets])
        if (self._params['processors'] == 1) or (sys.platform == 'win32'):
            fitnesses = map(comparator_eval, px_parameters)
        else:
            with futures.ProcessPoolExecutor(max_workers=self._params['processors']) as executor:
                fitnesses = executor.map(comparator_eval, px_parameters)
        for ind, fit in zip(targets, fitnesses):
            ind.fitness.values = (fit - (self._params['ref1'] + self._params['ref2']),)

    def parameters(self, value_means, value_ranges, arrangement):
        """Relates the individual to be evolved to the full parameter string for building the topology object.
        Special version for comparator type optimizers that doesn't require sequence. Should take up to six parameters
        defining the x, y, z rotations and x, y, z translations in that order. For example testing rotation of 60 +/- 20
        degrees around the z axis at a displacement of 20 +/- 10 Angstrom would require:
        value_means = [60, 20], value_ranges = [20, 10], arrangement = ['var0', 0, 0, 'var1', 0, 0]

        Parameters
        ----------
        value_means: list
            List containing mean values for parameters to be optimized.
        value_ranges: list
            List containing ranges for parameters to be optimized. Values must be positive.
        arrangement: list
            Full list of fixed and variable parameters for model building. Fixed values are the appropriate value.
            Values to be varied should be listed as 'var0', 'var1' etc, and must be in ascending numerical order.
            Variables can be repeated if required.
        """
        self._params['value_means'] = value_means
        self._params['value_ranges'] = value_ranges
        self._params['arrangement'] = arrangement
        if any(x <= 0 for x in self._params['value_ranges']):
            raise ValueError("range values must be greater than zero")
        self._params['variable_parameters'] = []
        for i in range(len(self._params['value_means'])):
            self._params['variable_parameters'].append("".join(['var', str(i)]))
        if len(set(arrangement).intersection(self._params['variable_parameters'])) != \
                len(self._params['value_means']):
            raise ValueError("argument mismatch!")
        if len(self._params['value_ranges']) != len(self._params['value_means']):
            raise ValueError("argument mismatch!")


class OptDE:
    """
    A differential evolution algorithm. Can use neighbourhood model to reduce chance of getting stuck in local optima.
    This is a very versatile algorithm, and its use is recommended.
    """
    def __init__(self, **kwargs):
        super().__init__()
        self._params.update(**kwargs)
        self._params.setdefault('cxpb', 0.75)
        self._params.setdefault('diff_weight', 1)
        self._params.setdefault('output_path', None)
        self._params.setdefault('neighbours', None)
        creator.create("Individual", list, fitness=creator.FitnessMin)

    def generate(self):
        """Generates a particle using the creator function. Position and speed are uniformly randomly seeded within
        allowed bounds. The particle also has speed limit settings taken from global values.

        Returns
        -------
        particle object
        """
        ind = creator.Individual([random.uniform(-1, 1) for _ in range(len(self._params['value_means']))])
        ind.ident = None
        ind.neighbours = None
        return ind

    def initialize_pop(self):
        """Assigns indices to individuals in population where neighbourhood model is used and assigns initial fitness"""
        self.toolbox.register("individual", self.generate)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.population = self.toolbox.population(n=self._params['popsize'])
        if self._params['neighbours']:
            for i in range(len(self.population)):
                self.population[i].ident = i
                self.population[i].neighbours = list(set([(i-x) % len(self.population)
                                                   for x in range(1, self._params['neighbours']+1)] +
                                                  [(i+x) % len(self.population)
                                                   for x in range(1, self._params['neighbours']+1)]))
        self.assign_fitnesses(self.population)

    def crossover(self, ind):
        """Used by the evolution process to generate a new individual. This is a tweaked version of the classical DE
        crossover algorithm, the main difference that candidate parameters are generated using a lognormal distribution.
        Bound handling is achieved by resampling where the candidate solution exceeds +/-1

        Parameters
        ----------

        Returns
        -------
        y: deap individual
            An individual representing a candidate solution, to be assigned a fitness.
        """
        if self._params['neighbours']:
            a, b, c = random.sample([self.population[i] for i in ind.neighbours], 3)
        else:
            a, b, c = random.sample(self.population, 3)
        y = self.toolbox.clone(a)
        y.ident = ind.ident
        y.neighbours = ind.neighbours
        del y.fitness.values
        # y should now be a copy of ind with the vector elements from a
        ident = random.randrange(len(self._params['value_means']))
        for i, value in enumerate(y):
            if i == ident or random.random() < self._params['cxpb']:
                entry = a[i] + random.lognormvariate(-1.2, 0.5)*self._params['diff_weight']*(b[i]-c[i])
                tries = 0
                while abs(entry) > 1.0:
                    tries += 1
                    entry = a[i] + random.lognormvariate(-1.2, 0.5)*self._params['diff_weight']*(b[i]-c[i])
                    if tries > 10000:
                        entry = a[i]
                y[i] = entry
        return y

    def update_pop(self):
        """
        Updates the population according to crossover and fitness criteria
        :return:
        """
        candidates = []
        for ind in self.population:
            candidates.append(self.crossover(ind))
        self._params['model_count'] += len(candidates)
        self.assign_fitnesses(candidates)
        for i in range(len(self.population)):
            if candidates[i].fitness > self.population[i].fitness:
                self.population[i] = candidates[i]


class OptPSO:
    """
    A particle swarm optimization algorithm, using the constriction factor method.
    This is good for avoiding bias and premature minimization, though it may struggle to find the ultimate optimum
    solution. Supports the neighbourhood model. Bound handling is achieved by allowing particles to exceed permitted
    bounds, but not assigning them a fitness in this case.
    """
    def __init__(self, **kwargs):
        self.population = None
        super().__init__()
        self._params.update(**kwargs)
        self._params.setdefault('output_path', None)
        self._params.setdefault('max_speed', 0.75)
        self._params.setdefault('neighbours', None)
        creator.create("Particle", list, fitness=creator.FitnessMin, speed=list, smin=None, smax=None, best=None)
        self.toolbox.register("particle", self.generate)
        creator.create("Swarm", list, gbest=None, gbestfit=creator.FitnessMin)  # can this pick up the global fitness?
        self.toolbox.register("swarm", tools.initRepeat, creator.Swarm, self.toolbox.particle)

    def initialize_pop(self):
        """
        Generates initial population with random positions and speeds.
        :return:
        """
        self.population = self.toolbox.swarm(n=self._params['popsize'])
        if self._params['neighbours']:
            for i in range(len(self.population)):
                self.population[i].ident = i
                self.population[i].neighbours = list(set([(i - x) % len(self.population)
                                                          for x in range(1, self._params['neighbours']+1)] + [i] +
                                                         [(i+x) % len(self.population)
                                                   for x in range(1, self._params['neighbours']+1)]))
        else:
            for i in range(len(self.population)):
                self.population[i].ident = i
                self.population[i].neighbours = [x for x in range(len(self.population))]
        self.assign_fitnesses(self.population)
        for part in self.population:
            part.best = creator.Particle(part)
            part.best.fitness.values = part.fitness.values
        # self.pop.gbestfit = max(part.fitness for part in self.pop)
        # self.pop.gbest = max(enumerate(self.pop), key=lambda x: self.pop[x[0]].fitness)[1]

    def generate(self):
        """Generates a particle using the creator function. Position and speed are uniformly randomly seeded within
        allowed bounds. The particle also has speed limit settings taken from global values.

        Returns
        -------
        particle object
        """
        part = creator.Particle([random.uniform(-1, 1) for _ in range(len(self._params['value_means']))])
        part.speed = [random.uniform(-self._params['max_speed'], self._params['max_speed'])
                      for _ in range(len(self._params['value_means']))]
        part.smin = -self._params['max_speed']
        part.smax = self._params['max_speed']
        part.ident = None
        part.neighbours = None
        return part

    def update_particle(self, part, chi=0.729843788, c=2.05):
        """Constriction factor update particle method that looks for a list of neighbours attached to a particle and
        uses the particle's best position and that of the best neighbour.
        """
        neighbour_pool = [self.population[i] for i in part.neighbours]
        best_neighbour = max(neighbour_pool, key=lambda x: x.best.fitness)
        ce1 = (c * random.uniform(0, 1) for _ in range(len(part)))
        ce2 = (c * random.uniform(0, 1) for _ in range(len(part)))
        ce1_p = map(operator.mul, ce1, map(operator.sub, part.best, part))
        ce2_g = map(operator.mul, ce2, map(operator.sub, best_neighbour.best, part))
        chi_list = [chi]*len(part)
        chi_list2 = [1 - chi]*len(part)
        a = map(operator.sub, map(operator.mul, chi_list, map(operator.add, ce1_p, ce2_g)),
                map(operator.mul, chi_list2, part.speed))
        part.speed = list(map(operator.add, part.speed, a))
        for i, speed in enumerate(part.speed):
            if speed < part.smin:
                part.speed[i] = part.smin
            elif speed > part.smax:
                part.speed[i] = part.smax
        part[:] = list(map(operator.add, part, part.speed))

    def update_pop(self):
        """
        Assigns fitnesses to particles that are within bounds.
        :return:
        """
        valid_particles = []
        invalid_particles = []
        for part in self.population:
            if any(x > 1 or x < -1 for x in part):
                invalid_particles.append(part)
            else:
                valid_particles.append(part)
        self._params['model_count'] += len(valid_particles)
        for part in valid_particles:
            self.update_particle(part)
        self.assign_fitnesses(valid_particles)
        for part in valid_particles:
            if part.fitness > part.best.fitness:
                part.best = creator.Particle(part)
                part.best.fitness = part.fitness
        for part in invalid_particles:
            self.update_particle(part)
        self.population[:] = valid_particles + invalid_particles
        self.population.sort(key=lambda x: x.ident)  # shouldn't need to sort?
        # for part in self.population:
        #     if part.best.fitness < part.fitness:
        #         part.best = creator.Particle(part)
        #         part.best.fitness.values = part.fitness.values
        # self.pop.gbestfit = max(part.fitness for part in self.pop) #this is the current best, not the all time best
        # self.pop.gbest = max(enumerate(self.pop), key=lambda x: self.pop[x[0]].fitness)[1]  #but these aren't used anyway


class OptGA:
    """
    A classic genetic algorithm optimization algorithm. Arguably the weakest of the algorithms available, but very good
    for eliminating unfavourable regions of the search space. Can be heavily customized in terms of mutation and
    crossover operators etc. Bound handling is achieved simply by amending any out of bounds parameters to the boundary
    value.
    """
    def __init__(self, **kwargs):
        super().__init__()
        self._params.update(**kwargs)
        self._params.setdefault('output_path', None)
        self._params.setdefault('cxpb', 0.5)
        self._params.setdefault('mutpb', 0.2)
        creator.create("Individual", list, fitness=creator.FitnessMin)
        self.toolbox.register("mate", tools.cxBlend, alpha=0.2)
        self.toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=0.2, indpb=0.4)
        self.toolbox.register("select", tools.selTournament)

    def generate(self):
        """
        Generates an individual with random parameters within allowed bounds
        """
        ind = creator.Individual([random.uniform(-1, 1) for _ in range(len(self._params['value_means']))])
        return ind

    def initialize_pop(self):
        """
        Assigns initial fitnesses
        :return:
        """
        self.toolbox.register("individual", self.generate)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.population = self.toolbox.population(n=self._params['popsize'])
        self.assign_fitnesses(self.population)
        self._params['model_count'] += len(self.population)

    def update_pop(self):
        offspring = list(map(self.toolbox.clone, self.population))
        # offspring.sort(reverse=True, key=lambda x: x.fitness)

        for _ in range(self._params['popsize']//2):
            if random.random() < self._params['cxpb']:
                child1, child2 = self.toolbox.select(self.population, 2, 6)
                temp1 = self.toolbox.clone(child1)
                temp2 = self.toolbox.clone(child2)
                self.toolbox.mate(temp1, temp2)
                del temp1.fitness.values
                del temp2.fitness.values
                offspring.append(temp1)
                offspring.append(temp2)

        # for child1, child2 in zip(offspring[::2], offspring[1::2]):
        #     if random.random() < self._params['cxpb']:
        #         self.toolbox.mate(child1, child2)
        #         del child1.fitness.values
        #         del child2.fitness.values

        for mutant in offspring:
            if random.random() < self._params['mutpb']:
                self.toolbox.mutate(mutant)
                del mutant.fitness.values

        #simple bound checking
        for i in range(len(offspring)):
            for j in range(len(offspring[i])):
                if offspring[i][j] > 1:
                    offspring[i][j] = 1
                if offspring[i][j] < -1:
                    offspring[i][j] = -1
        self._params['model_count'] += len([ind for ind in offspring if not ind.fitness.values])
        self.assign_fitnesses([ind for ind in offspring if not ind.fitness.valid])
        offspring.sort(reverse=True, key=lambda x: x.fitness)
        if len(self.halloffame) != 0:
            if offspring[0].fitness < self.halloffame[0].fitness:  # elitism- if none beat best so far it is reinserted
                offspring.insert(0, self.halloffame[0])
        self.population[:] = offspring[:self._params['popsize']]


class OptCMAES:
    """
    A covariance matrix adaptation evolutionary strategy optimization algorithm. Basically uses a covariance matrix
    at each step to identify the 'direction' of the optimal solution in the search space, and generates new individuals
    accordingly. Bound handling is achieved by moving any out of bounds parameters to the boundary condition. Other than
    that the implementation used here is as in the originating code from the deap module.
    """
    def __init__(self, **kwargs):
        super().__init__()
        self._params.update(**kwargs)
        self._params.setdefault('sigma', 0.3)
        self._params.setdefault('weights', 'superlinear')
        creator.create("Individual", list, fitness=creator.FitnessMin)

    def initialize_pop(self):
        """
        Generates the initial population and assigns fitnesses
        :return:
        """
        self.initialize_cma_es(sigma=self._params['sigma'], weights=self._params['weights'],
                               lambda_=self._params['popsize'], centroid=[0]*len(self._params['value_means']))
        self.toolbox.register("individual", self.make_individual)
        self.toolbox.register("generate", self.generate, self.toolbox.individual)
        self.toolbox.register("population", tools.initRepeat, list, self.initial_individual)
        self.toolbox.register("update", self.update)
        self.population = self.toolbox.population(n=self._params['popsize'])
        self.assign_fitnesses(self.population)
        self._params['model_count'] += len(self.population)

    def initial_individual(self):
        """
        Generates an individual with random parameters within bounds
        :return:
        """
        ind = creator.Individual([random.uniform(-1, 1) for _ in range(len(self._params['value_means']))])
        return ind

    def update_pop(self):
        self.toolbox.generate()
        #simple bound checking
        for i in range(len(self.population)):
            for j in range(len(self.population[i])):
                if self.population[i][j] > 1:
                    self.population[i][j] = 1
                if self.population[i][j] < -1:
                    self.population[i][j] = -1
        self.assign_fitnesses(self.population)
        self.toolbox.update(self.population)
        self._params['model_count'] += len(self.population)

    def make_individual(self, paramlist):
        part = creator.Individual(paramlist)
        part.ident = None
        return part

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
        self.centroid = numpy.array([0]*len(self._params['value_means']))

        self.dim = len(self.centroid)
        self.sigma = self.params.get("sigma", 0.5)
        self.pc = numpy.zeros(self.dim)
        self.ps = numpy.zeros(self.dim)
        self.chiN = numpy.sqrt(self.dim) * (1 - 1. / (4. * self.dim) + 1. / (21. * self.dim ** 2))

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
            self.weights = numpy.log(self.mu + 0.5) - numpy.log(numpy.arange(1, self.mu + 1))
        elif rweights == "linear":
            self.weights = self.mu + 0.5 - numpy.arange(1, self.mu + 1)
        elif rweights == "equal":
            self.weights = numpy.ones(self.mu)
        else:
            raise RuntimeError("Unknown weights : %s" % rweights)

        self.weights /= sum(self.weights)
        self.mueff = 1. / sum(self.weights ** 2)

        self.cc = params.get("ccum", 4. / (self.dim + 4.))
        self.cs = params.get("cs", (self.mueff + 2.) / (self.dim + self.mueff + 3.))
        self.ccov1 = params.get("ccov1", 2. / ((self.dim + 1.3) ** 2 + self.mueff))
        self.ccovmu = params.get("ccovmu", 2. * (self.mueff - 2. + 1. / self.mueff) /
                                 ((self.dim + 2.) ** 2 + self.mueff))
        self.ccovmu = min(1 - self.ccov1, self.ccovmu)
        self.damps = 1. + 2. * max(0, numpy.sqrt((self.mueff - 1.) / (self.dim + 1.)) - 1.) + self.cs
        self.damps = params.get("damps", self.damps)


class DE_Opt(OptDE, BaseScore):
    """
    Class for DE algorithm optimizing BUFF fitness
    """
    def __init__(self, topology, **kwargs):
        super().__init__(**kwargs)
        self._params['topology'] = topology


class DE_RMSD(OptDE, BaseRMSD):
    """
    Class for DE algorithm optimizing RMSD against target model
    """
    def __init__(self, topology, ref_pdb, **kwargs):
        super().__init__(**kwargs)
        self._params['topology'] = topology
        self._params['ref_pdb'] = ref_pdb


class DE_Comparator(OptDE, BaseComparator):
    """
    Class for DE algorithm optimizing BUFF fitness change on docking two AMPAL objects.
    """
    def __init__(self, top1, top2, params1, params2, seq1, seq2, **kwargs):
        super().__init__(**kwargs)
        self._params['top1'] = top1
        self._params['top2'] = top2
        self._params['params1'] = params1
        self._params['params2'] = params2
        self._params['seq1'] = seq1
        self._params['seq2'] = seq2
        obj1 = top1(*params1)
        obj1.pack_new_sequences(seq1)
        obj2 = top2(*params2)
        obj2.pack_new_sequences(seq2)
        self._params['ref1'] = obj1.buff_interaction_energy.total_energy
        self._params['ref2'] = obj2.buff_interaction_energy.total_energy


class PSO_Opt(OptPSO, BaseScore):
    """
    Class for PSO algorithm optimizing BUFF fitness
    """
    def __init__(self, topology, **kwargs):
        super().__init__(**kwargs)
        self._params['topology'] = topology


class PSO_RMSD(OptPSO, BaseRMSD):
    """
    Class for PSO algorithm optimizing RMSD against target model
    """
    def __init__(self, topology, ref_pdb, **kwargs):
        super().__init__(**kwargs)
        self._params['topology'] = topology
        self._params['ref_pdb'] = ref_pdb


class PSO_Comparator(OptPSO, BaseComparator):
    """
    Class for PSO algorithm optimizing BUFF fitness change on docking two AMPAL objects.
    """
    def __init__(self, top1, top2, params1, params2, seq1, seq2, **kwargs):
        super().__init__(**kwargs)
        self._params['top1'] = top1
        self._params['top2'] = top2
        self._params['params1'] = params1
        self._params['params2'] = params2
        self._params['seq1'] = seq1
        self._params['seq2'] = seq2
        obj1 = top1(*params1)
        obj1.pack_new_sequences(seq1)
        obj2 = top2(*params2)
        obj2.pack_new_sequences(seq2)
        self._params['ref1'] = obj1.buff_interaction_energy.total_energy
        self._params['ref2'] = obj2.buff_interaction_energy.total_energy


class GA_Opt(OptGA, BaseScore):
    """
    Class for GA algorithm optimizing BUFF fitness
    """
    def __init__(self, topology, **kwargs):
        super().__init__(**kwargs)
        self._params['topology'] = topology


class GA_RMSD(OptGA, BaseRMSD):
    """
    Class for GA algorithm optimizing RMSD against target model
    """
    def __init__(self, topology, ref_pdb, **kwargs):
        super().__init__(**kwargs)
        self._params['topology'] = topology
        self._params['ref_pdb'] = ref_pdb


class GA_Comparator(OptGA, BaseComparator):
    """
    Class for GA algorithm optimizing BUFF fitness change on docking two AMPAL objects.
    """
    def __init__(self, top1, top2, params1, params2, seq1, seq2, **kwargs):
        super().__init__(**kwargs)
        self._params['top1'] = top1
        self._params['top2'] = top2
        self._params['params1'] = params1
        self._params['params2'] = params2
        self._params['seq1'] = seq1
        self._params['seq2'] = seq2
        obj1 = top1(*params1)
        obj1.pack_new_sequences(seq1)
        obj2 = top2(*params2)
        obj2.pack_new_sequences(seq2)
        self._params['ref1'] = obj1.buff_interaction_energy.total_energy
        self._params['ref2'] = obj2.buff_interaction_energy.total_energy


class CMAES_Opt(OptCMAES, BaseScore):
    """
    Class for CMAES algorithm optimizing BUFF fitness
    """
    def __init__(self, topology, **kwargs):
        super().__init__(**kwargs)
        self._params['topology'] = topology


class CMAES_RMSD(OptCMAES, BaseRMSD):
    """
    Class for CMAES algorithm optimizing RMSD against target model
    """
    def __init__(self, topology, ref_pdb, **kwargs):
        super().__init__(**kwargs)
        self._params['topology'] = topology
        self._params['ref_pdb'] = ref_pdb


class CMAES_Comparator(OptCMAES, BaseComparator):
    """
    Class for CMAES algorithm optimizing BUFF fitness change on docking two AMPAL objects.
    """
    def __init__(self, top1, top2, params1, params2, seq1, seq2, **kwargs):
        super().__init__(**kwargs)
        self._params['top1'] = top1
        self._params['top2'] = top2
        self._params['params1'] = params1
        self._params['params2'] = params2
        self._params['seq1'] = seq1
        self._params['seq2'] = seq2
        obj1 = top1(*params1)
        obj1.pack_new_sequences(seq1)
        obj2 = top2(*params2)
        obj2.pack_new_sequences(seq2)
        self._params['ref1'] = obj1.buff_interaction_energy.total_energy
        self._params['ref2'] = obj2.buff_interaction_energy.total_energy


def buff_eval(params):
    """Builds and evaluates BUFF energy of model in parallelization

    Parameters
    ----------
    params: list
        Tuple containing the topology to be built, the sequence, and the parameters for model building.

    Returns
    -------
    model.bude_score: float
        BUFF score for model to be assigned to particle fitness value.
    """
    topology, sequence, parsed_ind = params
    model = topology(*parsed_ind)
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
    topology, sequence, parsed_ind, reference_pdb = rmsd_params
    model = topology(*parsed_ind)
    model.pack_new_sequences(sequence)
    ca, bb, aa = run_profit(model.pdb, reference_pdb, path1=False, path2=False)
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
