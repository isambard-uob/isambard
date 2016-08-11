from deap import algorithms, base, creator, tools
import random
import numpy
import matplotlib.pylab as plt
import datetime
import operator
from concurrent import futures

# TODO: Unit tests, tidy docstrings, compatibility with new polymer object


class OptPSO:
    """Parent class for running particle swarm optimization as implemented in DEAP. For most practical purposes it will
    be better to use the parallelized 'px' version.

    Notes
    -----
    This is a tweaked version of the simple example given in the DEAP documentation. Have updated the update particle
    algorithm to include an inertial term. The implementation here is a little different from the GA in that each
    particle is restricted to range -1 to +1, and is scaled afterwards to give the parameters for model building. This
    works out a little easier for speed ranges and bounds handling etc.

    Parameters
    ----------
    specification: specification class
        specification to be minimized- needs an evaluate method.
    output path: str
        path to save log and optimized model
    """
    def __init__(self, specification, output_path, bude_mode='average'):
        self.specification = specification
        self.output_path = output_path
        self.bude_mode = bude_mode
        self.toolbox = base.Toolbox()
        self.halloffame = tools.HallOfFame(maxsize=1)
        self.plot = None
        self.log = None
        self.value_ranges = None
        self.value_means = None
        self.minvalues = None
        self.maxvalues = None
        self.sequence = None
        self.smin = None
        self.smax = None
        self.pmin = None
        self.pmax = None
        self.arrangement = None
        self.keys = None
        self.size = None
        self.swarm = None
        self.stats = None
        self.gen = None
        self.numgen = None
        self.start_time = None
        self.end_time = None
        self.run_id = None
        self.logbook = None
        self.neighbours = None
        self.modelcount = None
        self.procs = None
        self.plot = None

        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Particle", list, fitness=creator.FitnessMin, speed=list, smin=None, smax=None, best=None)

    def parameters(self, sequence, value_means, value_ranges, arrangement):
        """Relates the individual to be evolved to the full parameter string for building the specification object
        Notes
        -----
        Note that the internal working of the PSO is slightly different in that the particles in the algorithm are all
        scaled to between -1 and +1, and then expanded/contracted to the appropriate range for model building. This
        makes no difference to the user, but changes how some things are handled internally.

        Parameters
        ----------
        sequence: str
            Full amino acid sequence for specification object to be optimized
        value_means: list
            List containing mean values for parameters to be optimized
        value_ranges: list
            List containing ranges for parameters to be optimized. Values must be positive.
        arrangement: list
            Full list of fixed and variable parameters for model building. Fixed values are the appropriate value.
            Values to be varied should be listed as 'var0', 'var1' etc, and must be in ascending numerical order.
            Variables can be repeated if required.

        Raises
        ------
        ValueError if parameter lists are of mismatching lengths
        """
        self.value_ranges = value_ranges
        if any(x <= 0 for x in self.value_ranges):
            raise ValueError("range values must be greater than zero")
        self.value_means = value_means
        if len(self.value_means) != len(self.value_ranges):
            raise ValueError("value_means and value_ranges must be of same length")
        self.sequence = sequence
        self.minvalues = [a - b for a, b in zip(self.value_means, self.value_ranges)]
        self.maxvalues = [a + b for a, b in zip(self.value_means, self.value_ranges)]
        self.smax = None
        self.smin = None
        self.pmin = -1
        self.pmax = 1
        self.arrangement = arrangement
        self.keys = []
        self.size = len(self.minvalues)
        for i in range(len(self.value_means)):
            self.keys.append("".join(['var', str(i)]))
        if len(set(arrangement).intersection(self.keys)) != len(self.value_means):
            raise ValueError("argument mismatch!")

        self.toolbox.register("particle", self.generate)
        self.toolbox.register("update", self.update_particle_chi_neighbours)
        creator.create("Swarm", list, gbest=None, gbestfit=creator.FitnessMin)
        self.toolbox.register("swarm", tools.initRepeat, creator.Swarm, self.toolbox.particle)

    def generate(self):
        """Generates a particle using the creator function. Position and speed are uniformly randomly seeded within
        allowed bounds. The particle also has speed limit settings taken from global values.

        Returns
        -------
        particle object
        """
        part = creator.Particle([random.uniform(self.pmin, self.pmax) for _ in range(self.size)])
        part.speed = [random.uniform(self.smin, self.smax) for _ in range(self.size)]
        part.smin = self.smin
        part.smax = self.smax
        part.index = None
        part.neighbours = None
        return part

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
        for i in range(len(self.value_means)):
            scaled_ind.append(self.value_means[i] + (individual[i] * self.value_ranges[i]))
        fullpars = list(self.arrangement)
        for k in range(len(self.keys)):
            for j in range(len(fullpars)):
                if fullpars[j] == self.keys[k]:
                    fullpars[j] = scaled_ind[k]
        return fullpars

    def update_particle_chi_neighbours(self, part, chi=0.729843788, c=2.05):
        """Constriction factor update particle method that looks for a list of neighbours attached to a particle and
        uses the particle's best position and that of the best neighbour.
        """
        neighbour_pool = [self.swarm[i] for i in part.neighbours]+[part]
        best_neighbour = min(neighbour_pool, key=lambda x: x.best.fitness.values[0])
        ce1 = (c * random.uniform(0, 1) for _ in range(len(part)))
        ce2 = (c * random.uniform(0, 1) for _ in range(len(part)))
        ce1_p = map(operator.mul, ce1, map(operator.sub, part.best, part))
        ce2_g = map(operator.mul, ce2, map(operator.sub, best_neighbour.best, part))
        # ce2_g = map(operator.mul, ce2, map(operator.sub, self.swarm.gbest, part))
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

    def run_pso(self, swarmsize, numgen, procs, neighbours, max_speed=0.75, run_id='PSOmodel', plot=False, log=False):
        """Method to start PSO running.

        Notes
        -----
        Each particle is connected to the others in a cyclic specification, and is allowed to read a best value only from
        neighbouring particles. This specification helps to prevent the swarm from becoming trapped in local minima.

        Parameters
        ----------
        swarmsize: int
            The size of the swarm. Size depends on how costly models are. Values of 10-50 are common.
        numgen: int
            Number of generations to run the PSO for. Values of 50-100 are common
        procs: int
            The number of processors to run the job on
        neighbours: int
            The number of neighbouring particles visible to each particle in either direction. I.e. a value of 4 will
            allow a particle to recognise eight total neighbours.
        max_speed: float
            Maximum speed that particles are allowed to reach. Typical value of 0.1.
        run_id: str
            String giving file identifier for best model and logfile.

        Modifies
        --------
        swarm: DEAP object
            Updates the swarm iteratively according to the update function chosen.
        logbook: DEAP object
            Records the stats for each generation
        """
        random.seed()
        self.plot = plot
        self.log = log
        self.toolbox.unregister("update")
        self.toolbox.register("update", self.update_particle_chi_neighbours)
        self.neighbours = neighbours
        if self.neighbours > swarmsize/2:
            raise ValueError("Number of neighbours in each direction should not be greater than the size of the swarm")
        self.smin = -max_speed
        self.smax = max_speed
        self.gen = 0
        self.run_id = run_id
        self.start_time = datetime.datetime.now()
        print ('Starting minimisation ({:%Y-%m-%d %H:%M:%S})'.format(self.start_time))

        self.swarm = self.toolbox.swarm(n=swarmsize)
        for i in range(len(self.swarm)):
            self.swarm[i].index = i
            self.swarm[i].neighbours = list(set([(i-x) % len(self.swarm) for x in range(1, self.neighbours+1)] +
                                                [(i+x) % len(self.swarm) for x in range(1, self.neighbours+1)]))
        # assigns self.neighbours number of (unique) neighbours to either side of particle i

        self.stats = tools.Statistics(lambda ind: ind.fitness.values)
        self.stats.register("avg", numpy.mean)
        self.stats.register("std", numpy.std)
        self.stats.register("min", numpy.min)
        self.stats.register("max", numpy.max)
        self.logbook = tools.Logbook()
        self.logbook.header = ["gen", "evals"] + self.stats.fields
        self.numgen = numgen
        self.modelcount = 0

        for g in range(self.numgen):
            self.gen = g
            # split out into groups to eval and not eval
            valid_particles = []
            invalid_particles = []
            for part in self.swarm:
                if any(x > 1 or x < -1 for x in part):
                    invalid_particles.append(part)
                else:
                    valid_particles.append(part)
            self.modelcount += len(valid_particles)
            topologies = [self.specification]*len(valid_particles)
            sequences = [self.sequence]*len(valid_particles)
            parameters = [self.parse_individual(part) for part in valid_particles]
            mode = [self.bude_mode] * len(valid_particles)
            px_parameters = zip(topologies, sequences, parameters, mode)
            with futures.ProcessPoolExecutor(max_workers=procs) as executor:
                fitnesses = executor.map(px_eval, px_parameters)
            for part, fit in zip(valid_particles, fitnesses):
                part.fitness.values = (fit,)
            self.swarm[:] = valid_particles + invalid_particles
            self.swarm.sort(key=lambda x: x.index)
            for part in self.swarm:
                if not part.best or part.best.fitness.values[0] > part.fitness.values[0]:
                    part.best = creator.Particle(part)
                    part.best.fitness.values = part.fitness.values
                if not self.swarm.gbest or self.swarm.gbestfit[0] > part.fitness.values[0]:
                    self.swarm.gbest = creator.Particle(part)
                    self.swarm.gbestfit = part.fitness.values
            self.halloffame.update(self.swarm)
            self.logbook.record(gen=g, evals=len(valid_particles), **self.stats.compile(self.swarm))
            print(self.logbook.stream)
            for part in self.swarm:
                self.toolbox.update(part)
        self.end_time = datetime.datetime.now()
        time_taken = self.end_time - self.start_time
        print("End of minimisation ({:%Y-%m-%d %H:%M:%S})".format(self.end_time))
        print("Run ID is {0}".format(self.run_id))
        print('Minimization time = {0}'.format(time_taken))
        print("Evaluated {0} models in total".format(self.modelcount))
        print("Best score is {0}".format(self.swarm.gbestfit[0]))
        print("Best parameters are {0}".format(self.parse_individual(self.swarm.gbest)))
        for i, entry in enumerate(self.swarm.gbest):
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
