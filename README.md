# ISAMBARD

Intelligent System for Analysis, Model Building And Rational Design.

[![CircleCI](https://circleci.com/gh/isambard-uob/isambard.svg?style=shield)](https://circleci.com/gh/isambard-uob/isambard)
[![Python Version](https://img.shields.io/badge/python-3.5%2C%203.6-lightgrey.svg)](https://woolfson-group.github.io/isambard/)
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/woolfson-group/isambard/blob/master/LICENSE.md)

ISAMBARD is a Python-based framework for structural analysis and rational
design of biomolecules, with a particular focus on parametric modelling of
proteins. It is developed and maintained by members of the [Woolfson group, University of Bristol](http://www.chm.bris.ac.uk/org/woolfson/index.html).

### Citing ISAMBARD
Any publication arising from use of the ISAMBARD software package should cite the following reference:

[Wood CW *et al* (2017) ISAMBARD: an open-source computational environment for biomolecular analysis, modelling and design. *Bioinformatics*, **33**, 3043-50](https://doi.org/10.1093/bioinformatics/btx352)

## Installation

ISAMBARD can be installed straight from PyPI using `pip`:

```
pip install isambard
```
Or if you want to try an experimental build (you'll need a C compiler), download
from GitHub either by downloading the zipped file or cloning, then navigate to
the ISAMBARD folder and type:

```
pip install .
```

## External Programs

If you want to add side chains to your designs, you need to have [Scwrl4](
http://dunbrack.fccc.edu/scwrl4/) installed and available on your system path. 

## Quick Start

> Note<br />
> If you're not sure what parametric modelling of proteins is, have a
> play with [CCBuilder 2.0](http://coiledcoils.chm.bris.ac.uk/ccbuilder2/builder).

Let's build a coiled-coil dimer with typical parameters:

```Python
import isambard.specifications as specifications
import isambard.modelling as modelling
import isambard.optimisation

my_dimer = specifications.CoiledCoil.from_parameters(2, 28, 5, 225, 283)
dimer_sequences = [
    'EIAALKQEIAALKKENAALKWEIAALKQ',
    'EIAALKQEIAALKKENAALKWEIAALKQ'
]
my_dimer = modelling.pack_side_chains_scwrl(my_dimer, dimer_sequences)
print(my_dimer.pdb)
# OUT: 
# HEADER ISAMBARD Model                                                                  
# ATOM      1  N   GLU A   1      -5.364  -1.566  -0.689  1.00  0.00           N  
# ATOM      2  CA  GLU A   1      -4.483  -2.220   0.308  1.00  0.00           C  
# ATOM      3  C   GLU A   1      -3.886  -1.143   1.216  1.00  0.00           C  
# ATOM      4  O   GLU A   1      -3.740  -1.337   2.425  1.00  0.00           O  
# ATOM      5  CB  GLU A   1      -3.389  -3.028  -0.392  1.00  0.00           C  
# ...
```

Don't know what your parameters might be? Let's optimise them then!

```Python
import budeff
import isambard.optimisation.evo_optimizers as ev_opts
from isambard.optimisation.evo_optimizers import Parameter

specification = specifications.CoiledCoil.from_parameters
sequences = [
    'EIAALKQEIAALKKENAALKWEIAALKQ',
    'EIAALKQEIAALKKENAALKWEIAALKQ'
]
parameters = [
    Parameter.static('Oligomeric State', 2),
    Parameter.static('Helix Length', 28),
    Parameter.dynamic('Radius', 5.0, 1.0),
    Parameter.dynamic('Pitch', 200, 60),
    Parameter.dynamic('PhiCA', 283, 27),  # 283 is equivalent a g position
]
def get_buff_total_energy(ampal_object):
    return budeff.get_internal_energy(ampal_object).total_energy
opt_ga = ev_opts.GA(specification, sequences, parameters, get_buff_total_energy)
opt_ga.run_opt(100, 5, cores=8)
# OUT:
# gen	evals	avg     	std    	min     	max     
# 0  	61   	-820.401	42.0119	-908.875	-750.001
# 1  	59   	-859.86 	31.4194	-950.15 	-807.265
# 2  	60   	-887.028	23.8683	-951.153	-847.346
# 3  	70   	-907.257	15.9615	-952.863	-882.028
# 4  	81   	-922.522	14.6206	-972.335	-903.444
# Evaluated 431 models in total in 0:00:29.523487
# Best fitness is (-972.3348571854714,)
# Best parameters are [2, 28, 4.678360526981807, 151.35365923229745, 277.2061538048508]
optimized_model = opt_ga.best_model
```

This quick example of parametric modelling with ISAMBARD, the next thing to do
is take a look at the [docs](https://isambard-uob.github.io/isambard/) from
tutorials on the tools available, or just take a look through the code base and
hack around. Feel free to contact us through email or the issues if you get
stuck.
