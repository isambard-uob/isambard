#ISAMBARD
###Intelligent System for Analysis, Model Building And Rational Design.
#### Version 2016.4.1 (Nov 29, 2016), Woolfson Group, University of Bristol.
[![Documentation Status](https://readthedocs.org/projects/isambard/badge/?version=latest)](http://isambard.readthedocs.io/en/latest/?badge=latest)
[![CircleCI](https://circleci.com/gh/woolfson-group/isambard.svg?style=shield&circle-token=27387ac82a6d30c7bd6a72ce3214fa57677e9d87)](https://circleci.com/gh/woolfson-group/isambard)
[![codecov](https://codecov.io/gh/woolfson-group/isambard/branch/master/graph/badge.svg)](https://codecov.io/gh/woolfson-group/isambard)
[![Gitter](https://img.shields.io/gitter/room/nwjs/nw.js.svg?maxAge=2592000)](https://gitter.im/woolfson-group/isambard?utm_source=share-link&utm_medium=link&utm_campaign=share-link)
[![PyPI version](https://badge.fury.io/py/isambard.svg)](https://badge.fury.io/py/isambard)
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/woolfson-group/isambard/blob/master/LICENSE.md)
![Python Version](https://img.shields.io/badge/python-3.5%2C%203.6-lightgrey.svg)

ISAMBARD (Intelligent System for Analysis, Model Building And Rational Design), is a Python-based framework for 
structural analysis and rational design of biomolecules. It is developed and maintained by members of the
[Woolfson group, University of Bristol](http://www.chm.bris.ac.uk/org/woolfson/index.html).

## Quick Start

### Basic Install

ISAMBARD can be installed straight from PYPI using pip:

```
pip install isambard
```
Or if you want to try an experimental build, download from GitHub either by downloading the zipped file or cloning, then navigate to the ISAMBARD folder and type:

```
pip install .
```

ISAMBARD has a few Python package requirements, just install them through pip if it asks for them. We recommend using the [Anaconda Python 3 distribution](https://www.continuum.io/downloads), it contains most of the dependencies. 

**Windows Users** - You'll need to download a C compiler, if you want to use ISAMBARD. This is because we use a package called Cython to make the code run fast. 
The exact version of Microsoft Visual C++ compiler depends on the python version that you are using. Here are the [installation instructions](https://wiki.python.org/moin/WindowsCompilers).  

*You can also use MinGW, but it's not currently supported. [Instructions here](http://cython.readthedocs.io/en/latest/src/tutorial/appendix.html).*

### External Programs

To get the most out of ISAMBARD, a couple of external programs are recommended:

1. [Scwrl4](http://dunbrack.fccc.edu/scwrl4/) - Used to pack sidechains, it's fast and accurate. Free for non-commercial use.
1. [DSSP](http://swift.cmbi.ru.nl/gv/dssp/) - Used to find secondary structure in models. Free to download.

When you first import ISAMBARD you'll be asked to give paths to the executables for the dependencies. These are not required to use ISAMBARD, if you don't want to use any of them, just leave the path blank. Once you've finished setting your paths and some other options, a small file called `.isambard_settings` will be made in your home directory, which contains your settings. If you ever want to rerun the configure script that runs when you first import, you can run the following command:

```python
import isambard
isambard.settings.configure()
```

Chat with us on Gitter if you get stuck (link above), or raise an issue!

### Once ISAMBARD is installed...

You might want to take a look at the [wiki](https://github.com/woolfson-group/isambard/wiki), there are a range of tutorials that demonstrate various aspects of ISAMBARD. These are IPython notebooks, so please download them and run through them (modify/hack/break them) on your own machine.

Wanting to delve a bit deeper? Take a look at the [docs](http://isambard.readthedocs.io/en/latest/) to find out more of the features in ISAMBARD, or just take a look through the code base and hack around. Remember, feel free to contact us on Gitter, email or through the issues if you get stuck.

### Unit tests 
To ensure that the package is correctly installed unit tests can be run by:

`python -m unittest discover unit_tests/`

##Principal Investigator
* Derek N. Woolfson
  * Email: d.n.woolfson@bristol.ac.uk

##Developers
###Core Dev Team
####Woolfson Group
* Gail J. Bartlett 
  * Email: g.bartlett@bristol.ac.uk
* Jack W. Heal
  * Email: jack.heal@bristol.ac.uk
  * Twitter: @JackHeal
* Kieran L. Hudson
  * Email: kieran.hudson@bristol.ac.uk
* Andrew R. Thomson
  * Email: drew.thomson@bristol.ac.uk
* Christopher W. Wood
  * Email: chris.wood@bristol.ac.uk
  * Twitter: @ChrisWellsWood

###Contributors
####Woolfson Group
* Caitlin Edgell
* Kathryn L. Porter Goff

###BUDE Dev Team
####Sessions Group
* Amaurys À. Ibarra
  * Email: amaurys.avilaibarra@bristol.ac.uk
* Richard B. Sessions
  * Email: r.sessions@bristol.ac.uk
