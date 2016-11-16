#ISAMBARD
###Intelligent System for Analysis, Model Building And Rational Design.
#### Version 2016.3 (Nov 11, 2016), Woolfson Group, University of Bristol.
[![Documentation Status](https://readthedocs.org/projects/isambard/badge/?version=latest)](http://isambard.readthedocs.io/en/latest/?badge=latest)
[![CircleCI](https://circleci.com/gh/woolfson-group/isambard.svg?style=shield&circle-token=27387ac82a6d30c7bd6a72ce3214fa57677e9d87)](https://circleci.com/gh/woolfson-group/isambard)
[![Gitter](https://img.shields.io/gitter/room/nwjs/nw.js.svg?maxAge=2592000)](https://gitter.im/woolfson-group/isambard?utm_source=share-link&utm_medium=link&utm_campaign=share-link)
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/woolfson-group/isambard/blob/master/LICENSE.md)

ISAMBARD (Intelligent System for Analysis, Model Building And Rational Design), is a Python-based framework for 
structural analysis and rational design of biomolecules. It is developed and maintained by members of the
[Woolfson group, University of Bristol](http://www.chm.bris.ac.uk/org/woolfson/index.html).

## Quick Start

### Basic Install

ISAMBARD can be installed straight from PYPI using pip:

```
pip install isambard
```

ISAMBARD has a few Python package requirements, just install them through pip if it asks for them. We recommend using the [Anaconda Python 3 distribution](https://www.continuum.io/downloads), it contains most of the dependencies. 

**Windows People** - You'll need to download [Visual Studio with the Visual C++ tools](https://www.visualstudio.com/vs/cplusplus/), if you want to use ISAMBARD. This is because we use a package called Cython to make the code run fast, and it needs a C compiler. *You can also use MinGW, but it's more complicated to [set up with Cython](http://cython.readthedocs.io/en/latest/src/tutorial/appendix.html).*

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

##Principal Investigator
* Derek N. Woolfson (d.n.woolfson@bristol.ac.uk)

##Developers
###Core Dev Team
####Woolfson Group
* Gail J. Bartlett (g.bartlett@bristol.ac.uk)
* Jack W. Heal (jack.heal@bristol.ac.uk)
* Kieran L. Hudson (kieran.hudson@bristol.ac.uk)
* Andrew R. Thomson (drew.thomson@bristol.ac.uk)
* Christopher W. Wood (chris.wood@bristol.ac.uk)

###Contributors
####Woolfson Group
* Caitlin Edgell
* Kathryn L. Porter Goff

###BUDE Dev Team
####Sessions Group
* Amaurys À. Ibarra (amaurys.avilaibarra@bristol.ac.uk)
* Richard B. Sessions (r.sessions@bristol.ac.uk)
