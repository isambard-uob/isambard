###############
Getting Started
###############

Installation
============

ISAMBARD can be installed straight from PYPI using pip:

.. code-block::

  pip install isambard

ISAMBARD has a few Python package requirements, just install them through pip if it asks for them. We recommend using the [Anaconda Python 3 distribution](https://www.continuum.io/downloads), it contains most of the dependencies. 

**Windows People** - You'll need to download `Visual Studio with the Visual C++ tools <https://www.visualstudio.com/vs/cplusplus/>`_, if you want to use ISAMBARD. This is because we use a package called Cython to make the code run fast, and it needs a C compiler. You can also use MinGW, but it's more complicated to `set up with Cython <http://cython.readthedocs.io/en/latest/src/tutorial/appendix.html>`_.

External Programs
=================

To get the most out of ISAMBARD, a couple of external programs are recommended:

1. [Scwrl4](http://dunbrack.fccc.edu/scwrl4/) - Used to pack sidechains, it's fast and accurate. Free for non-commercial use.
1. [DSSP](http://swift.cmbi.ru.nl/gv/dssp/) - Used to find secondary structure in models. Free to download.

When you first import ISAMBARD you'll be asked to give paths to the executables for the dependencies. These are not required to use ISAMBARD, if you don't want to use any of them, just leave the path blank. Once you've finished setting your paths and some other options, a wee file called `.isambard_settings` will be made in your home directory, which contains your settings. If you ever want to rerun the configure script that runs when you first import, you can run the following command:

.. code-block:: python

  import isambard
  isambard.settings.configure()

Chat with us on Gitter if you get stuck (link above), or raise an issue!

Other External Programs
-----------------------

ISAMBARD provides a convenient interface for running the following external programs:

- `Naccess <http://www.bioinf.manchester.ac.uk/naccess/>`_
- `ProFit <http://www.bioinf.org.uk/programs/profit/index.html>`_
- `Reduce <http://kinemage.biochem.duke.edu/software/reduce.php>`_

You will need to have working copies of these programs in order to enable some parts of ISAMBARD's functionality.
