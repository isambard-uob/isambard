Anyone contributions to ISAMBARD would be greatly appreciated, so please download the repo, make some changes and send us a pull request! If you do submit any new functionality, please try and document it well and add some unit tests to the testing suite. If you are working on something, please let us know about it (if you can), so that we can coordinate our development efforts. Finally, when submitting code, please try and adhere to style guide below.

# ISAMBARD Coding Style Guide

The style of ISAMBARD code should *almost* always adhere to [PEP8](https://www.python.org/dev/peps/pep-0008/). The few exceptions are as follows:

1. Lines of code should be no longer than 120 characters.
1. Comments should be no longer than 80 characters including the newline character, so they're effectively 79 chars.

## General Style Guidelines

### Doc Strings

Please use the [NumPy docstring style](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt) to document your code.

### Classes
1. Do not write getters and setters manually, please use the [property decorator and corresponding setter decorator](https://docs.python.org/2/library/functions.html#property).
1. Use [class attributes](https://docs.python.org/2/tutorial/classes.html) (section 9.3.5) to keep the `__init__` method lightweight, by only defining an attribute in the `__init__` method if it is **unique to the instance**.
1. Use [class methods](https://docs.python.org/2/library/functions.html#classmethod) to define alternate ways to instantiate a class. Try to avoid having multiple ways to instantiate a class in the `__init__` method as this can be difficult to read. Consider which should be the preferred method of instantiation and use that to define the `__init__`.
1. Class methods should have names that start with 'from'. For example `Quaternion.from_axis_and_angle`.

## ISAMBARD Specific Style Guidelines

### AMPAL Objects
1. Any class that inherits from `BaseAmpal`, or any of its subclasses, is required to have following methods:
    1. `get_atoms`
    1. `make_pdb`

##  External Program Interfaces
1. All strings produced by and fed into ISAMBARD *must* be in unicode, therefore if the external program requires byte strings, any decoding and decoding must occur in the body of the external program handler function.
1. Any external program handler function that interacts with structural data from ISAMBARD should take a PDB string as an input (or a path to a pdb file) and if it returns any structural data it should return a PDB string, **not** an AMPAL object.
