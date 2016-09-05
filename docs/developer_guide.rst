###############
Developer Guide
###############


General style guide
###################

The style of ISAMBARD code should *almost* always adhere to `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_.
The few exceptions are as follows:

- Lines of code should be no longer than 120 characters.
- Comments should be no longer than 80 characters including the newline character, so they're effectively 79 characters.


**Classes**

- Do not write getters and setters manually, please use the
  `property decorator and corresponding setter decorator <https://docs.python.org/2/library/functions.html#property>`_.
- Use `class attributes <https://docs.python.org/2/tutorial/classes.html>`_ (section 9.3.5)
  to keep the ``__init__`` method lightweight, by only defining an attribute in the __init__ method if it is
  **unique to the instance**.
- Use `class methods <https://docs.python.org/2/library/functions.html#classmethod>`_ to define alternate ways to
  instantiate a class. Try to avoid having multiple ways to instantiate a class in the ``__init__`` method as this can
  be difficult to read. Consider which should be the preferred method of instantiation and use that to define the
  ``__init__``.
- Class methods should have names that start with ``from``, for example: ``Quaternion.from_axis_and_angle``,
  ``CoiledCoil.from_parameters``.

ISAMBARD-specific style guide
#############################

**AMPAL objects**

- Any class that inherits from ``BaseAmpal``, or any of its subclasses, is required to have following methods:
    - ``get_atoms``
    - ``make_pdb``

**External program interfaces**

- All strings produced by and fed into ISAMBARD *must* be in unicode,therefore if the external program requires
  byte strings, any decoding and decoding must occur in the body of the external program handler function.
- Any external program handler function that interacts with structural data from ISAMBARD should take a PDB string
  as an input (or a path to a pdb file) and if it returns any structural data it should return a PDB string, *not* an
  AMPAL object.