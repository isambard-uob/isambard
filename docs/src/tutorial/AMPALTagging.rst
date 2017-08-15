
Tagging in the AMPAL framework
==============================

1. Introduction
---------------

Tagging allows us to store information related to AMPAL objects and
retrieve it easily later.

All AMPAL objects have a ``tags`` attribute, which is a dictionary: a
'useful place to put stuff'. When we add items to the ``tags``
dictionary, we refer to it as 'tagging' the AMPAL object. Any Python
object can be stored in ``tags``.

The ``tags`` attribute is used to store information internally so that
certain values only need to be calculated once.

Properties you may have already encountered, such as ``helices`` and
``strands``, write to the ``tags`` attribute.

When you use ``convert_pdb_to_ampal``, you tag the ``Atoms`` of your
``Protein`` with information from the PDB file.

``BaseAmpal`` has some in-built 'tagging functions' that calculate
things and add the results to the ``tags`` of the appropriate AMPAL
object.

2. Getting started
------------------

First, let's import some functions from ``isambard`` for making an
``Assembly`` object.

.. code:: python

    import isambard
    from isambard.ampal import convert_pdb_to_ampal

We'll use the structure
`2ebo <http://www.ebi.ac.uk/pdbe/entry/pdb/2ebo>`__ for this. It's a
viral protein containing coiled-coil packing that is involved in fusing
the membranes of the ebola virus and its host.

.. code:: python

    a = convert_pdb_to_ampal('2ebo.pdb')

Our new ``Protein`` object ``a`` contains 3 Chains, and its ``id``
attribute is the pdb code we used to create it.

.. code:: python

    a




.. parsed-literal::

    <Assembly (2ebo) containing 3 Polypeptides, 215 Ligands>



.. code:: python

    a.id




.. parsed-literal::

    '2ebo'



It comes with a ``tags`` attribute, which is an empty dictionary.

.. code:: python

    a.tags




.. parsed-literal::

    {}



The 3 ``Polypeptides`` in our ``Assembly`` are also AMPAL objects ...

.. code:: python

    a[0] # The first Polypeptide




.. parsed-literal::

    <Polypeptide containing 74 Residues. Sequence: GLRQLANETTQA...>



... and therefore have their own ``tags``:

.. code:: python

    a[0].tags




.. parsed-literal::

    {}



Similarly, all of the ``Residues`` in the ``Polypeptide`` are AMPAL
objects with empty ``tags``.

.. code:: python

    a[0][0] # The first Residue of the first chain of the Protein




.. parsed-literal::

    <Residue containing 4 Atoms. Residue code: GLY>



.. code:: python

    a[0][0].tags




.. parsed-literal::

    {}



Each ``Residue`` in our ``Assembly`` is made from ``Atoms``. Right down
to the ``Atoms``, AMPAL objects have ``tags``.

Let's look at the ``tags`` of the backbone Nitrogen ``Atom`` of the
first ``Residue`` of the first ``Polypeptide`` of our ``Assembly``.

.. code:: python

    a[0][0]['N']




.. parsed-literal::

    <Nitrogen Atom (N). Coordinates: (-14.780, 25.698, -6.988)>



.. code:: python

    a[0][0]['N'].tags




.. parsed-literal::

    {'bfactor': 71.51, 'charge': '', 'occupancy': 1.0, 'state': 'A'}



| Not an empty dictionary!
| That's because we've used a pdb file to generate our ``Protein``
  object.
| There was extra information in there at the ``Atom`` level that we did
  not want to throw away, so we've added it to ``tags`` automatically.

| **In summary: **
| AMPAL objects come with ``tags`` dictionaries for us to use if we want
  to.
| Many of these are empty on instantiation so that the AMPAL objects
  don't carry unnecessary baggage and are lightweight to begin with.

If you're familiar with adding items to Python dictionaries, then you're
already familiar with adding items to ``tags``.

.. code:: python

    a.tags['protein_description'] = 'viral membrane fusion protein'

.. code:: python

    a.tags['number_of_atoms'] = len(list(a.get_atoms()))

Our ``Protein`` object now carries some information around with it in
its ``tags``.

.. code:: python

    a.tags




.. parsed-literal::

    {'number_of_atoms': 2021,
     'protein_description': 'viral membrane fusion protein'}



Exercises
~~~~~~~~~

1. Look at the page for 2ebo on the
   `PDB <http://www.rcsb.org/pdb/explore/explore.do?structureId=2ebo>`__
   or `PDBE <http://www.ebi.ac.uk/pdbe/entry/pdb/2ebo>`__. Find the
   resolution of the structure and add that to ``tags``.
2. Add another tag to the ``Assembly`` object storing its
   '``number_of_residues``'.
   HINT: Use a similar expression to one used for ``number_of_atoms``.
3. Tag each ``Polypeptide`` in the ``Assembly`` with a
   '``number_of_atoms``' tag, like we did with the whole ``Assembly``
   earlier.
4. Tag each ``Residue`` in the first ``Polypeptide`` with a
   '``number_of_atoms``' tag.
5. Tag the second ``Polypeptide`` with a 'number\_of\_tryptophans' tag.
   HINT: You may want to use ``Counter`` for this. Use
   "``from collections import Counter``" to bring this into your
   namespace.

After completing Exercise 1. and 2., you should have 4 items in your
``a.tags`` dictionary.

.. code:: python

    a.tags




.. parsed-literal::

    {'number_of_atoms': 2021,
     'protein_description': 'viral membrane fusion protein'}



3. Tagging torsion angles
-------------------------

In an earlier tutorial we used the ``isambard.geometry.dihedral``
function to calculate the backbone torsion angles of a ``Residue`` in
the AMPAL framework.

| Suppose we wanted to calculate all of the torsion angles in an
  ``Assembly``.
| We could loop over all of the ``Residues`` in the ``Assembly`` using
  the ``get_monomers()`` method, call the ``isambard.geometry.dihedral``
  function at each stage of the loop and store the results somewhere
  convenient.

| This would be perfectly valid.
| However, we also have an in-built method for doing this.

The ``tag_torsion_angles`` method calculates the torsion angles for each
``Residue`` in an ``Assembly`` (or ``Polypeptide``) and adds them
``tags`` dictionary of the ``Residue``.

.. code:: python

    a.tag_torsion_angles()

The ``.tags`` dictionary for the first Residue now contains values for
its ``omega``, ``phi`` and ``psi`` angles.

Since it's the first ``Residue`` of the ``Polypeptide``, its ``omega``
and ``phi`` torsion angles are not defined and so are set to ``None``.

.. code:: python

    a[0][0].tags




.. parsed-literal::

    {'omega': None,
     'phi': None,
     'psi': -11.577463114977444,
     'tas': (None, None, -11.577463114977444)}



All three torsion angles are defined for the second ``Residue``, and
their values are now stored in its ``tags`` dictionary.

.. code:: python

    a[0][1].tags




.. parsed-literal::

    {'omega': 179.84483742099872,
     'phi': -173.24724224577457,
     'psi': 119.68790084554132,
     'tas': (179.84483742099872, -173.24724224577457, 119.68790084554132)}



If we look at a helical residue in the structure, we'll see the torsion
angles of the :math:`\alpha` helix with which we are familiar.

.. code:: python

    a[0][20].tags




.. parsed-literal::

    {'omega': 178.37366522259435,
     'phi': -66.08456374442004,
     'psi': -39.42791691118421,
     'tas': (178.37366522259435, -66.08456374442004, -39.42791691118421)}



Looking at these torsion angles gives us a decent guess that the 21st
``Residue`` of the first ``Polypeptide`` of our ``Assembly``
(``a[0][20]``) is part of an :math:`\alpha` helix.

To confirm this, we could use the built-in tagging function
``tag_secondary_structure``.

4. Tagging secondary structure
------------------------------

In an earlier tutorial, you will have been introduced to the ``helices``
and ``strands`` attributes of ``Assembly`` and ``Polypeptide`` objects.

Underlying these attributes is the tagging function
``tag_secondary_structure``.

Let's look at the first ``Residue`` in the first helix of ``a.helices``.

.. code:: python

    a.helices[0][0]




.. parsed-literal::

    <Residue containing 9 Atoms. Residue code: GLN>



We know that this ``Residue`` is already tagged with its torsion angles.

.. code:: python

    a.helices[0][0].tags




.. parsed-literal::

    {'omega': 177.09798623768685,
     'phi': -75.1600178143756,
     'psi': -27.827696876562598,
     'secondary_structure': 'H',
     'tas': (177.09798623768685, -75.1600178143756, -27.827696876562598)}



But there's an additional tag there too.

| That's because the ``helices`` method calls the function
  ``tag_secondary_structure``.
| You can see this by looking for the ``helices`` property within the
  ``Polypeptide`` class in ``isambard.ampal.protein``.

When we run ``tag_secondary_structure``, each ``Residue`` is tagged with
its secondary structure, as assigned by
`DSSP <http://www.cmbi.ru.nl/dssp.html>`__.

DSSP assigns each residue a secondary structure value using a single
character from the following list:

.. raw:: html

   <li>

S: Bend

.. raw:: html

   </li>

.. raw:: html

   <li>

H: Alpha helix (4-12)

.. raw:: html

   </li>

.. raw:: html

   <li>

I: pi helix

.. raw:: html

   </li>

.. raw:: html

   <li>

T: Turn

.. raw:: html

   </li>

.. raw:: html

   <li>

B: Isolated beta-bridge residue

.. raw:: html

   </li>

.. raw:: html

   <li>

E: Strand

.. raw:: html

   </li>

.. raw:: html

   <li>

G: 3-10 helix

.. raw:: html

   </li>

| A blank character ' ' is used if no secondary structure can be
  assigned.
| I think there is an acronym to help you remember which 7 letters are
  used as one-letter secondary structure assignments.

We call ``tag_secondary_structure`` just like we called
``tag_torsion_angles``.

.. code:: python

    a.tag_secondary_structure()

The ``tags`` dictionary of each ``Residue`` now contains its secondary
structure assignment in addition to its torsion angles.

.. code:: python

    a[0][0].tags




.. parsed-literal::

    {'omega': None,
     'phi': None,
     'psi': -11.577463114977444,
     'secondary_structure': ' ',
     'tas': (None, None, -11.577463114977444)}



| No secondary stricture has been assigned to the first ``Residue``.
| The 21st ``Residue`` of the ``Polypeptide`` is :math:`\alpha`-helical
  though:

.. code:: python

    a[0][20].tags




.. parsed-literal::

    {'omega': 178.37366522259435,
     'phi': -66.08456374442004,
     'psi': -39.42791691118421,
     'secondary_structure': 'H',
     'tas': (178.37366522259435, -66.08456374442004, -39.42791691118421)}



Exercises
~~~~~~~~~

Recall the ``helices`` method, which returns a new ``Assembly`` whose
``Polypeptides`` are the helices of the original ``Assembly``.

.. code:: python

    a.helices




.. parsed-literal::

    <Assembly (2ebo) containing 9 Polypeptides>



We can use the ``tag_secondary_structure`` method to verify that all of
the ``Residues`` in ``a.helices`` are indeed helical.

Read the code in the cell below before running it...

.. code:: python

    for res in a.helices.get_monomers():
        if res.tags['secondary_structure'] != 'H':
            print('{0} is not helical!'.format(res))

1. Make sure you understand what the code in the cell below is doing.
   Why does it not ``print`` anything?

2. What happens when you delete the ``.helices`` part of the first line
   (so that it reads "``for res in a.get_monomers():``"?) Why?

3. Write a list comprehension that will give you a list of all the
   helical residues in an ``Assembly``.

4. Do the exercise above without using the ``helices`` method.

Have a think about what is going on 'behind the scenes' when you run
``tag_secondary_structure``. - The program DSSP is run, using the
``pdb`` attribute of the ``Assembly`` object as its input. - The output
from DSSP is collected and parsed it for its secondary structure
assignments. - These assignments are added back into the ``Assembly``
object in the appropriate place. - Here, the appropriate place is the
``tags`` attributes of the ``Residue`` objects.

| Further, what is going on 'behind the scenes' when you use the
  ``helices`` method?
| - ``tag_secondary_structure`` is run (and therefore all of the steps
  above are carried out). - The ``seconadary_structure`` tag for each
  ``Residue`` is then looked at in turn. - Consecutive ``Residues`` with
  ``secondary_structure = 'H'`` are grouped together. - Each group of
  ``Residues`` is used to instantiate a ``Polypeptide`` object. - These
  ``Polypeptides`` are then used to instantiate a new ``Assembly``
  object. - It is this new ``Assembly`` object that is returned by the
  ``helices`` method.

5. Other tagging functions
--------------------------

Other tagging functions in ISAMBARD follow the pattern of
``tag_torsion_angles`` and ``tag_secondary_structure``. You call them on
your ``Assembly`` or ``Polypeptide`` object, they run some calculations
for you, and add the results to the ``tags`` atrribute of the relevant
object(s) in the AMPAL framework.

Try typing "``a.tag_``" into the code cell below, and then pressing the
``Shift`` key on your keyboard.


You should see a pop-up listing a small number of tagging functions, two
of which we have already explored. We'll cover two more of these very
briefly here, but feel free to look at / play around with the remaining
tagging functions.

tag\_ca\_geometry
~~~~~~~~~~~~~~~~~

Running the ``tag_ca_geometry`` function tags each ``Residue`` of the
``Protein`` with values for 'rise\_per\_residue',
'radius\_of\_curvature' and 'residues\_per\_turn'.

The concepts of 'rise\_per\_residue' and 'residues\_per\_turn' can be
calculated for any region of the protein structure. The
'radius\_of\_curvature' represents how bent or straight a part of the
protein backbone is, with the largest values occurring at the
straightest regions of the structure.

.. code:: python

    a.tag_ca_geometry()

Let's look at the ``tags`` of a helical ``Residue``. The values for
'residues\_per\_turn' and 'rise\_per\_residue' are what we'd expect to
see for something in an :math:`\alpha` helix.

.. code:: python

    a[0][20].tags




.. parsed-literal::

    {'omega': 178.37366522259435,
     'phi': -66.08456374442004,
     'psi': -39.42791691118421,
     'radius_of_curvature': 49.216212356782755,
     'residues_per_turn': 3.631817377474546,
     'rise_per_residue': 1.5240018755917517,
     'secondary_structure': 'H',
     'tas': (178.37366522259435, -66.08456374442004, -39.42791691118421)}



Exercise
~~~~~~~~

-  Take some time to have a look at the source code for each of the
   tagging functions that you saw in the pop-up.

-  Is there a general difference between the tagging functions at the
   ``Assembly`` level and those at the ``Polymer`` level?

-  | Now look at the tagging functions within the ``Polypeptide`` class.
     Notice their similar structure.
   | This similarity is **not** coincidental!
