
Tagging in the AMPAL framework
==============================

1. Introduction
---------------

Tagging allows us to store information related to AMPAL objects and
retrieve it easily later.

All AMPAL objects have a ``tags`` attribute, which is a dictionary: a
'useful place to put stuff'. When we add items to the ``tags``
dictionary, we refer to it as 'tagging' the AMPAL object. Any python
object can be stored in ``tags``.

You've already done some tagging. The ``tags`` attribute is also used to
store information internally so that certain values only need to be
calculated once.

It may not have been pointed out at the time, but you have already been
using ``tags`` in previous tutorial notebooks, when you used properties
such as ``helices`` and ``bude_score``.

When you use ``convert_pdb_to_assembly``, you tag the ``Atoms`` of your
``Protein`` with information from the pdb file. We'll look at that
shortly.

| In-built functions: ``isambard.ampal.base_ampal`` has some in-built
  'tagging functions' that calculate things and add the results to the
  ``tags`` of the appropriate AMPAL object.
| In this notebook, we'll look at some of these and (hopefully!)
  demonstrate their utility.

2. Getting started
------------------

First, let's import some functions from ``isambard`` for making a
``Protein`` object.

We'll use the function ``get_mmol`` from the ``filesystem`` module. If
you've not met this yet, don't worry - it's introduced fully in the
``IsmabardFileSystem`` tutorial. It's just a convenient way of
downloading the contents of a pdb file without having to actually write
any files.

.. code:: python

    from isambard.add_ons.filesystem import get_mmol
    from isambard.ampal.base_ampal import convert_pdb_to_assembly

We'll use the structure
`2ebo <http://www.ebi.ac.uk/pdbe/entry/pdb/2ebo>`__ for this. It's a
viral protein containing coiled-coil packing that is involved in fusing
the membranes of the ebola virus and its host.

.. code:: python

    code = '2ebo'
    a = convert_pdb_to_assembly(get_mmol(code), path=False, pdb_name=code)

Our new ``Protein`` object ``a`` contains 3 Chains, and its ``id``
attribute is the pdb code we used to create it.

.. code:: python

    a




.. parsed-literal::

    <Protein containing 3 Chains>



.. code:: python

    a.id




.. parsed-literal::

    '2ebo'



It comes with a ``tags`` attribute, which is an empty dictionary.

.. code:: python

    a.tags




.. parsed-literal::

    {}



The 3 Chains in our Protein are also AMPAL objects ...

.. code:: python

    a[0] # The first Chain of the Protein




.. parsed-literal::

    <Chain containing 74 Residues. Sequence: GLRQLANETTQA...>



... and therefore have their own ``tags``:

.. code:: python

    a[0].tags




.. parsed-literal::

    {}



Similarly, all of the Residues in the Chain are AMPAL objects with empty
``tags``.

.. code:: python

    a[0][0] # The first Residue of the first chain of the Protein




.. parsed-literal::

    <Residue containing 4 Atoms. Residue code: GLY>



.. code:: python

    a[0][0].tags




.. parsed-literal::

    {}



Each ``Residue`` in our ``Protein`` is made from ``Atoms``. Right down
to the ``Atoms``, AMPAL objects have ``tags``.

Let's look at the ``tags`` of the backbone Nitrogen ``Atom`` of the
first ``Residue`` of the first ``Chain`` of our ``Protein``.

.. code:: python

    a[0][0]['N']




.. parsed-literal::

    <Nitrogen Atom. Coordinates: (-14.780, 25.698, -6.988)>



.. code:: python

    a[0][0]['N'].tags




.. parsed-literal::

    {'bfactor': 71.51, 'charge': '', 'occupancy': 1.0}



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

If you're familiar with adding items to python dictionaries, then you're
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

    {'number_of_atoms': 1806,
     'protein_description': 'viral membrane fusion protein'}



Exercises
~~~~~~~~~

1. Look at the page for 2ebo on the
   `PDB <http://www.rcsb.org/pdb/explore/explore.do?structureId=2ebo>`__
   or `PDBE <http://www.ebi.ac.uk/pdbe/entry/pdb/2ebo>`__. Find the
   resolution of the structure and add that to ``tags``.
2. Add another tag to the ``Protein`` object storing its
   '``number_of_residues``'.
   HINT: Use a similar expression to one used for ``number_of_atoms``.
3. Tag each ``Chain`` in the ``Protein`` with a '``number_of_atoms``'
   tag, like we did with the whole ``Protein`` earlier.
4. Tag each ``Residue`` in Chain A with a '``number_of_atoms``' tag.
5. Tag Chain B with a 'number\_of\_tryptophans' tag.
   HINT: You may want to use the ``Counter`` class you met earlier in
   order to do this. Use "``from collections import Counter``" to bring
   this class into your namespace.

After completing Exercise 1. and 2., you should have 4 items in your
``a.tags`` dictionary.

.. code:: python

    a.tags




.. parsed-literal::

    {'number_of_atoms': 1806,
     'protein_description': 'viral membrane fusion protein'}



Let's add one more thing.

In the MakingModels tutorial, you will have used the ``bude_score``
attribute.

.. code:: python

    a.bude_score




.. parsed-literal::

    -3265.4395000000004



The bude\_score method automatically adds the output to tags:

.. code:: python

    a.tags




.. parsed-literal::

    {'bude_score': -3265.4395000000004,
     'number_of_atoms': 1806,
     'protein_description': 'viral membrane fusion protein'}



This makes it easy to retrieve later without re-running bude.

.. code:: python

    a.tags['bude_score']




.. parsed-literal::

    -3265.4395000000004



3. Tagging torsion angles
-------------------------

In an earlier tutorial we used the ``isambard.tools.geometry.dihedral``
function to calculate the backbone torsion angles of a ``Residue`` in
the AMPAL framework.

| Suppose we wanted to calculate all of the torsion angles in a
  ``Protein``.
| We could loop over all of the ``Residues`` in the ``Protein`` using
  the ``get_monomers()`` method, call the
  ``isambard.tools.geometry.dihedral`` function at each stage of the
  loop and store the results somewhere convenient.

| This would be perfectly valid.
| However, we also have an in-built method for doing this.

The ``tag_torsion_angles`` method calculates the torsion angles for each
``Residue`` in a ``Protein`` (or ``Chain``) and adds them ``tags``
dictionary of the ``Residue``.

.. code:: python

    a.tag_torsion_angles()

The ``.tags`` dictionary for the first Residue now contains values for
its ``omega``, ``phi`` and ``psi`` angles.

Since it's the first ``Residue`` of the ``Chain``, its ``omega`` and
``phi`` torsion angles are not defined - hence the ``'nan'`` values
(**n**\ ot **a** **n**\ umber).

.. code:: python

    a[0][0].tags




.. parsed-literal::

    {'omega': 'nan', 'phi': 'nan', 'psi': -11.577463114977443}



All three torsion angles are defined for the second ``Residue``, and
their values are now stored in its ``tags`` dictionary.

.. code:: python

    a[0][1].tags




.. parsed-literal::

    {'omega': 179.84483742099872,
     'phi': -173.24724224577466,
     'psi': 119.68790084554134}



If we look at a helical residue in the structure, we'll see the torsion
angles of the :math:`\alpha`-helix with which we are familiar.

.. code:: python

    a[0][20].tags




.. parsed-literal::

    {'omega': 178.37366522259435,
     'phi': -66.084563744420052,
     'psi': -39.427916911184212}



Looking at these torsion angles gives us a decent guess that the 21st
``Residue`` of the first ``Chain`` of our ``Protein`` (``a[0][20]``) is
part of an :math:`\alpha`-helix.

To confirm this, we could use the in-built tagging function
``tag_secondary_structure``.

4. Tagging secondary structure
------------------------------

In an earlier tutorial, you will have been introduced to the ``helices``
and ``strands`` attributes of ``Protein`` and ``Chain`` objects.

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

    {'omega': 177.09798623768697,
     'phi': -75.160017814375607,
     'psi': -27.827696876562612,
     'secondary_structure': 'H'}



But there's an additional tag there too.

| That's because the ``helices`` method calls the function
  ``tag_secondary_structure``.
| You can see this by looking for the ``helices`` property within the
  ``Chain`` class of ``isambard.amapl.base_ampal``.

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

    {'omega': 'nan',
     'phi': 'nan',
     'psi': -11.577463114977443,
     'secondary_structure': ' '}



| No secondary stricture has been assigned to the first ``Residue``.
| The 21st ``Residue`` of the ``Chain`` is :math:`\alpha`-helical
  though:

.. code:: python

    a[0][20].tags




.. parsed-literal::

    {'omega': 178.37366522259435,
     'phi': -66.084563744420052,
     'psi': -39.427916911184212,
     'secondary_structure': 'H'}



Exercises
~~~~~~~~~

Recall the ``helices`` method, which returns a new ``Protein`` whose
``Chains`` are the helices of the original ``Protein``.

.. code:: python

    a.helices




.. parsed-literal::

    <Protein containing 9 Chains>



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
   helical residues in a ``Protein``.

4. Do the exercise above without using the ``helices`` method.

Have a think about what is going on 'behind the scenes' when you run
``tag_secondary_structure``. - The program DSSP is run, using the
``pdb`` attribute of the ``Protein`` object as its input. - The output
from DSSP is collected and parsed it for its secondary structure
assignments. - These assignments are added back into the ``Protein``
object in the appropriate place. - Here, the appropriate place is the
``tags`` attributes of the ``Residue`` objects.

| Further, what is going on 'behind the scenes' when you use the
  ``helices`` method?
| - ``tag_secondary_structure`` is run (and therefore all of the steps
  above are carried out). - The ``seconadary_structure`` tag for each
  ``Residue`` is then looked at in turn. - Consecutive ``Residues`` with
  ``secondary_structure = 'H'`` are grouped together. - Each group of
  ``Residues`` is used to instantiate a ``Chain`` object. - These
  ``Chains`` are then used to instantiate a new ``Protein`` object. - It
  is this new ``Protein`` object that is returned by the ``helices``
  method.

5. Other tagging functions
--------------------------

Other tagging functions in ``ismabard`` follow the pattern of
``tag_torsion_angles`` and ``tag_secondary_structure``. You call them on
your ``Protein`` or ``Chain`` object, they run some calculations for
you, and add the results to the ``tags`` atrribute of the relevent
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

The concepts of 'rise\_per\_residue' and 'residues\_per\_turn' should be
familiar from Dek's tutorials on :math:`\alpha`-helical geometry. Their
definitions are generic and so these values can be calculated for any
region of the protein structure. The 'radius\_of\_curvature' represents
how bent or straight a part of the protein backbone is, with the largest
values occurring at the straightest regions of the ``Protein``
structure.

.. code:: python

    a.tag_ca_geometry()

Let's look at the ``tags`` of a helical ``Residue``. The values for
'residues\_per\_turn' and 'rise\_per\_residue' are what we'd expect to
see for something in an :math:`\alpha`-helix.

.. code:: python

    a[0][20].tags




.. parsed-literal::

    {'omega': 178.37366522259435,
     'phi': -66.084563744420052,
     'psi': -39.427916911184212,
     'radius_of_curvature': 49.216212356737174,
     'residues_per_turn': 3.6318173774745457,
     'rise_per_residue': 1.5240018755917517,
     'secondary_structure': 'H'}



tag\_socket
~~~~~~~~~~~

The ``tag_socket`` method runs socket to find knob-into-hole iteractions
within a protein structure. It then tags any knob residues with
information relating to the knob-into-hole interaction that they are
involved in. This information is packaged up into a dictionary called
'knob\_data'.

.. code:: python

    a.tag_socket()

We can find the knob residues using a list comprehension or a ``for``
loop and looking for the knob packages. :-)

First, we'll use a ``for`` loop to get a list of knob residues, called
``knobs``.

.. code:: python

    knobs = []
    for res in a.get_monomers():
        if 'knob_data' in res.tags.keys():
            knobs.append(res)

Let's have a look at the ``.tags`` of the first knob residue.

.. code:: python

    knobs[0].tags




.. parsed-literal::

    {'knob_data': {'h0': ('C', (' ', '622', ' ')),
      'h1': ('C', (' ', '625', ' ')),
      'h2': ('C', (' ', '626', ' ')),
      'h3': ('C', (' ', '629', ' ')),
      'helix': 0,
      'hole_helix': 8,
      'knob_type': 6,
      'max_cv_dist': 5.717654938871355,
      'packing_angle': 81.573},
     'omega': 179.58180595410454,
     'phi': -71.132872529577284,
     'psi': -34.858956625353734,
     'radius_of_curvature': 58.342097116582842,
     'residues_per_turn': 3.6118407776154617,
     'rise_per_residue': 1.5189562060319903,
     'secondary_structure': 'H'}



As stated above, the ``knob_data`` dictionary contains information about
the knob-into-hole interaction for which this ``Residue`` is the knob.
Don't worry if the things in ``knob_data`` are unfamiliar (that's not
the point of this tutorial) - just have a look to get the idea of the
things that can be stored in the ``tags`` dictionary.

| \*\* A(nother) list comprehension example \*\*
| Now, let's created the ``knobs`` list in a different way - using a
  list comprehension. We'll call this list knobs\_2.

.. code:: python

    knobs_2 = [res for res in a.get_monomers() if 'knob_data' in res.tags.keys()]

The ``==`` operator returns a value of ``True`` if the terms on either
side of it evaluate to the same thing, and ``False`` otherwise.

We can use it to show that ``knobs`` and ``knobs_2`` are equivalent.

.. code:: python

    knobs == knobs_2




.. parsed-literal::

    True



Exercise
~~~~~~~~

-  Take some time to have a look at the source code for each of the
   tagging functions that you saw in the pop-up. They are all located in
   the ``ismabard.ampal.base_ampal`` module, and they in *both* the
   ``Protein`` class and the ``Chain`` class.

-  Is there a general difference between the tagging functions in
   ``Protein`` and those in ``Chain``?

-  | Now look at the tagging functions within the ``Chain`` class.
     Notice their similar structure.
   | This similarity is **not** coincidental!

