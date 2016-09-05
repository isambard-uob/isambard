
.. code:: python

    import isambard

Advanced Selections and Analysis
================================

There are lots of tools in ``isambard`` for analysing protein structure,
most of them are baked into the AMPAL objects themselves.

1. Selecting Secondary Structure
--------------------------------

Let's load-up a PDB file:

.. code:: python

    my_protein_3qy1 = isambard.ampal.convert_pdb_to_assembly('3qy1.pdb')

You can select all the helices or strands simply by typing:

.. code:: python

    my_protein_3qy1.helices




.. parsed-literal::

    <Protein containing 20 Chains>



.. code:: python

    my_protein_3qy1.strands




.. parsed-literal::

    <Protein containing 10 Chains>



These attributes return ``Protein`` objects, where each of the
helices/strands are represented as a separate chain. All of the normal
attributes and methods associated with ``Protein`` and ``Chain`` objects
work just the same as before.

.. code:: python

    my_helices_3qy1 = my_protein_3qy1.helices

.. code:: python

    my_helices_3qy1.sequences




.. parsed-literal::

    ['IDTLISNNALWSKMLVEE',
     'GFFEKLA',
     'AERLT',
     'LNCLSVVQ',
     'GGIKAAVE',
     'INNWLLHIRDIWLK',
     'SSLLG',
     'RLDALYELNVMEQVYNLGH',
     'TIMQSAWK',
     'RETLENGYHKGISALSLKY',
     'IDTLISNNALWSKMLVEE',
     'GFFEKLAQ',
     'AERLT',
     'LNCLSVVQ',
     'GGIKAAVE',
     'INNWLLHIRDIWLK',
     'SSLLG',
     'RLDALYELNVMEQVYNLGH',
     'TIMQSAWK',
     'RETLENGYHKGISALSLKY']



.. code:: python

    my_helices_3qy1[0]




.. parsed-literal::

    <Chain containing 18 Residues. Sequence: IDTLISNNALWS...>



.. code:: python

    my_helices_3qy1[0].sequence




.. parsed-literal::

    'IDTLISNNALWSKMLVEE'



.. code:: python

    my_helices_3qy1[0].molecular_weight




.. parsed-literal::

    2076.3696800000002



.. code:: python

    my_helices_3qy1[0][15]




.. parsed-literal::

    <Residue containing 7 Atoms. Residue code: VAL>



.. code:: python

    my_helices_3qy1[0][15]['CA']




.. parsed-literal::

    <Carbon Atom. Coordinates: (23.340, -8.324, -33.970)>



It is worth noting that the AMPAL parent does **not** change when you
select the helices or strands, so ``ampal_parent`` returns the original
``Protein`` that they came from not the selection. This means that you
retain all the original information from the PDB file.

.. code:: python

    my_helix = my_helices_3qy1[2]

.. code:: python

    my_helix.ampal_parent == my_helices_3qy1




.. parsed-literal::

    False



.. code:: python

    my_helix.ampal_parent == my_protein_3qy1




.. parsed-literal::

    True



.. code:: python

    my_helix.ampal_parent




.. parsed-literal::

    <Protein containing 2 Chains>



2. Selecting All Residues or Atoms
----------------------------------

Sometimes it's convinient to select all of the ``Residues`` or ``Atoms``
in a ``Protein`` or ``Chain`` object:

.. code:: python

    my_protein_3qy1.get_monomers()




.. parsed-literal::

    <itertools.chain at 0x111af9128>



.. code:: python

    my_helices_3qy1.get_atoms()




.. parsed-literal::

    <itertools.chain at 0x111af7780>



As you can see an itertools object is returned. This might be slightly
confusing as you might expect a list, but this is what's known as an
``iterator``, you can loop over it like a list or string, but if you
want to use it repeatedly or examine its contents you'll need to convert
it to a list. The advantage of returning an ``iterator`` is that it's
much more memory efficient, and these lists could potentially be very
large. If you'd like to know more about ``iterables`` and ``iterators``,
as well as related objects called ``generators``, see the following link
(this is quite advanced Python and is not essential for using any of the
AMPAL framework: `Iterators and
Generators <http://anandology.com/python-practice-book/iterators.html>`__.

.. code:: python

    list(my_protein_3qy1.get_atoms())[:20]  # Showing the first 20 elements for clarity




.. parsed-literal::

    [<Nitrogen Atom. Coordinates: (14.714, -30.168, -26.423)>,
     <Carbon Atom. Coordinates: (15.518, -30.153, -25.207)>,
     <Carbon Atom. Coordinates: (16.111, -28.769, -24.931)>,
     <Oxygen Atom. Coordinates: (15.960, -27.855, -25.734)>,
     <Carbon Atom. Coordinates: (16.613, -31.220, -25.270)>,
     <Carbon Atom. Coordinates: (16.067, -32.624, -25.153)>,
     <Oxygen Atom. Coordinates: (14.899, -32.777, -24.743)>,
     <Oxygen Atom. Coordinates: (16.807, -33.576, -25.474)>,
     <Nitrogen Atom. Coordinates: (16.782, -28.637, -23.789)>,
     <Carbon Atom. Coordinates: (17.360, -27.364, -23.339)>,
     <Carbon Atom. Coordinates: (18.466, -26.876, -24.299)>,
     <Oxygen Atom. Coordinates: (18.513, -25.674, -24.586)>,
     <Carbon Atom. Coordinates: (17.842, -27.509, -21.862)>,
     <Carbon Atom. Coordinates: (16.643, -27.442, -20.889)>,
     <Carbon Atom. Coordinates: (18.956, -26.516, -21.456)>,
     <Carbon Atom. Coordinates: (15.873, -26.033, -20.777)>,
     <Nitrogen Atom. Coordinates: (19.310, -27.788, -24.836)>,
     <Carbon Atom. Coordinates: (20.370, -27.411, -25.786)>,
     <Carbon Atom. Coordinates: (19.801, -26.707, -27.030)>,
     <Oxygen Atom. Coordinates: (20.438, -25.781, -27.539)>]



.. code:: python

    list(my_protein_3qy1.get_monomers())[:20]  # Showing the first 20 elements for clarity




.. parsed-literal::

    [<Residue containing 8 Atoms. Residue code: ASP>,
     <Residue containing 8 Atoms. Residue code: ILE>,
     <Residue containing 8 Atoms. Residue code: ASP>,
     <Residue containing 7 Atoms. Residue code: THR>,
     <Residue containing 8 Atoms. Residue code: LEU>,
     <Residue containing 8 Atoms. Residue code: ILE>,
     <Residue containing 6 Atoms. Residue code: SER>,
     <Residue containing 8 Atoms. Residue code: ASN>,
     <Residue containing 8 Atoms. Residue code: ASN>,
     <Residue containing 5 Atoms. Residue code: ALA>,
     <Residue containing 8 Atoms. Residue code: LEU>,
     <Residue containing 14 Atoms. Residue code: TRP>,
     <Residue containing 6 Atoms. Residue code: SER>,
     <Residue containing 9 Atoms. Residue code: LYS>,
     <Residue containing 8 Atoms. Residue code: MET>,
     <Residue containing 8 Atoms. Residue code: LEU>,
     <Residue containing 7 Atoms. Residue code: VAL>,
     <Residue containing 9 Atoms. Residue code: GLU>,
     <Residue containing 9 Atoms. Residue code: GLU>,
     <Residue containing 8 Atoms. Residue code: ASP>]



3. Analysing Composition
------------------------

We can easily look at the composition of sequences and structures using
a ``Counter`` object. ``Counter`` can be fed any iterable (``lists`` and
``strings`` are the most commonly used) and will count the occurence of
each element inside. We can start by looking at the composition of amino
acids in a sequence:

.. code:: python

    from collections import Counter

.. code:: python

    my_protein_3qy1['A'].sequence




.. parsed-literal::

    'DIDTLISNNALWSKMLVEEDPGFFEKLAQAQKPRFLWIGCSDSRVPAERLTGLEPGELFVHRNVANLVIHTDLNCLSVVQYAVDVLEVEHIIICGHSGCGGIKAAVENPELGLINNWLLHIRDIWLKHSSLLGKMPEEQRLDALYELNVMEQVYNLGHSTIMQSAWKRGQNVTIHGWAYSINDGLLRDLDVTATNRETLENGYHKGISALSLKYI'



.. code:: python

    Counter(my_protein_3qy1['A'].sequence)




.. parsed-literal::

    Counter({'A': 13,
             'C': 4,
             'D': 11,
             'E': 16,
             'F': 4,
             'G': 16,
             'H': 9,
             'I': 16,
             'K': 9,
             'L': 29,
             'M': 4,
             'N': 14,
             'P': 6,
             'Q': 7,
             'R': 9,
             'S': 13,
             'T': 8,
             'V': 15,
             'W': 6,
             'Y': 6})



But as stated before, you can use Counters on any iterable, not just
strings. Let's make a list of all the pdb molecule codes of the ligands:

.. code:: python

    my_ligands = my_protein_3qy1.ligands

.. code:: python

    my_ligands[0]




.. parsed-literal::

    <Ligand containing 1 Atom. Ligand code: ZN>



.. code:: python

    my_ligands[0].mol_code  # This can be used to find the pdb molecule code of any residue




.. parsed-literal::

    'ZN'



There are two ways to generate a list of the codes, a ``for`` loop or a
``list`` comprehension, use whichever you are comfortable with. If you'd
like to know more about ``list`` comprehensions, please see the
following link (this is relatively advanced Python but while not
essential for using any of the AMPAL framework, it is very useful: `List
Comprehensions <https://docs.python.org/3.5/tutorial/datastructures.html>`__
(scroll down to the relevant section).

.. code:: python

    # With a for loop
    mol_codes_1 = []
    for lig in my_ligands:
        mol_codes_1.append(lig.mol_code)

.. code:: python

    mol_codes_1[:5]  # The first 5 elements




.. parsed-literal::

    ['ZN', 'ZN', 'HOH', 'HOH', 'HOH']



.. code:: python

    # A list comprehension
    mol_codes_2 = [lig.mol_code for lig in my_ligands]

.. code:: python

    mol_codes_2[:5]  # Showing the first 5 elements for clarity




.. parsed-literal::

    ['ZN', 'ZN', 'HOH', 'HOH', 'HOH']



You can use either of these methods, use whichever one you're more
comfortable with.

.. code:: python

    mol_codes_1 == mol_codes_2  # The lists that re produced are exactly the same




.. parsed-literal::

    True



Now the ``list`` of mol codes can be used to make a ``Counter`` object:

.. code:: python

    Counter(mol_codes_1)




.. parsed-literal::

    Counter({'HOH': 447, 'ZN': 2})



As you can see, there are 447 water molecules and 2 zinc ions.

4. Distance Analysis
--------------------

Now we can select all atoms in the protein and understand the structures
composition, we can perform some simple analysis. Let's try and find all
the residues that are close to the zinc ions.

.. code:: python

    zinc_1 = my_ligands[0]

.. code:: python

    zinc_1




.. parsed-literal::

    <Ligand containing 1 Atom. Ligand code: ZN>



All ligands are residues, even if they only contain a single atom. So we
use the zinc ``Atom`` itself when measuring distances:

.. code:: python

    zinc_1['ZN']




.. parsed-literal::

    <Zinc Atom. Coordinates: (-5.817, -20.172, -18.798)>



Measuring distances is simple, you can use the
``isambard.geometry.distance`` function. It takes two 3D vectors as an
input, these can be in list form, tuples, or even ``Atom`` objects:

.. code:: python

    isambard.geometry.distance(zinc_1['ZN'], (0, 0, 0))  # Distance from the origin




.. parsed-literal::

    28.179990720367528



.. code:: python

    first_ca = my_protein_3qy1['A'][0]['CA']  # CA of the first residue in chain A

.. code:: python

    first_ca




.. parsed-literal::

    <Carbon Atom. Coordinates: (15.518, -30.153, -25.207)>



.. code:: python

    isambard.geometry.distance(zinc_1['ZN'], first_ca)  # Distance in angstroms




.. parsed-literal::

    24.41060972200408



Now we need to loop over all the atoms and find which are close (<= 3 Å)
to the zinc. We can use the distance function in geometry to do this:

.. code:: python

    atoms_close_to_zinc = []
    for at in my_protein_3qy1.get_atoms():
        if isambard.geometry.distance(zinc_1['ZN'], at) <= 3.0:
            atoms_close_to_zinc.append(at)

.. code:: python

    atoms_close_to_zinc




.. parsed-literal::

    [<Sulfur Atom. Coordinates: (-4.322, -18.933, -17.640)>,
     <Oxygen Atom. Coordinates: (-4.771, -22.057, -19.213)>,
     <Nitrogen Atom. Coordinates: (-6.209, -19.569, -20.787)>,
     <Sulfur Atom. Coordinates: (-7.753, -20.619, -17.709)>]



There are 4 atoms within 3 Å of the zinc, 2 sulphur atoms, an oxygen and
a nitrogen. Let's find the residues that are coordinating the zinc:

.. code:: python

    atoms_close_to_zinc[0]




.. parsed-literal::

    <Sulfur Atom. Coordinates: (-4.322, -18.933, -17.640)>



.. code:: python

    atoms_close_to_zinc[0].ampal_parent




.. parsed-literal::

    <Residue containing 6 Atoms. Residue code: CYS>



We can get all the residues using a for loop or a list comprehension:

.. code:: python

    residues_close_to_zinc = []
    for at in atoms_close_to_zinc:
        residues_close_to_zinc.append(at.ampal_parent)

.. code:: python

    residues_close_to_zinc




.. parsed-literal::

    [<Residue containing 6 Atoms. Residue code: CYS>,
     <Residue containing 8 Atoms. Residue code: ASP>,
     <Residue containing 10 Atoms. Residue code: HIS>,
     <Residue containing 6 Atoms. Residue code: CYS>]



.. code:: python

    # The list comprehension is much more concise
    residues_close_to_zinc_2 = [at.ampal_parent for at in atoms_close_to_zinc]

.. code:: python

    residues_close_to_zinc_2




.. parsed-literal::

    [<Residue containing 6 Atoms. Residue code: CYS>,
     <Residue containing 8 Atoms. Residue code: ASP>,
     <Residue containing 10 Atoms. Residue code: HIS>,
     <Residue containing 6 Atoms. Residue code: CYS>]



It looks like the zinc is coordinated by two cysteines, an aspartate and
a histidine residue.

5. Is Within
------------

This kind of operation is very common when analysing proteins. So we
have some built-in methods for handling this on ``Protein`` and
``Chain`` objects:

.. code:: python

    my_protein_3qy1.is_within(3, zinc_1['ZN'])  # It takes a distance and a point




.. parsed-literal::

    [<Sulfur Atom. Coordinates: (-4.322, -18.933, -17.640)>,
     <Oxygen Atom. Coordinates: (-4.771, -22.057, -19.213)>,
     <Nitrogen Atom. Coordinates: (-6.209, -19.569, -20.787)>,
     <Sulfur Atom. Coordinates: (-7.753, -20.619, -17.709)>]



.. code:: python

    my_protein_3qy1.is_within(3, (-10, -10, -10))




.. parsed-literal::

    [<Nitrogen Atom. Coordinates: (-7.278, -11.112, -10.445)>,
     <Carbon Atom. Coordinates: (-7.525, -10.089, -8.339)>,
     <Oxygen Atom. Coordinates: (-8.676, -9.655, -8.177)>,
     <Carbon Atom. Coordinates: (-11.364, -9.364, -12.255)>,
     <Carbon Atom. Coordinates: (-10.337, -10.223, -12.972)>,
     <Carbon Atom. Coordinates: (-10.752, -12.047, -11.188)>,
     <Oxygen Atom. Coordinates: (-10.046, -11.618, -10.265)>,
     <Nitrogen Atom. Coordinates: (-11.667, -9.798, -7.931)>]



.. code:: python

    my_protein_3qy1['A'].is_within(3, zinc_1['ZN'])




.. parsed-literal::

    [<Sulfur Atom. Coordinates: (-4.322, -18.933, -17.640)>,
     <Oxygen Atom. Coordinates: (-4.771, -22.057, -19.213)>,
     <Nitrogen Atom. Coordinates: (-6.209, -19.569, -20.787)>,
     <Sulfur Atom. Coordinates: (-7.753, -20.619, -17.709)>]



.. code:: python

    my_protein_3qy1['B'].is_within(3, zinc_1['ZN']) 
    # This list is empty as nothing on chain B is close to the zinc




.. parsed-literal::

    []



There is a partner method to is ``is_within``, every monomer (this
includes ``Residues`` and ``Ligands``) has an ``environment`` method.
This returns all ``Residues`` within a given cutoff value. This means
that you can either select all atoms within a distance of a point for a
given ``Protein`` or ``Chain`` or find the residues that surround a
particular ``Residue`` or ``Ligand``.

.. code:: python

    zinc_1.environment()




.. parsed-literal::

    [<Residue containing 6 Atoms. Residue code: CYS>,
     <Residue containing 8 Atoms. Residue code: ASP>,
     <Residue containing 10 Atoms. Residue code: HIS>,
     <Residue containing 6 Atoms. Residue code: CYS>]



.. code:: python

    zinc_1.environment(cutoff=6)




.. parsed-literal::

    [<Residue containing 6 Atoms. Residue code: CYS>,
     <Residue containing 6 Atoms. Residue code: SER>,
     <Residue containing 8 Atoms. Residue code: ASP>,
     <Residue containing 6 Atoms. Residue code: SER>,
     <Residue containing 11 Atoms. Residue code: ARG>,
     <Residue containing 7 Atoms. Residue code: VAL>,
     <Residue containing 5 Atoms. Residue code: ALA>,
     <Residue containing 8 Atoms. Residue code: ASN>,
     <Residue containing 10 Atoms. Residue code: HIS>,
     <Residue containing 4 Atoms. Residue code: GLY>,
     <Residue containing 6 Atoms. Residue code: CYS>,
     <Residue containing 4 Atoms. Residue code: GLY>,
     <Residue containing 4 Atoms. Residue code: GLY>,
     <Residue containing 8 Atoms. Residue code: ILE>]



.. code:: python

    zinc_1.environment(include_self=True)




.. parsed-literal::

    [<Residue containing 6 Atoms. Residue code: CYS>,
     <Residue containing 8 Atoms. Residue code: ASP>,
     <Residue containing 10 Atoms. Residue code: HIS>,
     <Residue containing 6 Atoms. Residue code: CYS>,
     <Ligand containing 1 Atom. Ligand code: ZN>]



6. Geometry in ISAMBARD
-----------------------

There are a range of tools in ISAMBARD for performing geometric
operations. We've already covered distance, but other commonly used
functions include ``angle_between_vectors``, ``dihedral``,
``unit_vector``, ``find_foot``, ``radius_of_circumcircle``. Be sure to
check out the source code if you need a specific geometric function or
have a look through the documentation.

The ``dihedral`` function is probably the most useful of these for
analysing proteins, so let's use it to measure some torsion angles. It
requires 4 3D vectors to calculate the dihedral, again these can be
``lists``, ``tuples``, ``numpy.arrays`` or ``Atoms``. **Note:** This
method of calculating torsion angles is only as an example, see the
Tagging tutorial for the proper, low-effort method!

.. code:: python

    r1 = my_protein_3qy1['B'][4]
    r2 = my_protein_3qy1['B'][5]
    r3 = my_protein_3qy1['B'][6]

.. code:: python

    omega = isambard.geometry.dihedral(r1['CA'], r1['C'], r2['N'], r2['CA'])

.. code:: python

    phi = isambard.geometry.dihedral(r1['C'], r2['N'], r2['CA'], r2['C'])

.. code:: python

    psi = isambard.geometry.dihedral(r2['N'], r2['CA'], r2['C'], r3['N'])

.. code:: python

    print(omega, phi, psi)


.. parsed-literal::

    -179.800200103 -64.1086116305 -45.8437396848


We can use it to calculate the :math:`\chi` torsion angles too. R2 is
leucine, so we can calculate the :math:`\chi_1` and :math:`\chi_2`
angles:

.. code:: python

    r2.atoms




.. parsed-literal::

    OrderedDict([('N', <Nitrogen Atom. Coordinates: (-5.186, -2.004, -31.807)>),
                 ('CA', <Carbon Atom. Coordinates: (-4.911, -3.362, -31.310)>),
                 ('C', <Carbon Atom. Coordinates: (-5.985, -4.346, -31.786)>),
                 ('O', <Oxygen Atom. Coordinates: (-5.650, -5.434, -32.255)>),
                 ('CB', <Carbon Atom. Coordinates: (-4.788, -3.418, -29.770)>),
                 ('CG', <Carbon Atom. Coordinates: (-3.838, -2.437, -29.061)>),
                 ('CD1', <Carbon Atom. Coordinates: (-3.653, -2.831, -27.613)>),
                 ('CD2', <Carbon Atom. Coordinates: (-2.478, -2.359, -29.736)>)])



.. code:: python

    chi1 = isambard.geometry.dihedral(r2['N'], r2['CA'], r2['CB'], r2['CG'])

.. code:: python

    chi2 = isambard.geometry.dihedral(r2['CA'], r2['CB'], r2['CG'], r2['CD1'])

.. code:: python

    print(chi1, chi2)


.. parsed-literal::

    -51.3804036635 -168.700062066


Our simple analysis shows that the leucine residue is in the
gauche-/trans conformation.

7. Summary and activities
-------------------------

There are lots of tools for making complex selections in ISAMBARD.
Combined with the tools for geometry, detailed analysis can be performed
on these selections.

1. Find all the residues that are:

   1. within 5 Å of crystal water.
   2. *not* within 5 Å of crystal water.

2. Find how many cis-peptide bonds there are in this structure.
3. Perform these activities on another PDB file.
