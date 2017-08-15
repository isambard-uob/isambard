
An Introduction to the AMPAL Framework
======================================

ISAMBARD represents proteins with Python objects that are part of the
AMPAL module. AMPAL stands for Atom, Monomer, Polymer, Assembly, Ligand;
which are the main functional elements of the framework. These objects
are designed to represent the hierarchical nature of proteins in an
intuitive way.

All protein models produced by ISAMBARD are AMPAL objects, but you can
also load in crystal structures.

    **Q. So why aren't these objects just called Protein, Chain and
    Residue?**

    Well, there are ``Polypeptide``, and ``Residue`` objects, but they
    are just protein-specific versions of ``Polymer`` and ``Monomer``.
    We wanted to keep the base objects as generic as possible to allow
    other biomolecules (like DNA or RNA), or even unnatural polymers
    (like :math:`\beta`-amino acids) to be represented using this
    architecture. While some these features are not currently
    implemented in ISAMBARD, this will lead to more scalable code with
    reduced duplication in the future.

1. Converting PDB files to AMPAL Objects
----------------------------------------

.. code:: python

    import isambard

Any PDB file can be parsed into an AMPAL object, which allows you to
easily analyse the structure. The only function you need for this is
``isambard.ampal.convert_pdb_to_ampal``. It takes a file path string as
the input argument:

.. code:: python

    isambard.ampal.convert_pdb_to_ampal('3qy1.pdb')




.. parsed-literal::

    <Assembly (3qy1) containing 2 Polypeptides, 449 Ligands>



The object that is returned is an ``Assembly``. We can assign it to a
variable and look inside it.

.. code:: python

    my_protein = isambard.ampal.convert_pdb_to_ampal('3qy1.pdb')

Remember if you are using the object in Jupyter Notebook once it's
assigned to a variable, you can have a look at its attributes and
methods by typing ``my_protein.`` and then pressing tab.

2. Basic Analysis
-----------------

This ``Assembly`` contains two ``Polypeptides``, and 449 ``Ligands``. It
is easy to check the amino acid sequences of the ``Polypeptides``:

.. code:: python

    my_protein.sequences




.. parsed-literal::

    ['DIDTLISNNALWSKMLVEEDPGFFEKLAQAQKPRFLWIGCSDSRVPAERLTGLEPGELFVHRNVANLVIHTDLNCLSVVQYAVDVLEVEHIIICGHSGCGGIKAAVENPELGLINNWLLHIRDIWLKHSSLLGKMPEEQRLDALYELNVMEQVYNLGHSTIMQSAWKRGQNVTIHGWAYSINDGLLRDLDVTATNRETLENGYHKGISALSLKYI',
     'KDIDTLISNNALWSKMLVEEDPGFFEKLAQAQKPRFLWIGCSDSRVPAERLTGLEPGELFVHRNVANLVIHTDLNCLSVVQYAVDVLEVEHIIICGHSGCGGIKAAVENPELGLINNWLLHIRDIWLKHSSLLGKMPEEQRLDALYELNVMEQVYNLGHSTIMQSAWKRGQNVTIHGWAYSINDGLLRDLDVTATNRETLENGYHKGISALSLKYI']



The ``.sequences`` attribute is a list of sequence strings, one for each
``Polymer``. We can determine other basic properties of the
``Assembly``:

.. code:: python

    my_protein.molecular_weight




.. parsed-literal::

    48508.93158



.. code:: python

    my_protein.molar_extinction_280




.. parsed-literal::

    83640



.. code:: python

    my_protein.isoelectric_point




.. parsed-literal::

    5.4000000000000039



.. code:: python

    my_protein.id




.. parsed-literal::

    '3qy1'



3. Selecting Chains
-------------------

Items inside each ``Assembly`` object can be accessed analogously to
accessing items in a standard Python list:

.. code:: python

    my_protein[0]  # The first chain




.. parsed-literal::

    <Polypeptide containing 215 Residues. Sequence: DIDTLISNNALW...>



You can also select a ``Polymer`` using a string of the chain id from
the PDB file. In this case there are two chains 'A' and 'B'.

.. code:: python

    my_protein['A']




.. parsed-literal::

    <Polypeptide containing 215 Residues. Sequence: DIDTLISNNALW...>



.. code:: python

    my_protein['B']




.. parsed-literal::

    <Polypeptide containing 216 Residues. Sequence: KDIDTLISNNAL...>



The ``Polypeptide`` object has a lot of the same functionality as the
``Assembly``:

.. code:: python

    my_chain_a = my_protein['A']

.. code:: python

    my_chain_a.molecular_weight




.. parsed-literal::

    24199.38728



.. code:: python

    my_chain_a.molar_extinction_280




.. parsed-literal::

    41820



.. code:: python

    my_chain_a.isoelectric_point




.. parsed-literal::

    5.4000000000000039



.. code:: python

    my_chain_a.id




.. parsed-literal::

    'A'



4. Selecting Residues
---------------------

Each ``Polypeptide`` object is made from one or more ``Residue``
objects. You can access the ``Residues`` using square brackets:

.. code:: python

    my_chain_a[0]




.. parsed-literal::

    <Residue containing 8 Atoms. Residue code: ASP>



.. code:: python

    my_chain_a[4]




.. parsed-literal::

    <Residue containing 8 Atoms. Residue code: LEU>



.. code:: python

    my_chain_a[20]




.. parsed-literal::

    <Residue containing 7 Atoms. Residue code: PRO>



You can use a string of a residue id from the PDB file to select a
``Residue``:

.. code:: python

    my_chain_a['23']




.. parsed-literal::

    <Residue containing 7 Atoms. Residue code: PRO>



.. code:: python

    my_chain_a['40']




.. parsed-literal::

    <Residue containing 8 Atoms. Residue code: ILE>



If you use a residue number that isn't defined in the PDB a ``KeyError``
will be raised:

.. code:: python

    my_chain_a['2']


::


    ---------------------------------------------------------------------------

    KeyError                                  Traceback (most recent call last)

    <ipython-input-22-e01d62d6ac3a> in <module>()
    ----> 1 my_chain_a['2']
    

    /Users/jackheal/Projects/isambard/isambard/ampal/protein.py in __getitem__(self, item)
        169         if isinstance(item, str):
        170             id_dict = {str(m.id): m for m in self._monomers}
    --> 171             return id_dict[item]
        172         elif isinstance(item, int):
        173             return self._monomers[item]


    KeyError: '2'


.. code:: python

    my_residue_A23 = my_chain_a['23']

``Residues`` contain an ``OrderedDict`` (a special type of
``dictionary`` that retains the order you add elements) which has atom
identifiers and ``Atom`` objects all the atoms that make up the
``Residue``.

.. code:: python

    my_residue_A23.atoms




.. parsed-literal::

    OrderedDict([('N',
                  <Nitrogen Atom (N). Coordinates: (22.124, -4.140, -35.654)>),
                 ('CA',
                  <Carbon Atom (CA). Coordinates: (22.664, -3.954, -34.292)>),
                 ('C', <Carbon Atom (C). Coordinates: (21.911, -2.875, -33.515)>),
                 ('O', <Oxygen Atom (O). Coordinates: (21.863, -2.926, -32.283)>),
                 ('CB',
                  <Carbon Atom (CB). Coordinates: (24.120, -3.555, -34.534)>),
                 ('CG',
                  <Carbon Atom (CG). Coordinates: (24.124, -2.964, -35.917)>),
                 ('CD',
                  <Carbon Atom (CD). Coordinates: (23.118, -3.764, -36.681)>)])



5. Selecting Atoms
------------------

Atoms can be selected using a string of their PDB atom type, for example
the C\ :math:`\alpha` atom of the residue can be selected like this:

.. code:: python

    my_residue_A23['CA']




.. parsed-literal::

    <Carbon Atom (CA). Coordinates: (22.664, -3.954, -34.292)>



.. code:: python

    my_residue_A23['CG']




.. parsed-literal::

    <Carbon Atom (CG). Coordinates: (24.124, -2.964, -35.917)>



.. code:: python

    my_residue_A23['N']




.. parsed-literal::

    <Nitrogen Atom (N). Coordinates: (22.124, -4.140, -35.654)>



.. code:: python

    my_atom_A23ca = my_residue_A23['CA']

The individual coordinates can be selected using square brackets:

.. code:: python

    my_atom_A23ca[0]




.. parsed-literal::

    22.664000000000001



.. code:: python

    my_atom_A23ca[2]




.. parsed-literal::

    -34.292000000000002



Or with the ``x``, ``y`` and ``z`` properties:

.. code:: python

    my_atom_A23ca.x




.. parsed-literal::

    22.664000000000001



.. code:: python

    my_atom_A23ca.y




.. parsed-literal::

    -3.9540000000000002



.. code:: python

    my_atom_A23ca.z




.. parsed-literal::

    -34.292000000000002



The ``Atom`` object contains some useful attributes:

.. code:: python

    my_atom_A23ca.id  # The atom number from the PDB file




.. parsed-literal::

    162



.. code:: python

    my_atom_A23ca.element  # The element of the atom




.. parsed-literal::

    'C'



6. AMPAL Parents
----------------

Hopefully you can see that it's easy to traverse down the AMPAL
framework from ``Assembly`` level to the ``Atom`` level, but it's just
as easy to work your way back up. With any AMPAL object you can use the
``ampal_parent`` attribute to find the AMPAL object that it is contained
inside.

.. code:: python

    my_atom_A23ca.ampal_parent




.. parsed-literal::

    <Residue containing 7 Atoms. Residue code: PRO>



.. code:: python

    my_residue_A23.ampal_parent




.. parsed-literal::

    <Polypeptide containing 215 Residues. Sequence: DIDTLISNNALW...>



.. code:: python

    my_chain_a.ampal_parent




.. parsed-literal::

    <Assembly (3qy1) containing 2 Polypeptides, 449 Ligands>



This attribute returns the original object itself, meaning you can
access all its methods and functions, including its own
``ampal_parent``!

.. code:: python

    my_atom_A23ca.ampal_parent == my_residue_A23




.. parsed-literal::

    True



.. code:: python

    my_residue_A23.ampal_parent == my_chain_a




.. parsed-literal::

    True



.. code:: python

    my_atom_A23ca.ampal_parent.ampal_parent




.. parsed-literal::

    <Polypeptide containing 215 Residues. Sequence: DIDTLISNNALW...>



.. code:: python

    my_atom_A23ca.ampal_parent.ampal_parent.ampal_parent




.. parsed-literal::

    <Assembly (3qy1) containing 2 Polypeptides, 449 Ligands>



.. code:: python

    my_residue_A23.ampal_parent.ampal_parent




.. parsed-literal::

    <Assembly (3qy1) containing 2 Polypeptides, 449 Ligands>



.. code:: python

    my_atom_A23ca.ampal_parent.id




.. parsed-literal::

    '23'



.. code:: python

    my_residue_A23.ampal_parent.sequence




.. parsed-literal::

    'DIDTLISNNALWSKMLVEEDPGFFEKLAQAQKPRFLWIGCSDSRVPAERLTGLEPGELFVHRNVANLVIHTDLNCLSVVQYAVDVLEVEHIIICGHSGCGGIKAAVENPELGLINNWLLHIRDIWLKHSSLLGKMPEEQRLDALYELNVMEQVYNLGHSTIMQSAWKRGQNVTIHGWAYSINDGLLRDLDVTATNRETLENGYHKGISALSLKYI'



.. code:: python

    my_chain_a.ampal_parent.sequences




.. parsed-literal::

    ['DIDTLISNNALWSKMLVEEDPGFFEKLAQAQKPRFLWIGCSDSRVPAERLTGLEPGELFVHRNVANLVIHTDLNCLSVVQYAVDVLEVEHIIICGHSGCGGIKAAVENPELGLINNWLLHIRDIWLKHSSLLGKMPEEQRLDALYELNVMEQVYNLGHSTIMQSAWKRGQNVTIHGWAYSINDGLLRDLDVTATNRETLENGYHKGISALSLKYI',
     'KDIDTLISNNALWSKMLVEEDPGFFEKLAQAQKPRFLWIGCSDSRVPAERLTGLEPGELFVHRNVANLVIHTDLNCLSVVQYAVDVLEVEHIIICGHSGCGGIKAAVENPELGLINNWLLHIRDIWLKHSSLLGKMPEEQRLDALYELNVMEQVYNLGHSTIMQSAWKRGQNVTIHGWAYSINDGLLRDLDVTATNRETLENGYHKGISALSLKYI']



7. Ligands
----------

The last AMPAL objects to discuss are ``Ligand`` and ``Ligands``. These
are intended to store non-protein elements from the PDB file. The
ligands can be extracted from the ``Assembly``:

.. code:: python

    my_protein.get_ligands()




.. parsed-literal::

    <Ligands chain containing 449 Ligands>



``Ligands`` is a special ``Polymer`` object, with none of the
protein-specific ``Polypeptide`` functionality. It contains one or more
``Ligand`` objects which you can select in exactly the same way as
selecting ``Residues`` from ``Polypeptides``:

.. code:: python

    my_ligands = my_protein.get_ligands()

.. code:: python

    my_ligands[0]




.. parsed-literal::

    <Ligand containing 1 Atom. Ligand code: ZN>



.. code:: python

    my_ligands['221']




.. parsed-literal::

    <Ligand containing 1 Atom. Ligand code: ZN>



The ``Ligand`` objects are ``Monomer`` objects, without the
protein-specific functionality that is present for ``Residues``. Since
``Ligand`` and ``Residue`` are both examples of ``Monomer``, they have a
lot of the same functionality:

.. code:: python

    my_ligand_zinc = my_ligands[0]

.. code:: python

    my_ligand_zinc.atoms




.. parsed-literal::

    OrderedDict([('ZN',
                  <Zinc Atom (ZN). Coordinates: (-5.817, -20.172, -18.798)>)])



.. code:: python

    my_ligand_zinc['ZN']




.. parsed-literal::

    <Zinc Atom (ZN). Coordinates: (-5.817, -20.172, -18.798)>



This zinc atom is associated with one of the ``Polypeptide`` chains, and
this is reflected in its ``ampal_parent``.

.. code:: python

    my_ligand_zinc.ampal_parent




.. parsed-literal::

    <Polypeptide containing 215 Residues. Sequence: DIDTLISNNALW...>



8. Summary and activities
-------------------------

With these simple methods you can load in a PDB file and select various
different parts of the protein. Please try playing around with the
example code and try to select different parts of the protein.

1. Try loading in a PDB file of your own and select various parts of the
   protein and ligands.
2. Find the other builtin functions either by:

   1. Tabbing the object in Jupyter Notebook
   2. Looking at the documentation
   3. Finding the ``base_ampal`` code in the ISAMBARD folder and looking
      through it (tip: you can do this with the IPython file browser)

In the next section we'll look at how we can perform more complex
selections and more detailed analysis on these objects.
