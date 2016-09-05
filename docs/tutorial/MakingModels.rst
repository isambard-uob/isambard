
Making Models With Isambard
===========================

1. Introduction
---------------

Topologies available
~~~~~~~~~~~~~~~~~~~~

A number of topologies are currently available through Isambard. These
include: alpha-helices, coiled-coil assemblies of different flavours,
deltaprots, pi-helices and poly-proline helices. To make a structure
based on a topology, you need to know the parameters of that topology.
We will demonstrate this with a simple alpha-helical coiled-coil dimer.

Parameters required to make a coiled-coil dimer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: https://cloud.githubusercontent.com/assets/16082702/12421221/08613f3e-beb9-11e5-9bb0-652fd40e33c1.png
   :alt: dimer\_parameters

   dimer\_parameters

To make a coiled-coil dimer, you need to know

1. the number of residues on each helix,
2. the radius of the assembly,
3. the pitch of the assembly and
4. the interface angle at which helices meet (Phi-C:math:`\alpha`)

Reading in a structure and identifying parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We will read a known coiled-coil dimer structure (CC-DI, PDB
`4DZM.pdb <https://drive.google.com/open?id=0B2sZ5She4fA2QzdsUENNa2hmTE0>`__),
into the AMPAL framework, and use tools within Isambard to calculate the
parameters we need. The following code imports isambard so it is
available to the ipython notebook, reads in the PDB file of interest and
converts it into an AMPAL object. [N.B. This is a slightly modified form
of 4DZM where we've converted iodo-phenylalanine to phenylalanine, and
got rid of water molecules to make life a bit easier. We've put it in
the tutorial directory for today, but you can get another copy from the
link above.]

.. code:: python

    import isambard
    
    my_cc = isambard.ampal.convert_pdb_to_ampal("4DZM.pdb")
    
    my_cc




.. parsed-literal::

    <Assembly (4DZM) containing 2 Polypeptides>



This has given us an Assembly object, which contains two Polypeptides in
the AMPAL framework. We can access the sequences and each chain
individually using the ``.sequences`` attribute. This returns a list,
and you can access individual list elements as you would a regular
python list.

.. code:: python

    my_cc.sequences




.. parsed-literal::

    ['GEIAALKQEIAALKKENAALKFEIAALKQGY', 'GEIAALKQEIAALKKENAALKFEIAALKQGY']



.. code:: python

    my_cc_a = my_cc[0]
    my_cc_b = my_cc[1]
    
    my_cc_a




.. parsed-literal::

    <Polypeptide containing 31 Residues. Sequence: GEIAALKQEIAA...>



From the Polypeptide object, which is a special type of Polymer within
the AMPAL framework, we can access the number of residues, and the chain
identifier if needed. You can access the length of the chain either from
the chain object itself, or from the length of the sequence.

.. code:: python

    len(my_cc_a.sequence)




.. parsed-literal::

    31



.. code:: python

    len(my_cc_a)




.. parsed-literal::

    31



.. code:: python

    my_cc_a.id




.. parsed-literal::

    'A'



Before we extract parameters from the structure, we might want to know
how good BUDE thinks the structure is, so that when we build a structure
*de novo* we have something to compare it to. There are two ways of
doing this. The first is via a dedicated method for the Assembly, and
returns the total interaction energy between all Polypeptides in the
Assembly.

.. code:: python

    my_cc.bude_score




.. parsed-literal::

    -580.8352



You can also ask for the average bude score:

.. code:: python

    my_cc.average_bude_score




.. parsed-literal::

    -290.4176



This returns the average interaction energy (i.e. half the value
returned by the ``my_cc.bude_score`` method for this case.). You can get
the same value as ``my_cc.bude_score`` by calling the
``run_bude_additive`` method.

2. Measuring geometric parameters
---------------------------------

In order to measure radius, pitch and Phi-C\ :math:`\alpha` of the
assembly, we need to define a ***reference axis***, that is the line
that runs down the centre of the assembly. For a coiled-coil dimer, it
runs between the two helices; for a barrel, this would be at the centre
of the barrel, and is a list of points in 3D space. We will use this to
calculate the other parameters. The reference axis is defined as a
Primitive chain object populated with PseudoMonomers which represent the
points in space of the axis.

.. code:: python

    reference_axis = isambard.analyse_protein.reference_axis_from_chains(my_cc)

.. code:: python

    reference_axis




.. parsed-literal::

    <Primitive chain containing 31 PseudoMonomers>



Interface angles (Phi-C:math:`\alpha`)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In Isambard, interface angles are accessed via the
``analyse_protein.crick_angles`` method, and are calculated for each
residue in each helix individually.

.. code:: python

    crangles_a = isambard.analyse_protein.crick_angles(my_cc_a,reference_axis)
    crangles_b = isambard.analyse_protein.crick_angles(my_cc_b,reference_axis)

.. code:: python

    crangles_a




.. parsed-literal::

    [165.46311695313284,
     -73.7731235966416,
     26.54672679839078,
     128.47141485160603,
     -127.72847955562587,
     -29.230898305107225,
     75.77896871985133,
     179.13811106898652,
     -78.91543487166474,
     24.619715382288263,
     124.86576008079678,
     -132.1659822842109,
     -33.094770958544274,
     71.5612908445053,
     172.66909579579703,
     -87.21396295102991,
     17.62150247217829,
     122.41786905701451,
     -132.30491551944263,
     -25.580896174877395,
     77.39202166004726,
     178.68641715037646,
     -81.08703803043669,
     21.378411765253517,
     124.46637604063143,
     -137.48525842491375,
     -38.516041308102864,
     65.8051730938853,
     173.8634442139843,
     -69.52238541573473,
     None]



To get an average Phi-C\ :math:`\alpha` value, we take every seventh
value from the list, and take the mean of this. To get the mean we need
to use a module called ``numpy``, which is *Numerical Python* and gives
us access to lots of mathematical functions. In the following cell, we
are going to use it to calculate mean Phi-C\ :math:`\alpha` values.

.. code:: python

    import numpy
    
    phica_a_list = [crangles_a[x] for x in range(0,len(crangles_a),7) if crangles_a[x] is not None]
    phica_b_list = [crangles_b[x] for x in range(0,len(crangles_b),7) if crangles_b[x] is not None]
    
    phica_a = numpy.mean(phica_a_list)
    phica_b = numpy.mean(phica_b_list)
    
    phica_a, phica_b




.. parsed-literal::

    (173.96403703645541, 173.96403703645541)



Radius
^^^^^^

Radius is calculated by looking at successive distances from the centre
of the :math:`\alpha` helix to the reference axis we calculated earlier.
Again, this is calculated for each helix, and we take the mean value to
use for model building.

.. code:: python

    radius_a_list = isambard.analyse_protein.polymer_to_reference_axis_distances(my_cc_a, reference_axis)
    radius_b_list = isambard.analyse_protein.polymer_to_reference_axis_distances(my_cc_b, reference_axis)
    
    radius = numpy.mean(radius_a_list+radius_b_list)
    
    radius




.. parsed-literal::

    5.12974619218142



Pitch
^^^^^

Pitch is calculated on a per-helix basis using alpha angles, which
measures the tilt of a helix in a coiled-coil assembly wtith respect to
the central reference axis that we already calculated. The radius of the
helices is also required. For model building we take a mean value.

.. code:: python

    alpha_a_list = isambard.analyse_protein.alpha_angles(my_cc_a, reference_axis)
    alpha_b_list = isambard.analyse_protein.alpha_angles(my_cc_b, reference_axis)
    
    pitch_a_list = [(2* numpy.pi * radius) / numpy.tan(numpy.deg2rad(x)) for x in alpha_a_list if x is not None]
    pitch_b_list = [(2* numpy.pi * radius) / numpy.tan(numpy.deg2rad(x)) for x in alpha_b_list if x is not None]
    
    pitch = numpy.mean(pitch_a_list + pitch_b_list)
    
    pitch




.. parsed-literal::

    220.28472172509427



3. Building a model
-------------------

We now know all the parameters we need to rebuild this structure, and it
can be done by generating a topology object from the CoiledCoil class.
We will do this by 'subclassing', i.e. inheriting all the things we need
from the general CoiledCoil class, and adding specific parameters for
our case, into a class called ``SimpleDimer2Phi``.

.. code:: python

    class SimpleDimer2Phi(isambard.topology.CoiledCoil):
        def __init__(self, aa, r, p, phica1,phica2,n=2):
            super().__init__(n, auto_build=False)
            self.aas = [aa]*n
            self.major_radii = [r]*n
            self.major_pitches = [p]*n
            self.phi_c_alphas = [phica1,phica2]
            self.orientations = [1]*n 
            self.build()

We can now make a model by calling the ``SimpleDimer2Phi`` class and
filling in the parameters obtained earlier in the workbook. Enter these
parameters in the cell below before pressing Ctrl+Enter to run it.

.. code:: python

    my_model = SimpleDimer2Phi(len(my_cc.sequences[0]), radius, pitch, phica_a, phica_b)
    my_model.build()

Viewing models
~~~~~~~~~~~~~~

We have built a visualizer for AMPAL objects into Isambard. It is very
basic at the moment (written by Chris in a very short space of time this
week!) First of all, you have to create a view using an AMPAL object,
then call the ``.view()`` method on it to visualize it in your browser.
Finally, the ``.control_panel()`` method allows you to change the
colors, etc. We expect to make changes and improvements to this over the
coming weeks. Execute the three cells below to view your model.

.. code:: python

    my_view = isambard.add_ons.AMPALViewer(my_model)



.. parsed-literal::

    <IPython.core.display.Javascript object>


.. code:: python

    my_view.view()

.. code:: python

    my_view.control_panel()




.. parsed-literal::

    <function ipywidgets.widgets.interaction.interact.<locals>.<lambda>>



Now write out this model to a PDB file for later reference.

.. code:: python

    with open('my_model.pdb','w') as outf:
        outf.write(my_model.pdb)

.. figure:: https://cloud.githubusercontent.com/assets/16082702/12421222/0861c33c-beb9-11e5-919d-67168da10135.png
   :alt: my\_model\_bb

   my\_model\_bb

Hopefully your model looks something like the picture above. It has made
a poly-alanine version of a coiled-coil helical dimer. Now we can model
the sidechains using SCWRL, output the PDB file for future reference,
and view the structure in the browser. You should be able to see the
sidechains.

.. code:: python

    sequences = my_cc.sequences
    my_model.pack_new_sequences(sequences)
    with open('my_model_sc.pdb','w') as outf:
        outf.write(my_model.pdb)
    
    my_view = isambard.add_ons.AMPALViewer(my_model)
    my_view.view()



.. parsed-literal::

    <IPython.core.display.Javascript object>


.. code:: python

    my_view.control_panel()

.. figure:: https://cloud.githubusercontent.com/assets/16082702/12421307/85a3b2f6-beb9-11e5-848a-cd43b081f2be.png
   :alt: my\_model\_sc

   my\_model\_sc

Your model should look something like the picture above. Now we can
score this model using BUDE, and also we can calculate an RMSD from this
model to the original structure.

.. code:: python

    my_model.bude_score




.. parsed-literal::

    -600.249



.. code:: python

    rmsds,overlay = isambard.external_programs.run_profit_output_pdbstr(my_cc.pdb, my_model.pdb)

``rmsds`` is a list of C-alpha, backbone and all-atom RMSDs calculated
by ProFit.

.. code:: python

    rmsds




.. parsed-literal::

    [0.913, 0.899, 2.173]



Write out the overlaid structure for future reference.

.. code:: python

    with open('my_model_overlaid.pdb','w') as outf:
        outf.write(overlay)

.. figure:: https://cloud.githubusercontent.com/assets/16082702/12421220/0860a862-beb9-11e5-8682-13974df51711.png
   :alt: my\_cc\_my\_model\_overlay

   my\_cc\_my\_model\_overlay

The overlay is a simple PDB string returned from ProFit. We can convert
this back to an AMPAL object, and combine it with our original structure
AMPAL object, and view both overlaid. [NB This is a very experimental
feature, only just implemented, so it may break!!]

.. code:: python

    overlay_ampal = isambard.ampal.convert_pdb_to_ampal(overlay,path=False) ## convert overlay to AMPAL object
    both = my_cc + overlay_ampal ## combine AMPAL objects
    both_view = isambard.add_ons.AMPALViewer(both) ## set up viewer



.. parsed-literal::

    <IPython.core.display.Javascript object>


.. code:: python

    both_view.view()

.. code:: python

    both_view.control_panel()

If you have a look closely at the asparagine pair at the middle of the
dimeric structure (shown below), you'll see that the rotamers that SCWRL
has picked for asparagine (right-hand side, dark blue) are more sensible
than those in the original structure (left-hand side, cyan). The
original 4DZM structure has multiple occupancies for the asparagine
residues, and for simplicity we picked the first pair. The structure
that SCWRL has returned has a hydrogen bond between the asparagine
residues.

.. figure:: https://cloud.githubusercontent.com/assets/16082702/12421660/29935b86-bebb-11e5-9746-9b9f6588113f.png
   :alt: asn\_hbonds

   asn\_hbonds

Tweaking the parameters
~~~~~~~~~~~~~~~~~~~~~~~

Try rebuilding the model again, this time changing the parameters as you
like. What happens if you change the phi-C\ :math:`\alpha` values of one
helix (or both)? Try varying the radius. Each time, score your model
using BUDE and make a note of the RMSD to the original structure, and
have a look at the models you produce on PyMOL.

4. Going antiparallel
---------------------

We are going to remake the coiled-coil dimer as an antiparallel
structure. To do this, we need to modify the
``isambard.topology.SimpleDimer2Phi`` class, two extra parameters to
specify firstly that the orientation of the second helix is
antiparallel, and secondly, the z-shift, how far the helices can slide
past one another. We'll make this using another subclass of CoiledCoil:

.. code:: python

    class SimpleDimer2PhiAP(isambard.topology.CoiledCoil):
        def __init__(self, aa, r, p, phica1,phica2,zshift,n=2):
            super().__init__(n=2, auto_build=False)
            self.aas = [aa]*n
            self.major_radii = [r]*n
            self.major_pitches = [p]*n
            self.phi_c_alphas = [phica1,phica2]
            self.z_shifts = [0,zshift]
            self.orientations = [1,-1] 
            self.build()

.. code:: python

    my_ap_cc = SimpleDimer2PhiAP(len(my_cc.sequences[0]),radius,pitch,phica_a,phica_b,0.0)
    my_ap_cc.build()
    my_ap_cc.pack_new_sequences(my_cc.sequences)

Write the model out to PDB for future reference and view in the viewer.

.. code:: python

    with open('my_ap_cc_model.pdb','w') as outf:
        
        outf.write(my_ap_cc.pdb)

.. code:: python

    my_ap_view = isambard.add_ons.AMPALViewer(my_ap_cc)
    my_ap_view.view()



.. parsed-literal::

    <IPython.core.display.Javascript object>


.. code:: python

    my_ap_view.control_panel()




.. parsed-literal::

    <function ipywidgets.widgets.interaction.interact.<locals>.<lambda>>



Score the model using BUDE:

.. code:: python

    my_ap_cc.bude_score




.. parsed-literal::

    -101.5608



5. Making it better - an introduction to optimisation
-----------------------------------------------------

It's not a good score :-) (obviously?) We could improve this score by
'minimizing' the structure, i.e. trying to find new parameter for the
sequence that improve the BUDE score. This is, in real life, a very bad
idea, because we know that this sequence comes from a crystal structure
of a parallel coiled-coil dimer, and it is highly unlikely ever to go
antiparallel, but we are going to pretend that we can make it so, as a
good exercise in learning how an optimiser works, and also to give you
the heads up that although you can get an improved BUDE score through
this process, if you make an unrealistic assumption at the beginning,
BUDE and the optimiser will not help you out.

We'll carry out this optimisation using an algorithm called
**differential evolution**. This is a metaheuristic, which means it
makes no assumptions about the optimisation (this means it knows nothing
about your parameters), and is able to search large spaces of candidate
solutions (i.e. a large parameter space). It optimizes the parameters by
keeping a population of candidate solutions (sets of model parameters)
and creates new solutions (sets of parameters) by combining existing
ones and then keeping whichever candidate has the best score.

The optimiser is actually quite easy to set up. You need to give it an
isambard topology (in our case ``SimpleDimer2Phi``) and an output path
(which can be your home directory for the moment)

.. code:: python

    optimiser = isambard.optimisation.OptDE(SimpleDimer2PhiAP,output_path='/Users/chgjb/')

The next step is to give the optimiser the sequences to model, a set of
parameter means and variances within which to sample. Again, you can see
what the optimiser is expecting by using SHIFT+TAB inside the brackets
of ``optimiser.parameters``.

The order of the parameters inside the value means and value ranges must
match the order required by ``SimpleDimer2Phi``, which is
``radius, pitch, phica1, phica2``. We will allow the radius to vary by
up to 1 :math:`\unicode{x212B}` either side of 5
:math:`\unicode{x212B}`, the pitch by 100 :math:`^\circ` and the
Phi-C\ :math:`\alpha`\ s by a full circle. We'll keep the z-shift to
zero for this instance.

.. code:: python

    optimiser.parameters(my_ap_cc.sequences, [5.0,180,180,180], [1.0,100,180,180], 
                        [len(my_ap_cc.sequences[0]), 'var0','var1','var2','var3', 0])

Now, run it! You need to specify the size of each generation and the
number of generations. We will use 20 models per generation and 30
generations.

.. code:: python

    %matplotlib inline
    optimiser.run_de(20,30,1,plot=True,log=True)


.. parsed-literal::

    Starting minimisation (2016-07-25 11:30:53)
    gen	evals	avg     	std    	min     	max    
    0  	20   	-70.3353	128.851	-163.287	448.254
    1  	20   	-114.398	35.6882	-163.287	-48.2052
    2  	20   	-135.619	26.2276	-174.117	-73.2537
    3  	20   	-142.132	24.928 	-174.117	-73.2537
    4  	20   	-149.806	26.282 	-204.855	-88.3584
    5  	20   	-156.847	27.5945	-204.855	-88.3584
    6  	20   	-168.784	24.0334	-212.613	-111.862
    7  	20   	-169.648	24.1884	-212.613	-111.862
    8  	20   	-178.925	22.9477	-235.511	-134.796
    9  	20   	-183.644	20.0573	-235.511	-141.32 
    10 	20   	-187.08 	18.3236	-235.511	-151.413
    11 	20   	-192.341	18.1955	-235.511	-167.127
    12 	20   	-195.731	17.0495	-235.511	-167.127
    13 	20   	-200.526	14.8479	-235.511	-178.863
    14 	20   	-204.964	15.9058	-238.203	-178.863
    15 	20   	-206.692	15.4156	-238.203	-178.863
    16 	20   	-207.442	15.151 	-238.203	-178.863
    17 	20   	-207.442	15.151 	-238.203	-178.863
    18 	20   	-210.471	14.6684	-238.203	-178.863
    19 	20   	-211.328	14.4909	-238.203	-178.863
    20 	20   	-212.308	12.8293	-238.203	-186.033
    21 	20   	-213.954	12.6931	-238.203	-186.033
    22 	20   	-217.335	14.6989	-249.95 	-186.033
    23 	20   	-219.56 	12.7335	-249.95 	-201.021
    24 	20   	-222.445	12.524 	-249.95 	-201.021
    25 	20   	-225.887	10.6563	-249.95 	-207.1  
    26 	20   	-230.533	11.6984	-251.719	-210.016
    27 	20   	-232.563	11.8331	-251.719	-210.016
    28 	20   	-233.251	12.9418	-260.737	-210.016
    29 	20   	-235.425	14.4558	-260.737	-210.016
    End of minimisation (2016-07-25 11:35:26)
    Run ID is DE_model
    Minimization time = 0:04:32.875241
    Evaluated 620 models in total
    Best score is -260.7365778820732
    Best parameters are [31, 5.257943612728227, 82.14285642448519, 120.1487038846328, 285.68371939341034, 0]
    Warning! Parameter 2 is at or near minimum allowed value
    ----Minimisation plot:



.. image:: MakingModels_files/MakingModels_80_1.png


The output above shows you the progress of the optimisation process,
which is reproduced in a log file in the location we requested the
output to be. Additionally the best model pdb file will be written to
the same location. Open this up and have a look at the optimised model.
What has the optimiser done to the model?

There are different optimization methods available which will be
explained in further tutorials.

6. Making bigger things
-----------------------

Now you know how to subclass the CoiledCoil AMPAL topology object, see
if you can make an alpha-helical barrel by increasing the value of
``n``. Don't worry about what sequence to use initially, make some
models and get comfortable with the parameters - for example try a few
different radii, pitch and phi-c\ :math:`\alpha` values. As an aide,
here are some parameters for the coiled-coil barrels in Science paper.
Score the models you make using BUDE.

+-----------+---------+----------+---------+------------------------+-------------------+
| Name      | Oligo   | Radius   | Pitch   | PhiC\ :math:`\alpha`   | Sequence (g->f)   |
+===========+=========+==========+=========+========================+===================+
| CC-Pent   | 5       | 8.62     | 182.8   | 14.44                  | KIEQILQ           |
+-----------+---------+----------+---------+------------------------+-------------------+
| CC-Hex2   | 6       | 9.48     | 131.7   | 18.22                  | EIAKSLK           |
+-----------+---------+----------+---------+------------------------+-------------------+
| CC-Hex3   | 6       | 9.72     | 162.2   | 13.07                  | EIAQSIK           |
+-----------+---------+----------+---------+------------------------+-------------------+
| CC-Hept   | 7       | 9.80     | 328.6   | 15.1                   | EIAQALK           |
+-----------+---------+----------+---------+------------------------+-------------------+
