########
Overview
########

ISAMBARD (Intelligent System for Analysis, Model Building And Rational Design), is a Python-based framework for 
structural analysis and rational design of biomolecules. It is developed and maintained by members of the
`Woolfson group <http://www.chm.bris.ac.uk/org/woolfson/index.html>`_, University of Bristol.

The ISAMBARD API provides a suite of tools for biomolecular stucture analysis, protein design, model building and 
evaluation.

Goals
#####

ISAMBARD is intended to provide
 
* tools for the study and analysis of biomolecular structure(s).
* a standard representation of these biomolecules, allowing ease of conversion between natural and designed structures.
* a rapid development environment for *de novo* deisgn and bioinformatics studies.
* modifiable pipelines to assist the deisgn of feasible *de novo* biomolecules in tractable computational time.   

Analysis, model building and design
###################################

All biomolecules in ISAMBARD, including natural structures and designed models, are represented using the AMPAL 
(Atom, Monomer, Polymer, Assembly, Ligand) framework. 
This is a formal representation of biomolecules in a hierarchical structure of lightweight Python objects. 

A range of geometric and bioinformatic tools are attached to this framework, facilitating the interrogation of natural 
and designed structures.

There are a number of specifications provided to aid the design of entirely new proteins using simple sets of 
parametric equations.

The Bristol University Force Field (BUFF) is a stand-alone implementation of the all-atom force field
`BUDE <http://comjnl.oxfordjournals.org/content/early/2011/09/12/comjnl.bxr091>`_ is provided for assessing 
structural quality.  

Further, we provide a set of optimisation strategies for quickly finding the 'best' solutions for a given 
sequence-structure problem.

A more in-depth look at ISAMBARD's functionality is provided through the tutorials and the API reference.

Python, Cython and C
####################
The majority of ISAMBARD's code is written in Python, which allows for clear, readable and expandable code. 
Parts of the code-base, including BUFF and the geomtric tools, are written in C++ and Python. 
Communication between these layers is achieved using `Cython <http://cython.org/>`_. 
This allows us to achieve a speed-up in the code implementation, whilst retaining the interactivity and readability of 
Python.

Contribute
##########

ISAMBARD users are encouraged to contribute to the code, improving and expanding existing functionality. 
Please read our Developer Guide for more information.


License
#######

ISAMBARD uses The MIT License (MIT)

Copyright (c) 2016 Woolfson Group, University of Bristol.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

