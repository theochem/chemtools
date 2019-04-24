..
    : ChemTools is a collection of interpretive chemical tools for
    : analyzing outputs of the quantum chemistry calculations.
    :
    : Copyright (C) 2014-2015 The ChemTools Development Team
    :
    : This file is part of ChemTools.
    :
    : ChemTools is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : ChemTools is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --


.. _tutorial_conceptual_global:

Global Descriptive Tools
========================

To refresh your mind on theoretical background, please refer to
:ref:`Scientific Documentation on Global Descriptive Tools <global_tools>`.

To learn about the API documentation, please refer to :mod:`chemtools.toolbox.conceptual`.

To compute global conceptual reactivity descriptors of your molecule of interest, please follow
the steps bellow:


1. Import Global Conceptual DFT Class
-------------------------------------

To use global conceptual DFT class in ChemTools, you should first import it. So,


   .. ipython:: python

      from chemtools import GlobalConceptualDFT


2. Initialize Global Conceptual DFT Class
-----------------------------------------

There are three ways to initialize an instance of ``GlobalConceptualDFT`` class for computing the global
reactivity descriptors for the selected energy model.
In all these, one needs to select the energy model. The supported energy models include:
``'linear'``, ``'quadratic'``, ``'rational'``, ``'exponential'`` and ``'general'``.
The energy model is specified as a sting in place of ``energy_model`` argument in what follows.
To lean more about these models, please refer to
:ref:`Scientific Documentaion on Global Energy Models <global_energy_models>`.

- **Using Output File(s) of Electronic Structure Computations:**

  To learn about supported output file formats in ChemTools, please refer to ``to_be_added``.

  The supported output files can be used directly to initialize an instance of ``GlobalConceptualDFT`` class
  for the selected energy model. These are denoted by ``filenames`` and ``energy_model`` arguments, respectively.

    .. code-block:: python

       model = GlobalConceptualDFT.from_file(filenames, energy_model)

  If ``filenames`` represents one output file, the frontier molecular orbital (FMO) theory approach is used to compute
  global reactivity descriptors. However, if ``filenames`` represents a list of output files, the finite difference (FD)
  approach is taken (the order of output files in the list is not important).


- **Using Dictionary of Energy Values:**

  The energy values corresponding to various number of electrons can also be used to directly initialize an instance
  of ``GlobalConceptualDFT``. These values should be provided as a dictionary (with number of electrons as key and
  corresponding energy as value) alongside the selected energy model. These are denoted by ``dict_values`` and
  ``energy_model`` arguments, respectively, and passed directly to ``GlobalConceptualDFT`` class.

    .. code-block:: python

       model = GlobalConceptualDFT(dict_values, energy_model)

In this tutorial, we will compute formaldehyde's global reactivity descriptors using quadratic energy model
within the FMO approach. So,

  .. ipython:: python

     # obtain path to the formatted checkpoint file for a ub3lyp/aug-cc-pvtz
     # calculation on formaldehyde
     filename = 'ch2o_q+0.fchk'
     # initialize quadratic global conceptual DFT class from one output file
     model = GlobalConceptualDFT.from_file(filename, 'quadratic')

The ``model`` object is an instance of ``GlobalConceptualDFT``, and in this example contains all
quadratic global reactivity descriptors of formaldehyde computed within FMO framework.


3. Get Global Conceptual DFT Reactivity Descriptors
---------------------------------------------------

The ``model`` instance, disregarding of how it has been initialized, contains all the
global reactivity descriptors. To get an overview of the content of this instance,

  .. ipython:: python

     print model

The attributes (i.e. variables) denote available global reactivity descriptors. The methods (i.e. functions)
denote available functions that can be evaluated when given required arguments.
To specifically obtain any of the available attributes listed above, just add its name with a ``.`` after the
``model`` instance. For example,

  .. ipython:: python

     print model.n0                   # reference number of electrons
     print model.n_max                # maximum number of electrons
     print model.softness             # chemical softness
     print model.electronegativity    # electronegativity
     print model.electrophilicity     # electrophilicity
     print model.nucleofugality       # nucleofugality

 Some of these descriptors are stored with two names: a longer (and more clear) name and a shorter (and more convenient) name.
 Both of these will return the same value. For example,

  .. ipython:: python

     print model.ip, model.ionization_potential    # ionization potential
     print model.ea, model.electron_affinity       # electron affinity
     print model.mu, model.chemical_potential      # chemical potential
     print model.eta, model.chemical_hardness      # chemical hardness


4. Compute Higher-Order Reactivity Descriptors
----------------------------------------------

The nth-order ``hyper_hardness`` and ``hyper_softness`` values, for :math:`n\geq2`, can be computed by specifying
the order as an argument for these methods (i.e. functions):

  .. ipython:: python

     print model.hyper_hardness(order=2)
     print model.hyper_hardness(order=3)
     print model.hyper_softness(order=2)
     print model.hyper_softness(order=3)

In this example, as expected for quadratic energy model, the ``hyper_hardness`` values are zero and the
``hyper_softness`` cannot be defined (``None`` is returned for descriptors that are not defined).


5. Compute Energy and its Derivatives
-------------------------------------

The interpolated energy model and its nth-order derivatives (with respect to number of electrons) can be evaluated
for a given number of electrons.
To compute energy, the number of electrons should be provided as an argument to the ``energy`` method
(i.e. functions):

  .. ipython:: python

     print model.energy(15.8)
     print model.energy(model.n0)
     print model.energy(model.n_max)
     print model.energy(model.n_max + 0.1)

To compute energy derivatives, the number of electrons alongside the order of derivative should be provided as an
argument to the ``energy_derivative`` method (i.e. functions):

  .. ipython:: python

     print model.energy_derivative(15.8, order=1)
     print model.energy_derivative(model.n0, order=1)
     print model.energy_derivative(model.n_max, order=1)
     print model.energy_derivative(model.n_max + 0.1, order=1)
     print model.energy_derivative(model.n0, order=2)
     print model.energy_derivative(model.n_max, order=2)
     print model.energy_derivative(model.n0, order=3)

In this example, as expected for quadratic energy model, the 2nd derivative of energy is a constant.
(i.e. its value does not depend on the number of electrons). Also, 3rd and higher derivatives of energy
are all zero.


6. Compute Grand Potential and its Derivatives
----------------------------------------------

The grand potential (corresponding to the interpolated energy model) and its nth-order derivatives (with respect to
number of electrons) can be evaluated for a given number of electrons.
To compute grand potential, the number of electrons should be provided as an argument to ``grand_potential``
method (i.e. functions):

  .. ipython:: python

     print model.grand_potential(15.8)
     print model.grand_potential(model.n0)
     print model.grand_potential(model.n_max)
     print model.grand_potential(model.n_max + 0.1)

To compute grand potential derivatives, the number of electrons alongside the order of derivative should be
provided as an argument to ``grand_potential_derivative`` method (i.e. functions):

  .. ipython:: python

     print model.grand_potential_derivative(15.8, order=1)
     print model.grand_potential_derivative(model.n0, order=1)
     print model.grand_potential_derivative(model.n_max, order=1)
     print model.grand_potential_derivative(model.n_max + 0.1, order=1)
     print model.grand_potential_derivative(model.n0, order=2)
     print model.grand_potential_derivative(model.n_max, order=2)
     print model.grand_potential_derivative(model.n0, order=3)


7. Get Parameters of Energy Model
---------------------------------

The parameters of the interpolated energy model can be obtained by:

  .. ipython:: python

     print model.params

In this example, these denote parameters ``a``, ``b`` and ``c`` of the quadratic energy model.


 .. todo::
    Add more examples?
    Example of general energy model

    .. code-block:: python
       :linenos:

       # define symbols used in the energy expression
       n, a, b, c = sympy.symbols('N, a, b, c')
       # define the energy expression
       expression = a + b * n + c * (n**2)
       # dictionary {N : E(N)}
       energies = {}
       # parametrize energy model
       model = GeneralizedGlobalTool(expression, energies, n)
       # ready to retrieve any global tool
       print model.mu
