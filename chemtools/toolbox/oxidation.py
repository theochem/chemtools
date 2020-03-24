# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2019 The ChemTools Development Team
#
# This file is part of ChemTools.
#
# ChemTools is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# ChemTools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
"""Module for Oxidation State."""


from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.grid import MolecularGrid
from chemtools.wrappers.part import DensPart

from horton import ProAtomDB, BeckeMolGrid

import numpy as np
import scipy.linalg as la
import glob as glob
import itertools

from operator import itemgetter




class EOS(object):
    def __init__(self, molecule, part, grid):
        self.molecule = molecule
        self.part = part
        self.grid = grid

    @classmethod
    def from_molecule(cls, molecule, part, grid):
        """Initialize class from `Molecule` object.

        Parameters
        ----------
        molecule : `Molecule`
            Instance of `Molecular` class.
        part : `DensPart`
            Instance of `DensPart` class.
        grid : `MolecularGrid`
            Molecular numerical integration grid.

        """
        return cls(molecule, part, grid)

    @classmethod
    def from_file(cls, fname, part, grid):
        """Initialize class using wave-function file.

        Parameters
        ----------
        fname : str
            A string representing the path to a molecule's fname.
        part : `DensPart`
            Instance of `DensPart` class.
        grid : `MolecularGrid`
            Molecular numerical integration grid.

        """
        molecule = Molecule.from_file(fname)
        return cls.from_molecule(molecule, part, grid)

    def compute_fragment_overlap(self, fragments=None, spin='ab'):
        # compute MO overlap matrix for fragments
    

      nalpha=int(self.molecule.mo.nelectrons[0])
      nbeta= int(self.molecule.mo.nelectrons[1])

      orbitals_a=[self.molecule.compute_molecular_orbital(self.grid.points, "a", index=i).ravel() for i in range(1,nalpha+1)]
      orbitals_b=[self.molecule.compute_molecular_orbital(self.grid.points, "b", index=i).ravel() for i in range(1,nbeta+1)]


   
      s_ov_a=[self.grid.integrate(i[0]*i[1]) for i in itertools.product(orbitals_a, repeat=2)]
      s_ov_b=[self.grid.integrate(i[0]*i[1]) for i in itertools.product(orbitals_b, repeat=2)] 
 
      return orbitals_a, orbitals_b

    def compute_oxidation_state(self, fragments=None):

      nalpha=int(self.molecule.mo.nelectrons[0])
      nbeta= int(self.molecule.mo.nelectrons[1])

      orbitals_a, orbitals_b = self.compute_fragment_overlap()


     #defining fragments/atoms
      if fragments is None:
        frags =  [[item] for item in self.molecule.numbers]
      else:
        frags = fragments


#      #generating qij for each atom/fragment
      qij_alpha=np.zeros((len(frags),nalpha,nalpha))
      qij_beta=np.zeros((len(frags), nbeta, nbeta))

      for x in range(len(frags)):
        qij_a=[self.part.condense_to_fragments(i[0]*i[1], fragments, w_power=2)[x] for i in itertools.product(orbitals_a, repeat=2)]
        qij_b=[self.part.condense_to_fragments(i[0]*i[1], fragments, w_power=2)[x] for i in itertools.product(orbitals_b, repeat=2)]
        qij_a=np.reshape(qij_a, (nalpha,nalpha))
        qij_b=np.reshape(qij_a, (nbeta, nbeta))
        qij_alpha[x]=qij_a
        qij_beta[x]=qij_b

     #qij_a diagonalization
      u_a, s_a, vt_a = np.linalg.svd(qij_alpha)
      u_b, s_b, vt_b = np.linalg.svd(qij_beta)

      
     #compute oxidation state for fragments

      occupations_alpha=[]
      occupations_beta=[]


      for a in range(len(frags)):
        for i in range(nalpha):
          occupations_alpha.append((s_a[a][i], a))

      for a in range(len(frags)):
        for i in range(nbeta):     
          occupations_beta.append((s_b[a][i], a))



      sorted_alpha = sorted(occupations_alpha, key=itemgetter(0), reverse=True)
      sorted_beta = sorted(occupations_beta, key=itemgetter(0), reverse=True)


      s_a_occ=[[] for i in range(len(frags))]
      s_b_occ=[[] for i in range (len(frags))]

      s_a_uncc=[[] for i in range(len(frags))]
      s_b_uncc=[[] for i in range (len(frags))]


      for index, e in enumerate(sorted_alpha):
        if index < nalpha:
          s_a_occ[e[1]].append(e)
          lo_a = index
        else:
          s_a_uncc[e[1]].append(e)


      for index, e in enumerate(sorted_beta):
        if index < nbeta:
          s_b_occ[e[1]].append(e)
          lo_b = index
        else:
          s_b_uncc[e[1]].append(e)


     #Reliability index
      if len(frags) == 1:
        r_alpha = 100.000
      else:      
        fu_a = lo_a + 1
        while True:
          if sorted_alpha[lo_a][1] == sorted_alpha[fu_a][1]:
            fu_a = fu_a + 1
          else:
            break

        fu_b = lo_b + 1
        while True:
          if sorted_beta[lo_b][1] == sorted_beta[fu_b][1]:
            fu_b = fu_b + 1
          else:
            break

        r_alpha = 100 * (sorted_alpha[lo_a][0] - sorted_alpha[fu_a][0] + 0.5)
        r_beta = 100 * (sorted_beta[lo_b][0] - sorted_beta[fu_b][0] + 0.5)


     
      for a in range(len(frags)):
        print 'Fragment', a, 'net occupations'
        print 'alpha occupied ', s_a_occ[a]
        print 'alpha unoccupied ', s_a_uncc[a]
        print
        print 'beta occupied ', s_b_occ[a]
        print 'beta unocupied ', s_b_uncc[a] 
        print 
 
        


#     print(sorted_alpha[lo_a], sorted_alpha[fu_a], r_alpha)
#     print(sorted_beta[lo_b], sorted_beta[fu_b], r_beta)


      print 'Reliability index R(%) =', r_alpha
     
      print'Fragment', '    ' , 'oxidation state'
      for a in range(len(frags)):
        occ = len(s_a_occ[a]) + len(s_b_occ[a])
        z = 0
        if fragments is None:
          z = self.molecule.numbers[a]
        else:
          for elem in frags[a]:
            z = z + self.molecule.numbers[elem]

        oxidation = z - occ
        print a,'    ' , oxidation
      
      

   


    def compute_effective_orbital(self):
      pass








