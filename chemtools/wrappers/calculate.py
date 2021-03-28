#
#
# The main idea is to use HORTON for initializing molecules and getting any basic
# information necessary (e.g. coordinates, chargers, atoms), HOTRON supports many input_file formats (as opposed to PySCF)
# We then use PySCF to perform quantum chemistry calculations
# Different theories (i.e. SCF, DFT etc) could be set up as different classes, although
# some could overlap
#
#
import os
import numpy as np
from pyscf import gto, scf
from chemtools.wrappers.molecule import Molecule


__all__ = ["RHF"]

#TODO: account for other optional parameters like symmetry. Also see if inheritance could be used instade of using 2 separate classes
class RHF(object):
    """ Intialize the molecule via MolSetup class and perform RHF on it"""
    def __init__(self, input_file):
        self.molecule = MolSetup(input_file)
        self.geometry = self.molecule.geometry

    def get_energy(self, basis, charge=0, spin=0, verbose=4, gen_output=False):
        mol = gto.Mole()
        mol.basis = basis
        mol.atom = self.geometry # atom geometry
        mol.charge = charge
        mol.spin = spin
        if gen_output:
            mol.output = os.path.splitext(self.molecule.input_file)[0] + '.log'
            mol.verbose = verbose # how much you want the output to contain (0 to 5)
        mol.build()
        rhf = scf.RHF(mol)
        return rhf.kernel()

#TODO: add properties?
class MolSetup(object):
    """ Class initializes from a molecule input_file to be set up for SCF calculations """

    def __init__(self, input_file):
        self.input_file = input_file
        # Initialize the molecule from a input_file
        self.mol = Molecule.from_file(input_file)
        self.coordinates = self.mol.coordinates
        self.charges = self.mol.numbers
        self.atoms = self.get_atom_list(self.charges)
        self.geometry = list(zip(self.atoms, self.coordinates))

    @staticmethod
    def charge_to_atom(num):
        """ Convert an integer to an atom string.
            If integer does not correspond to any atom, return None"""
        atoms = {1: 'H', 6: 'C', 7: 'N', 8: 'O'}
        if num in atoms:
            return atoms[num]
        else:
            raise ValueError('{0} does not correspond to any atom in the atom dict'.format(num))

    @classmethod
    def get_atom_list(cls, charge_list):
        """ Get atom list from a charge list """
        transform = np.vectorize(cls.charge_to_atom)
        return transform(charge_list)
