#
#
# The main idea is to use HORTON for initializing molecules and getting any basic
# information necessary (e.g. coordinates, chargers, atoms), HORTON supports many file formats (as opposed to PySCF)
# We then use PySCF to perform quantum chemistry calculations
# Different theories (i.e. SCF, DFT etc) could be set up as different classes, although
# another option is to write a main class and then write multiple subclasses (for each theoretical method)
# that would inherit from the main class.
#
#
import os
import numpy as np
from pyscf import gto, scf
from chemtools.wrappers.molecule import Molecule


__all__ = ["RHF"]

#TODO: account for other optional parameters like symmetry.
# Also see if inheritance could be used instade of using 2 separate classes
class RHF(object):
    """ Intialize the molecule via MolSetup class and perform RHF on it"""
    def __init__(self, input_file):
        self.molecule = MolSetup(input_file)
        self.geometry = self.molecule.geometry

    def get_energy(self, basis, charge=0, spin=0, verbose=4, gen_output=False):
        """ Initialize PySCF molecule and perform RHF calculation given certain parameters.
        If opted, the output will be dumped into a .log file. The verbosite can range from 0 (no output)
        to 5 (all of the output, which will a lot of information that may not be useful).
        """
        mol = gto.Mole()
        mol.basis = basis
        mol.atom = self.geometry # atom geometry
        mol.charge = charge
        mol.spin = spin
        if gen_output:
            mol.output = os.path.splitext(self.molecule.input_file)[0] + '.log'
            mol.verbose = verbose
        mol.build()
        rhf = scf.RHF(mol)
        return rhf.kernel()

#TODO: add properties?
class MolSetup(object):
    """ Class initializes from a molecule input_file to be set up for SCF calculations"""

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
        """ Convert an integer to an atom string"""
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
