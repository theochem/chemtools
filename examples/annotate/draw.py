import sys,os
from rdkit import Chem
from rdkit.Chem import AllChem 
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdMolTransforms 
from rdkit.Chem import PeriodicTable

mol = Chem.MolFromSmiles('Cl[C@H](F)NC\C=C\C')
image = rdMolDraw2D.MolDraw2DCairo(500,400) 
outName = "out.png"


def visualizeChiralCenter(image):
    image.drawOptions().addStereoAnnotation = True

def atomicNumber(mol,image): 
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetProp("molAtomMapNumber", str(atom.GetAtomicNum()+1))

def imageGenerate(image):
  #  image.drawOptions().addAtomIndices = True
    image.DrawMolecule(mol) 
    image.FinishDrawing()
    image.WriteDrawingText(outName)


# THIS CODE ALLOWS MULTIPLE FLAGS 
i = 1
if len(sys.argv) >= 2:
    while True:
        if sys.argv[i] == '-o' or sys.argv[i] == '--output':
           outName = sys.argv[i+1]  
        if sys.argv[i] == '-c' or sys.argv[i] == '--chirality': 
            visualizeChiralCenter(image) 
        if sys.argv[i] == '-n' or sys.argv[i] == '--atom-number': 
            atomicNumber(mol,image) 
        if sys.argv[i] == '-h' or sys.argv[i] == '--help':
            print(" This script visualizes atomic properties")
            print(" Input currently must be in SMILES") 
            print(" ")
            print("-h   --help      obtain this message")
            print("-c   --chirality visualize E-Z and Stereoisomers") 
            print("-n   --atom-number visualize atomic numbers") 
            quit() 

        i = i + 1
        if len(sys.argv) <= i: 
            imageGenerate(image)
            quit() 
