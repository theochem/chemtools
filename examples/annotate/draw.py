import sys,os

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


mol = Chem.MolFromSmiles('Cl[C@H](F)NC\C=C\C')
image = rdMolDraw2D.MolDraw2DCairo(500,400) 

def visualizeChiralCenter(image):
    image.drawOptions().addStereoAnnotation = True

#def visFormalCharge(mol,image): 

    





def imageGenerate(image):
    image.drawOptions().addAtomIndices = True
    image.DrawMolecule(mol) 
    image.FinishDrawing()
    image.WriteDrawingText('test.png')

if len(sys.argv) >= 2:
    if sys.argv[1] == '-c' or sys.argv[1] == '--chirality': 
        visualizeChiralCenter(image) 

imageGenerate(image)
