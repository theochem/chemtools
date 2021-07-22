import sys,os

from rdkit import Chem
from rdkit.Chem import AllChem 
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdMolTransforms 

mol = Chem.MolFromSmiles('Cl[C@H](F)NC\C=C\C')
image = rdMolDraw2D.MolDraw2DCairo(500,400) 

def visualizeChiralCenter(image):
    image.drawOptions().addStereoAnnotation = True

def dihedralAngles(mol,image): 
    # NONFUNCTIONAL  
    """
    for q in mol.GetBonds():
        bonds = (q.GetBeginAtomIdx(),q.GetEndAtomIdx()) 
    AllChem.Compute2DCoords(mol) 
    isomerize = mol.GetConformer()
    
    #TORSION
    rdMolTransforms.GetAngleDeg(isomerize,bonds) 
    # DIHEDRAL 
    """
def imageGenerate(image):
    image.drawOptions().addAtomIndices = True
    image.DrawMolecule(mol) 
    image.FinishDrawing()
    image.WriteDrawingText('test.png')

if len(sys.argv) >= 2:
    if sys.argv[1] == '-c' or sys.argv[1] == '--chirality': 
        visualizeChiralCenter(image) 


dihedralAngles(mol,image) 
imageGenerate(image)
