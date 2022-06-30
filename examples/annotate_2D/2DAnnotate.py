from iodata import load_one, dump_one
    
import glob
import os
    
import numpy as np
import pandas as pd
    
from chemml.chem import Molecule

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import Atom
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import MolToImage
from rdkit.Chem.Draw import DrawingOptions
    
def annotation2D(inputName):
    """
    This function is designed to consolidate what is necessary to use RDKit to visualize molecular data onto a 2D representation of a molecule
    
    Conda Environment
    -----------------
        This function requires a conda environment containing
            1. Pandas
            2. Numpy
            3. Glob
            4. OS
            5. Chemml
            6. RDKit 
    Yields
    ------
        This function produces a PNG image of our input molecule (inFile), annotated with our assigned values (annotateValue)
    Parameters
    ----------
        inputName, str
            The file name of the Gaussian Checkpoint *.fchk file that will be processed
        
        annotateValue : IOData Call
            This represents the value that we want to annotate our molecule with.
            While annotateValue does not necessarily need to be an IOData Call, this example uses IOData's atcharges['esp'][:]
            IOData Examples can be found here
            
                https://iodata.readthedocs.io/en/latest/index.html
            
        xyzName, str
            Derived from inputName
            This represents the name of the intermediate XYZ File in which the SMILES format is extracted from
            
        pdbName, str
            Derived from inputName
            This represents the name of the intermediate PDB File in which the atomic indices are preserved for annotation
      
        molName, str
            Derived from xyzName
            This represents the name of the molecule
        
        imageName, str
            Derived from molName
            This represents the name of the annotated image. 
            
        xyz_files, glob
            This is meant to find all relevant XYZ Files
        
        mol and molecule_list, chemml object
            This is the loaded XYZ file, and will be tabulated in order to find the SMILES String
       
        table, Pandas DataFrame
            Tabulated format allowing us to extract the SMILES String from the chemml object
            
        smiles, str
            This represents the extracted SMILES string from our XYZ intermediate file
            
        chem, RDKit mol object 
            This represents our intermediate PDB loaded into RDKit 
        
        template, RDKit mol object
            This represents our SMILES string (originally from XYZ intermediate) loaded into RDKit 
       
        tempMol, corrected RDKit mol object
            This represents our PDB Corrected using the XYZ SMILES String
            This is accomplished by using AllChem.AssignBondOrdersFromTemplate
     
        drawMolecule, RDKit mol object visualization
            This represents our RDKit object being visualized
            MolDraw2DCairo is used to output to PNG
            
    -------------- 
    End Parameters
    """
    
    loadInput = load_one(inputName) 

    xyzName = f"{inputName[:-4]}xyz"    # EDIT THIS TO CHANGE INTERMEDIATE XYZ FILE
    pdbName = f"{inputName[:-4]}pdb"    # EDIT THIS TO CHANGE INTERMEDIATE PDB FILE 
    molName = xyzName[:-4]
    imageName = f"{molName}.png" 

    print("Molecule Name : ")
    print(molName)
    print(xyzName)

    dump_one(loadInput,pdbName)     # DUMP XYZ
    dump_one(loadInput, xyzName)    # DUMP PDB
   
    # EDIT annotateValue TO CHANGE ANNOTATIONS 
    annotateValue = loadInput.atcharges['esp'][:] #EDITABLE 

    # Extract SMILES String from XYZ Intermediate
    xyz_files = glob.glob(xyzName)
    mol = Molecule('dichloropyridine26_q+0.xyz', 'xyz')

    molecule_list = [Molecule(xyzName, input_type='xyz') for file in xyz_files]
    for molecule in molecule_list:
            molecule.to_smiles()
    print(molecule_list)
    table = pd.DataFrame(data={'Name':xyzName,  
                         'SMILES':[molecule.smiles for molecule in molecule_list]}, 
                   columns=['Name','SMILES'])

    print(" ")
    print(table)
    print(" ")
        
    smiles = table.iloc[0,1]
    print("Our Molecule in SMILES format: ")
    print(smiles)
    print(" ")
    
    # Force PDB Intermediate Bond Orders to SMILES String 
    chem = Chem.MolFromPDBFile(pdbName, sanitize=True)
    template = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(chem)
    tempMol = AllChem.AssignBondOrdersFromTemplate(template,chem)

    # Visualize Molecule
    
    drawMolecule = rdMolDraw2D.MolDraw2DCairo(300, 300)
    rdMolDraw2D.PrepareAndDrawMolecule(drawMolecule, chem)
    drawMolecule.drawOptions().addStereoAnnotation = True
    
    # Annotate Molecule
    for i, atom in enumerate(tempMol.GetAtoms()):
        atom.SetProp("atomNote", f"{np.round(float(annotateValue[i]),3)}") 
    
    # Output Visualization to File
    Chem.Draw.MolToFile(tempMol,imageName) 
    print(" ")
    print(f"The molecule {molName} has been annotated and output as file {imageName}")
    print(" ")
    #END FUNCTION 


# Example using 2,6-Dichloropyridine
inputName = 'dichloropyridine26_q+0.fchk'
annotation2D(inputName)



