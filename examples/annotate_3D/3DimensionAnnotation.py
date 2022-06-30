import iodata
from iodata import load_one, dump_one 

def Annotate3D(inFile):
    
    """
    This script allows one to take an input file, generates an intermediate PDB file, then writes data values to the B Factor Column of the intermediate to produce the final output PDB
    
    To View Annotations In ChimeraX 
    -------------------------------
        Download Chimerax
        -----------------
            Nightly and Stable Releases of ChimeraX can be found below, courtesy of University of California San Francisco
            
            https://www.cgl.ucsf.edu/chimerax/download.html
            
        Visualizing Annotation
        ----------------------
            1. Open ChimeraX
            2. Either enter into the ChimeraX command line `open outFile.pdb` or use File > Open from the uppermost toolbar 
            3. Once the PDB has been loaded, run from ChimeraX command line `label #1 atoms attribute bfactor` or use Actions > Label > Atoms > Other > Bfactor, from the afforementioned toolbar 
    
    Yields 
    ------
    This function produces a PDB File (outFile), which can then be viewed in ChimeraX
    
    Parameters
    ----------
    
        inFile : str
            This represents the input file, which should be a Gaussian CheckPoint File, also known as *.fchk
          
        workingBFactor : IOData call
            This represents the data we wish to write to PDB B Factor Column and annotate.
            While this example uses the IOData Call atcharges['esp'][:] to input electrostatic potential charges, it can be other values
            IOData calls can be found here:
                
                https://iodata.readthedocs.io/en/latest/index.html
                
        interFile : str
            This is the intermediate PDB File in which B Factor Values are written into, it is named 'intermediate.pdb'
        
        outFile : str
            This is the final PDB File in which our B Factor Values are written into, originating from interFile
            This is the file we open in ChimeraX
            
        title : str
            Title of the molecule
        
    --------------            
    End Parameters 
     """
    
    inFile = 'dichloropyridine26_q+0.fchk'
    interFile = 'intermediate.pdb'
    outFile = f"{inFile[:-5]}.pdb" 
    title = inFile[:-5]

    print("Filenames: ")
    print(inFile)
    print(interFile)
    print(outFile)
    print(" ")
    print("Name of Molecule:")
    print(title)
    print(" ")

    loadInput = load_one(inFile)
    
    # EDIT workingBFactor TO CHANGE ANNOTATION VALUES
    workingBFactor = loadInput.atcharges['esp'][:] 
    
    
    print("Our Working B Factor to replace within IOData Dump : ")
    print(workingBFactor)
    print(f"\nWill be written to {outFile} and dumped")
    dump_one(loadInput,interFile)
    
    # END FUNCTION 
    
# EXAMPLE OF Annotate3D(inFile) using 2,6-Dichloropyridine
inFile = 'dichloropyridine26_q+0.fchk' 

Annotate3D(inFile) 


