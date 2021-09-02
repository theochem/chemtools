
from chimerax.core.commands import run
from chimerax.core import errors
from chimerax import io

def print_chimerax_isosurfaces(session,outFile,outSuffix, isoFile, colorFile, isoSurf, material='shiny', scalemin='compute',scalemax='compute',colorscheme='rainbow',representation='surface', lighting='full', shadows='False', negative=False): 
    """
     This script is designed to visualize the electrostatic potential isosurfaces generated from Chemtools (original file is a Gaussian Checkpoint, then processed to cube files *_esp.cube and *_rho.cube, representing the surface and map files respectively)
     
    Execution
    ----------
          To execute, run:
    
            chimerax <scriptname> 
    Download ChimeraX
    -----------------
        ChimeraX can be found here, with stable and nightly releases : 

            https://www.cgl.ucsf.edu/chimerax/download.html

    Parameters
    ----------
        session : as is
            ChimeraX uses this to recognize to run this within the same session
            
        outFile : str
            name of output file
            
        outSuffix " str
            Suffix of output file
            Options: 
                xyz
                pdb
                cxs  (ChimeraX save file) )
                png  (default) 
                cube
                
        isoFile : str
            Input file to create the volume isosurfaces, 
            Should have suffix *_esp.cube

        colorFile : str
            Input file to use as a color map
            should have suffix *_rho.cube
            NOTE: If colorizing is not needed or colorFile is not present, please use colorFile = None 
        
        isoSurf : float
            This float instructs ChimeraX to display a certain isosurface level
            NOTE: There are disparities between appropriate values for VMD visualization and ChimeraX Visualization
            
        material : str
            This defines the material for the surface volume rendering
            Options: 
                shiny (default) 
                dull 
                
        lighting : str
            This defines the lighting of the render window
            Options:
                simple 
                full (default) 
                
        shadows : str
            Define if we want to include shadows in our render
            Options:
                True
                False (default)
                NOTE: cannot be boolean 
        
        scalemin : float
            Defines the minimum electrostatic potential for colorization
            NOTE: If both  scalemin and scalemax are set to 'compute', ChimeraX can determine the best minima and maxima to use
        
        scalemax : float
            Defines the maximum electrostatic potential for colorization
            NOTE: If both scalemin and scalemax are set to 'compute,' ChimeraX can determine the best minima and maxima to use 
    
        representation : str
            This defines how the isosurfaces should be rendered
            Options: 
                surface (default and recommended) 
                mesh 
        colorscheme : str or color1:color2:color3:color(n+1)... format for custom schemes
            This defines the colorization of the color gradient function, AKA palette in ChimeraX terminology 
            Options: 
                custom, formatted as color1:color2:color3:color(n+1) ...  
                rainbow (default) 
                red-white-blue
                grayscale
                
                please see the following for more colorization options : 
                    https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/color.html#palette-options 
    
    --------------
    End Parameters
    """
    
    run(session, 'open %s' % (isoFile))     # Open ESP
    run(session, 'volume #1 style %s level %s' % (representation , isoSurf))        # Setup Volume  
    
    # colorFile specified as None, will render without colorization 
    if colorFile == None :
        print("No ColorFile")
        if negative == True: 
            isoSurfNegative = f"-{isoSurf}" 
            run(session, 'open %s' %(isoFile)) # Open a second isoSurface file 
            run(session, 'volume #2 style %s level %s' %(representation,isoSurfNegative))

    # ColorFile Present 
    else:
        run(session, 'open %s' % (colorFile))   # Open RHO
        run(session, 'hide #2')     # Hides Colorfile from Rendering Window

        if scalemin != 'compute' and scalemax != 'compute':
            run(session, 'color gradient #1 map #2 palette %s range %s,%s' %(colorscheme, scalemin, scalemax))
        else:
            run(session, 'color gradient #1 map #2 palette %s' % (colorscheme))
        if negative == True: 
            isoSurfNegative = f"-{isoSurf}" 
            run(session, 'open %s' %(isoFile)) # Open a second isoSurface file 
            run(session, 'volume #3 style %s level %s' %(representation,isoSurfNegative))

    # Rendering Options 
    run(session, 'material %s' % (material))        # Establish surface material 
    run(session, 'lighting %s' % (lighting))        #Establish lighting and shadows
    run(session, 'lighting shadows %s' % (shadows))  # Set Shadows
    
    # SAVE  OUTFILE.OUTSUFFIX 
    run(session, 'save %s.%s' % (outFile,outSuffix)) 

    #END FUNCTION

# This Example uses 2,6-Dichloropyridine
# INPUT FILE VALUES HERE
isoFile = 'dichloropyridine26_q+0_rho.cube'
colorFile = 'dichloropyridine26_q+0_esp.cube'

# OUTPUT FILENAME AND SUFFIX HERE
outFile = isoFile[:-9]  # Removes '_esp.cube' from esp cube file to generate name
outSuffix = 'png'

#ISOSURFACE VALUES
isoSurf = .005

print_chimerax_isosurfaces(session,outFile,outSuffix,isoFile,colorFile,isoSurf,material='shiny',scalemin='compute',scalemax='compute',colorscheme='rainbow',representation='surface',lighting='full',shadows='False', negative=False)
