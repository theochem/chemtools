from chimerax.core.commands import run
from chimerax.core import errors
from chimerax.map import Volume,VolumeSurface
from chimerax import io
"""
THIS SCRIPT IS DESIGNED TO TAKE GAUSSIAN CUBE FILES AND VISUALIZE THE FILES IN CHIMERAX

TO EXECUTE SCRIPT RUN : 

chimerax <scriptname> 

ChimeraX can be found here, with stable and nightly releases : 

https://www.cgl.ucsf.edu/chimerax/download.html


PARAMETERS : 
    session, as is
        ChimeraX uses this to recognize to run this within the same session

    outFile, str
        name of output file

    outSuffix, str
        Suffix of output file
        Options: 
            XYZ
            PDB
            CXS (DEFAULT CHIMERAX SAVE STATE) 
            PNG
            CUBE

    isoFile, str
        Input file to create the volume isosurfaces, 
        Should have suffic *_esp.cube

    colorFile
        Input file to use as a color map
        should have suffix *_rho.cube

    isoSurf, float
       This float instructs ChimeraX to display a certain isosurface level

    material, str
        define the material for the surface volume rendering
        OPTIONS: 
            shiny
            dull 

    lighting, str
        define the lighting of the render window
        OPTIONS:
            simple
            full

    shadows, str
        Define if we want to include shadows in our render
        OPTIONS
            True
            False
            NOTE: cannot be boolean 

    scalemin, float
        Defines the minimum electrostatic potential for colorization

    scalemax, flaot
        Defines the maximum electrostatic potential for colorization

    representation, str
        defines how the isosurfaces should be rendered
        OPTIONS: 
            surface (default and recommended) 
            mesh 

    colorscheme, str, or color1:color2:color3:color(n+1)... format for custom schemes
        defines the colorization of the color gradient function, AKA palette in ChimeraX terminology 
        OPTIONS: 
            Custom
            rainbow
            red-white-blue
            grayscale
            please see https://www.rbvi.ucsf.edu/chimerax/docs/user/commands/color.html#palette-options for more options 
    
    END PARAMETERS

"""


def print_chimerax_isosurfaces(session,outFile,outSuffix, isoFile, colorFile, isoSurf, material, scalemin,scalemax,colorscheme,representation, lighting, shadows,silhouettes): 
    run(session, 'open %s' % (isoFile))     # Open ESP
    run(session, 'open %s' % (colorFile))   # Open RHO 
    run(session, 'hide #2')     # Hides Colorfile from Rendering Window
    run(session, 'volume #1 style %s level %s' % (representation , isoSurf))        # Setup Volume  
    run(session, 'color gradient #1 map #2 palette %s range %s,%s' % (colorscheme , scalemin , scalemax)) 
    run(session, 'material %s' % (material))        # establish surface material 
    run(session, 'lighting %s' % (lighting))        #establish lighting and shadows
    run(session, 'lighting shadows %s' % (shadows))
    
    # SAVE  OUTFILE.OUTSUFFIX 
    run(session, 'save %s.%s' % (outFile,outSuffix)) 

    #END FUNCTION

# INPUT FILE VALUES HERE
isoFile = 'dichloropyridine26_q+0_esp.cube'
colorFile = 'dichloropyridine26_q+0_rho.cube'

# OUTPUT FILENAME AND SUFFIX HERE
outFile = isoFile[:-9]  # Removes '_esp.cube' from esp cube file to generate name
outSuffix = 'png'      

#THESE ARE DEFAULT VALUES
isoSurf = .055
material = 'shiny' 
colorscheme = '#b800ff:#7500ff:#6500ff:#5405ff:#3000ff:#2934ff:#31ff38:#ca1313:#ca1111:#ca1818' # THIS FORMAT IS USERMADE AND NOT FOUND IN DOCUMENTATION 
representation = 'surface'
lighting = 'full'
shadows = 'False'
scalemin = -.02 
scalemax = .02


# CALL FUNCTION HERE
print_chimerax_isosurfaces(session,outFile,outSuffix,isoFile,colorFile,isoSurf,material,scalemin,scalemax,colorscheme,representation,lighting,shadows,silhouettes)
