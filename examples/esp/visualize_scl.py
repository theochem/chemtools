import os,sys

from chimerax.core.commands import run, Command, CmdDesc, RestOfLine
from chimerax.core import errors
 
from chimerax.map import Volume 
from chimerax.map_data.gaussian.gaussian_grid import GaussianGrid 
from chimerax.surface import color_electrostatic 

# END MODULES
"""
THIS SCRIPT FIRST OPENS ESP AND RHO CUBE FILES, VISUALIZES THE ISOSURFACES OF THE ESP CUBE FILE USING VOLUME VIEWER

I RETURN ERRORS WHEN TRYING TO IMPLEMENT STRING LITERALS SO THAT CUBEFILE NAMES CAN BE VARIABLES

""" 
# NONFUNCTIONALL
espName = 'dichloropyridine26_q+0_esp.cube'
rhoName = 'dichloropyridine26_q+0_esp.cube'


def VisualizeMolecule(session,espName,rhoName): 

    run(session, "open scl2_esp.cube")
    run(session, "open scl2_rho.cube")
    run(session, "volume #1 style surface level .038") 
    run(session, "lighting simple")
    run(session, "graphics silhouettes false")
    run(session, "lighting shadows false")  
    #run(session, "color electrostatic #1.1#1.2 map #2 palette 0,#ff0000:0.025,#ffffff:0.05,#0000ff")  
    run(session, "hide #2 model") 
# END VISUALIZEMOLECULE 

VisualizeMolecule(session,espName,rhoName)


"""
TODO:

COLOR MAPPING FUNCTION
SYSARGS AND STRING LITERALS 
"""
