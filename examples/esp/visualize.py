import os,sys
from chimerax.core.commands import run, Command, CmdDesc, RestOfLine
from chimerax.core import errors
 
# END MODULES
"""
THIS SCRIPT FIRST OPENS ESP AND RHO CUBE FILES, VISUALIZES THE ISOSURFACES OF THE ESP CUBE FILE USING VOLUME VIEWER

I RETURN ERRORS WHEN TRYING TO IMPLEMENT STRING LITERALS SO THAT CUBEFILE NAMES CAN BE VARIABLES


WE NEED TO FIND A MEANS OF IMPLEMENTING FSTRING{rhoName} etc.  
""" 

espName = 'dichloropyridine26_q+0_esp.cube'

rhoName = 'dichloropyridine26_q+0_esp.cube'

colorFile = rhoName # ?????? 
isosurfaceValue = '.0555'


def VisualizeMolecule(session,espName,rhoName): 

    run(session, "open dichloropyridine26_q+0_esp.cube")
    run(session, "open dichloropyridine26_q+0_rho.cube")
    run(session, "volume #1 style surface level .0555") 
    run(session, "lighting simple")
    run(session, "graphics silhouettes false")
    run(session, "lighting shadows false")  
    run(session, "hide #2 model") 
# END VISUALIZEMOLECULE 

def ColorElectrostaticChimera(session,espName,rhoName,colorFile): 
    run(session, "hide #2 model") 
    run(session, "color gradient #1 map #2 palette rainbow")



def SaveMoleculeChimera(session):
    run(session, "hide #2 model") 
    run(session, "save outFile.cxs") #CHIMERAX SAVE STATE 
    run(session, "save outFile.png") # SAVE AS PNG IMAGE    


#END SAVEMOLECULECHIMERA 
VisualizeMolecule(session,espName,rhoName)
ColorElectrostaticChimera(session,espName,rhoName,colorFile)
SaveMoleculeChimera(session) 

"""
TODO:

COLOR MAPPING FUNCTION
SYSARGS AND STRING LITERALS 
"""
