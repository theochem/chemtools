from chimerax.core.commands import run
from chimerax.core import errors
 
# END MODULES
"""
THIS SCRIPT FIRST OPENS ESP AND RHO CUBE FILES, VISUALIZES THE ISOSURFACES OF THE ESP CUBE FILE USING VOLUME VIEWER

WE NEED TO FIND A MEANS OF IMPLEMENTING FSTRING{rhoName} etc.  


TO EXECUTE THIS SCRIPT RUN 
' chimerax visualize.py 

CURRENTLY THE INPUT *_esp.cube and *_rho.cube inputs, ISOSURFACE LEVELS, AND OUTPUT FILENAMES ARE HARD CODED 

""" 


def VisualizeMolecule(session): 

    run(session, "open dichloropyridine26_q+0_esp.cube")
    run(session, "open dichloropyridine26_q+0_rho.cube")
    run(session, "volume #1 style surface level .0555") 
    run(session, "lighting simple")
    run(session, "graphics silhouettes false")
    run(session, "lighting shadows false")  
    run(session, "hide #2 model") 
# END VISUALIZEMOLECULE 

def ColorElectrostaticChimera(session): 
    run(session, "hide #2 model") 
    run(session, "color gradient #1 map #2 palette rainbow")



def SaveMoleculeChimera(session):
    run(session, "hide #2 model") 
    run(session, "save outFile.cxs") #CHIMERAX SAVE STATE 
    run(session, "save outFile.png") # SAVE AS PNG IMAGE    


#END SAVEMOLECULECHIMERA 
VisualizeMolecule(session)
ColorElectrostaticChimera(session)
SaveMoleculeChimera(session) 
