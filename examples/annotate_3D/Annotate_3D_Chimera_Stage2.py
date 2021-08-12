from chimerax.core.commands import run

# END DEPENDENCIES
"""
THIS SCRIPT IS STAGE 2 OF A 3D MOLECULAR ANNOTATION PIPLINE
TO GENERATE THE INTERMEDIATE ATTRIBUTE FILE AND PDB FILE, RUN STAGE 1 SCRIPT FIRST

THIS SCRIPT USES CHIMERAX AS A VISUALIZATION PLATFORM DISTRIBUTED BY THE UNIVERSITY OF CALIFORNIA SAN FRANCISCO

CHIMERAX DAILY AND STABLE BUILDS AVAILABLE HERE:


TO EXECUTE THIS SCRIPT, RUN IN COMMAND LINE : 
        
    chimerax <filename> 


PARAMETERS

    session
        This instructs ChimeraX to run all commands within the same environment, and does not need to be initialized
    
    moleculeFile, string
        This string represents the input molecular file
    
    annotationFile, string ,must have *.defattr as suffix
        This represents then name of the input annotation file
    
    attributeName, string
        This is the name of the attribute factor being used, which can be found in the first line of the attributeFile

    outputFileName, string
        This represents the output file name 

    outputFileType, string
        This represents the output file type
            OPTIONS
                CXS - CHIMERAX SAVE STATE
                PNG 
    style, string
        This represents how the molecule should be rendered
            OPTIONS: 
                stick (default)
                sphere
                ball-and-stick

    units, string
        This specifies that we are working with atoms, 

    lighting, string
        This specifies the lighting of the ChimeraX render window 
            OPTIONS: 
                simple (default) 
                full

END PARAMETERS 
"""

def annotateChimerax(session,moleculeFile,annotationFile,attributeName,outputFileName,outputFileType,style,units,lighting):
    run(session,'open %s' % (moleculeFile))
    run(session,'open %s' % (annotationFile)) 
    run(session,'label #1 %s attribute %s' %(units, attributeName) )
    run(session,'lighting %s' %(lighting)) 
    run(session,'save %s.%s' %(outputFileName, outputFileType)) 
 
# INPUT FILES
moleculeFile = 'dichloropyridine26_q+0.pdb' 
annotationFile = 'test.defattr' 

# ATTRIBUTE NAME 
attributeName = 'myFactor' 

#OUTPUT FILES
outputFileName = moleculeFile[:-4]
outputFileType = 'cxs' 
 
# VISUALIZATION PARAMETERS
style = 'stick' 
units = 'atoms'
lighting = 'simple' 


# ACTIVATE COMMAND
annotateChimerax(session,moleculeFile,annotationFile,attributeName,outputFileName,outputFileType,style,units,lighting) 
