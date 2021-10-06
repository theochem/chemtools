def chimerax_isosurfaces(outFile,outSuffix, isoFile, colorFile, isoSurf, material='shiny', scalemin='compute',scalemax='compute',colorscheme='rainbow',representation='surface',lighting='full',shadows='False',negative=False): 
    """
    """
    print("Writing ChimeraX Visualization Script")
    outScript = f"{outFile}_chimeraxScript.py"                                          # output script filename initialize
    file = open(f"{outScript}","w")                                                     # open outputScript 

    file.write("from chimerax.core.commands import run \n") 
    file.write("from chimerax.core import errors \n") 
    file.write("from chimerax import io \n") 
    
    file.write(f"run(session, 'open {isoFile}') \n") 
    file.write(f"run(session, 'volume #1 style {representation} level {isoSurf}') \n") 

    
    # colorFile specified as None, will render without colorization 
    if colorFile == None :
        print("\n",
            "*********************************\n",
            "WARNING: NO COLORFILE PROVIDED!!!\n",
            "THIS IS NOT NECESSARILY AN ERROR\n",
            "*********************************\n")
        if negative == True: 
            isoSurfNegative = f"-{isoSurf}" 
            file.write(f"run(session, 'open {isoFile}') \n")                                                    # Open a second isoSurface file 
            file.write(f"run(session, 'volume #2 style {representation} level {isoSurfNegative}')\n")
    else:
        file.write(f"run(session, 'open {colorFile}')\n")  # Open RHO
        file.write(f"run(session, 'hide #2')\n")      # Hides Colorfile from Rendering Window

        if scalemin != 'compute' and scalemax != 'compute':
            file.write(f"run(session, 'color gradient #1 map #2 palette {colorscheme} range %{scalemin},{scalemax}') \n") 
        else:
            file.write(f"run(session, 'color gradient #1 map #2 palette {colorscheme}') \n")
        if negative == True: 
            isoSurfNegative = f"-{isoSurf}" 
            file.write(f"run(session, 'open {isoFile}') \n") # Open a second isoSurface file 
            file.write(f"run(session, 'volume #3 style {representation} level {isoSurfNegative}') \n")
    # Rendering Options 
    file.write(f"run(session, 'material {material}') \n")       # Establish surface material 
    file.write(f"run(session, 'lighting {lighting}') \n")         # Establish lighting and shadows
    file.write(f"run(session, 'lighting shadows {shadows}') \n")   # Set Shadows
    
    # SAVE  OUTFILE.OUTSUFFIX   
    file.write(f"run(session, 'save {outFile}.{outSuffix}') \n") 

    print(f"ChimeraX Visualization Script {outScript} Written !") 
    file.close
    print(f"You can now execute {outScript} by running 'chimerax {outScript}'")
    print("Exiting...")
    # END FUNCTION

# BEGIN EXAMPLE 

isoFile= 'dichloropyridine26_q+0_rho.cube'
colorFile = 'dichloropyridine26_q+0_esp.cube'
# OUTPUT FILENAME AND SUFFIX HERE
outFile = isoFile[:-9]  # Removes '_esp.cube' from esp cube file to generate name
outSuffix = 'png'

#ISOSURFACE VALUES
isoSurf = .005


chimerax_isosurfaces(outFile,outSuffix, isoFile, colorFile, isoSurf, material='shiny', scalemin='compute',scalemax='compute',colorscheme='rainbow',representation='surface',lighting='full',shadows='False',negative=False)
# EOF 
