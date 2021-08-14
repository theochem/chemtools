#!/usr/bin/env python
# coding: utf-8

# In[15]:


import pandas as pd
import numpy as np 

import iodata
from iodata import load_one, dump_one 

from biopandas.pdb import PandasPdb
import io
from typing import TextIO, Iterator


# In[23]:


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


workingBFactor = loadInput.atcharges['esp'][:]
print("Our Working B Factor to replace within IOData Dump : ")
print(workingBFactor)
print(f"\nWill be written to {outFile} and dumped")


# In[27]:



workingBfactor = np.array(loadInput.atcharges['esp'][:])
print(workingBfactor)
test = load_one('intermediate.pdb')

for i,item in enumerate(test.extra['bfactors']):
test.extra['bfactors'][i] = workingBfactor[i]
print(" ")
print(f"Written {test.extra['bfactors'][:]} to {outFile}")
dump_one(test,'test22.pdb')


# In[ ]:





# In[ ]:




