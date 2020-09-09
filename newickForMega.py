# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 15:29:35 2020

@author: 
    
Description:
    Correct the Newick tree apostrophes from "\'" to "''" to be readable with
    Mega (otherwise it raise a error).
"""

import re

def correct_newick(file_path):
    
    # open file
    with open(file_path, mode='r') as f:
        file = f.read()
        # replace the  apostrophes
        newfile = re.sub(r"\\'", "''", file)
    
    # overwrite the file
    with open(file_path, mode='w') as j:
        j.write(newfile)
    