#!/usr/bin/env python3
import bpy
import ensurepip
ensurepip.sys.path.append("/home/mark/anaconda3/envs/amber/lib/python3.7/site-packages")
ensurepip.sys.path.append("/home/mark/GH_Repositories/BlenderMolecularOrbitals/Libraries")
from DrawingFunctions import *
from ElementColors import *
from Molden import *
from QM_Structures import *
import parmed
import numpy as np
import math


test_data_dir = "/home/mark/GH_Repositories/BlenderMolecularOrbitals/test_data/"

coordlines = extract_coordinate_lines(test_data_dir + "optim.molden.1")

for line in coordlines:
    subline = line.split()
    print(subline)

