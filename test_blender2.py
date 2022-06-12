import bpy
import ensurepip
ensurepip.sys.path.append("/home/mark/anaconda3/envs/blender/lib/python3.7/site-packages")
ensurepip.sys.path.append("/home/mark/GitHubRepositories/BlenderMolecularOrbitals/Libraries")
from Initialize_Blender import *


from QM_Structures import *
#from DrawingFunctions import *
from ElementColors import *
#from Molden import *
#import parmed
#import numpy as np
import pytraj as pt
#import math

Initialize()

### Trajectory Data and target mask
prmtop = "/home/mark/GitHubRepositories/BlenderMolecularOrbitals/test_data/BAP_sphere_nobox.prmtop"
trajec = "/home/mark/GitHubRepositories/BlenderMolecularOrbitals/test_data/coors.dcd"
target_mask = ":1"

### Process trajectory information
traj = pt.load(trajec,prmtop)
all_atoms = [atom for atom in traj.top.atoms]
all_bonds = [bond for bond in traj.top.bond_indices]


### Get atomic indices for desired mask
atom_indices = traj.top.atom_indices(target_mask)
### Get xyz coordinates for all atoms in first frame
coords = traj[0].xyz
### Add Atoms to space
AtomicSpheres = ProcessAtoms(atom_indices,all_atoms,coords)
### Add Bonds to space
BondCylinders = ProcessBonds(atom_indices,all_bonds,coords)
### Track Camera to center of molecule.
TrackCameraToTarget()
#LockMoleculeToCameraTarget()


#### BEGIN ADDING IN TRAJECTORY DATA ###
for frame in range(len(traj)):
    for Atomic in AtomicSpheres:
        index = Atomic.index
        [x,y,z] = traj[frame].xyz[index]
        Atomic.Update_Location(x,y,z,frame)
    for BondCyl in BondCylinders:
        [idx1,idx2] = BondCyl.name.split(".")[1:]
        p1 = traj[frame].xyz[int(idx1)]
        p2 = traj[frame].xyz[int(idx2)]
        BondCyl.Update_Location(p1,p2,frame)