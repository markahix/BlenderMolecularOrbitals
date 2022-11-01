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

## testing
scene = bpy.data.scenes["Scene"]
mycube = bpy.data.objects['CameraTarget']


bpy.ops.curve.primitive_nurbs_circle_add(radius=30, enter_editmode=False, align='WORLD', location=mycube.location, scale=(1, 1, 1))



#scene.frame_start = 1
#scene.frame_end = 360

#for frame in range(0,360):
#    mycube.rotation_euler = (0,0, (np.pi*2/300) * frame)
#    mycube.keyframe_insert('rotation_euler', index=2 ,frame=frame+1)



