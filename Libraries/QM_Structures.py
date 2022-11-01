import bpy 
import math
import numpy as np

def atomic_number_to_symbol(number):
    atom_list = ['XX','H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']
    return atom_list[number]

def SubdivideSurface(objectname,subdiv=2,render=3):
    obj = bpy.data.objects[objectname]
    bpy.context.view_layer.objects.active = obj    # 'obj' is the active object now
    bpy.ops.object.modifier_add(type='SUBSURF')
    bpy.context.object.modifiers["Subdivision"].levels = subdiv
    bpy.context.object.modifiers["Subdivision"].render_levels = render
    bpy.ops.object.modifier_apply(modifier="Subdivision")

class LightObject():
    def __init__(self,name,energy=100):
        self.light_data = bpy.data.lights.new(name=f"{name}.light_data", type='POINT')
        self.light_data.energy = energy
        self.light_object = bpy.data.objects.new(name=f"{name}.light_object", object_data=self.light_data)
        bpy.context.collection.objects.link(self.light_object)
        bpy.context.view_layer.objects.active = self.light_object

class Atom():
    def __init__(self,atomic_number,index,x,y,z):
        self.element = atomic_number_to_symbol(atomic_number)
        self.index   = index
        self.x       = x
        self.y       = y
        self.z       = z
        self.r       = 0.35
        if self.element == "H":
            self.r = 0.2
        self.name = f"{self.element}.{self.index}"
        self._object = bpy.ops.mesh.primitive_uv_sphere_add(location=self.location(),radius=self.r)
        self.name_self()
        SubdivideSurface(self.name,subdiv=1,render=3)
        self.color()
        self.add_internal_light()
        self.add_to_collection("Atoms")
    def location(self):
        return [self.x,self.y,self.z]
    def color(self):
        if bpy.context.active_object.data.materials:
            bpy.context.active_object.data.materials[0] = bpy.data.materials.get( f"Elements.{self.element}")
        else:
            bpy.context.active_object.data.materials.append( bpy.data.materials.get(f"Elements.{self.element}"))
    def add_internal_light(self):
        LightObject(self.name,5)
        light_object = bpy.data.objects[f"{self.name}.light_object"]
        light_object.parent = bpy.data.objects[self.name]
    def name_self(self):
        ctx = bpy.context.copy()
        ctx['active_object'].name = self.name
    def add_to_collection(self,coll_name):
        col = bpy.data.collections.get("Collection")
        if col:
            print("found Collection")
            col.name = "SceneSupportItems"
        self.master_collection = bpy.data.collections["SceneSupportItems"]
        self.layer_collection = bpy.data.collections[coll_name]
        self.layer_collection.objects.link(bpy.context.object)
        self.master_collection.objects.unlink(bpy.context.object)
    def Update_Location(self,x,y,z,frame):
        self.x = x
        self.y = y
        self.z = z
        bpy.data.objects[f"{self.element}.{self.index}"].location=self.location()
        bpy.data.objects[f"{self.element}.{self.index}"].keyframe_insert(data_path="location", frame=frame)
#         bpy.data.lights[f"{self.element}.{self.index}.light_object"].location=self.location()
        

class Bond():
    def __init__(self,p1,p2,name):
        self.x1 = p1[0]
        self.y1 = p1[1]
        self.z1 = p1[2]
        self.x2 = p2[0]
        self.y2 = p2[1]
        self.z2 = p2[2]
        self.r  = 0.1
        self.name = name
        self.make_bond()
        self.name_self()
        self.color()
        self.add_to_collection("Bonds")
    def make_bond(self):
        dx = self.x2 - self.x1
        dy = self.y2 - self.y1
        dz = self.z2 - self.z1    
        lx = (dx/2) + self.x1
        ly = (dy/2) + self.y1
        lz = (dz/2) + self.z1
        dist = math.sqrt(dx**2 + dy**2 + dz**2)
        phi = math.atan2(dy, dx) 
        theta = math.acos(dz/dist) 
        bpy.ops.mesh.primitive_cylinder_add(radius = self.r, depth = dist, location = (lx,ly,lz) ) 
        bpy.context.object.rotation_euler[1] = theta 
        bpy.context.object.rotation_euler[2] = phi
    def name_self(self):
        ctx = bpy.context.copy()
        ctx['active_object'].name = self.name
    def color(self):
        if bpy.context.active_object.data.materials:
            bpy.context.active_object.data.materials[0] = bpy.data.materials.get( f"Bond")
        else:
            bpy.context.active_object.data.materials.append( bpy.data.materials.get(f"Bond"))
    def add_to_collection(self,coll_name):
        col = bpy.data.collections.get("Collection")
        if col:
            col.name = "SceneSupportItems"
        self.master_collection = bpy.data.collections["SceneSupportItems"]
        self.layer_collection = bpy.data.collections[coll_name]
        self.layer_collection.objects.link(bpy.context.object)
        self.master_collection.objects.unlink(bpy.context.object)
    def Update_Location(self,p1,p2,frame):
        self.x1 = p1[0]
        self.y1 = p1[1]
        self.z1 = p1[2]
        self.x2 = p2[0]
        self.y2 = p2[1]
        self.z2 = p2[2]
        dx = self.x2 - self.x1
        dy = self.y2 - self.y1
        dz = self.z2 - self.z1    
        lx = (dx/2) + self.x1
        ly = (dy/2) + self.y1
        lz = (dz/2) + self.z1
        dist = math.sqrt(dx**2 + dy**2 + dz**2)
        phi = math.atan2(dy, dx) 
        theta = math.acos(dz/dist) 
        bpy.data.objects[self.name].rotation_euler[1] = theta 
        bpy.data.objects[self.name].rotation_euler[2] = phi
        bpy.data.objects[self.name].location = [lx,ly,lz]
        bpy.data.objects[self.name].keyframe_insert(data_path="location",frame=frame)
        bpy.data.objects[self.name].keyframe_insert(data_path="rotation_euler",frame=frame)
        

def ProcessAtoms(indices,atoms,coords):
    AllAtomSpheres=[]
    coord_sum = []
    for atom in indices:
        [x,y,z] = coords[atom]
        elem = atoms[atom].atomic_number
        AllAtomSpheres.append( Atom(elem,atom,x,y,z) )
        coord_sum.append( coords[atom] )
    coord_avg = np.average( coord_sum, axis=0 )
    bpy.ops.object.empty_add( type='PLAIN_AXES', align='WORLD', location=coord_avg, scale=(0, 0, 0) )
    empty = bpy.context.object
    empty.name = "CameraTarget"
    return AllAtomSpheres

def TrackCameraToTarget():
    target_obj = bpy.data.objects['CameraTarget']
    camera_obj = bpy.data.objects['Camera']
    constraint = camera_obj.constraints.new(type='TRACK_TO')
    constraint.target = target_obj
       

def ProcessBonds(indices,bonds,coords):
    AllBondCylinders = []
    for bond in bonds:
        if bond[0] in indices:
            if bond[1] in indices:
                p1 = coords[bond[0]]
                p2 = coords[bond[1]]
                name = f"Bond.{bond[0]}.{bond[1]}"
                AllBondCylinders.append( Bond(p1,p2,name) )
    return AllBondCylinders
                
def LockMoleculeToCameraTarget():
    for obj in bpy.data.collections["Atoms"].all_objects:
        obj.parent = bpy.data.objects["CameraTarget"]
        obj.location = obj.location - bpy.data.objects["CameraTarget"].location
    for obj in bpy.data.collections["Bonds"].all_objects:
        obj.parent = bpy.data.objects["CameraTarget"]
        obj.location = obj.location - bpy.data.objects["CameraTarget"].location
                
