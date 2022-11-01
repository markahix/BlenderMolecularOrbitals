import math

def cylinder_between(x1, y1, z1, x2, y2, z2, r):
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1    
    dist = math.sqrt(dx**2 + dy**2 + dz**2)
    phi = math.atan2(dy, dx) 
    theta = math.acos(dz/dist) 
    return_string = f"""
bpy.ops.mesh.primitive_cylinder_add(
  radius = {r}, 
  depth = {dist},
  location = ({dx/2+x1}, {dy/2+y1}, {dz/2+z1})   
) 
bpy.context.object.rotation_euler[1] = {theta} 
bpy.context.object.rotation_euler[2] = {phi}
"""
    return return_string

def draw_atoms(traj,frame_index):
    elements = ['XX','H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo']
    atoms = [atom for atom in traj.top.atoms]
    coords = traj[frame_index].xyz
    return_string = """class Atom():
    def __init__(self,element,index,x,y,z):
        self.element = element
        self.index   = index
        self.x       = x
        self.y       = y
        self.z       = z
        self._object = bpy.ops.mesh.primitive_uv_sphere_add(location=self.location(),radius=1.0)
        self.color()
    def location(self):
        return [self.x,self.y,self.z]
    def color(self):
        if bpy.context.active_object.data.materials:
            bpy.context.active_object.data.materials[0] = bpy.data.materials.get( f"Elements.{self.element}")
        else:
               bpy.context.active_object.data.materials.append( bpy.data.materials.get(f"Elements.{self.element}"))
"""
    for atom in atoms:
        element = elements[int(atom.atomic_number)]
        index = int(atom.index)
        [x,y,z] = coords[index]
        return_string += f"Atom('{element}',{index},{x},{y},{z})\n"
        return_string += f"print('{index}.{element} @ [{x},{y},{z}]')\n"
    return return_string

def draw_bonds(traj,frame_index):
    bond_pairs = [list(bond.indices) for bond in traj.top.bonds]
    coords = traj[frame_index].xyz
    full_string = ""
    for bond in bond_pairs:
        p1 = coords[bond[0]]
        p2 = coords[bond[1]]
        full_string += cylinder_between(p1[0],p1[1],p1[2],p2[0],p2[1],p2[2],0.1)
        full_string += f"print('{bond}')"
    return full_string
