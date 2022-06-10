import bpy

class Atom():
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
        if self._object.data.materials:
            self._object.data.materials[0] = bpy.data.materials.get(f"Elements.{self.element}")
        else:
            self._object.data.materials.append(bpy.data.materials.get(f"Elements.{self.element}"))
