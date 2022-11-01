import bpy

def Create_Or_Find_Collection(coll_name):
    collectionFound = False
    for myCol in bpy.data.collections:
        print(myCol.name)
        if myCol.name == coll_name:
            collectionFound = True
            break
    if collectionFound == False:
        myCol = bpy.data.collections.new(coll_name)
        bpy.context.scene.collection.children.link(myCol) #Creates a new collection

def Initialize():
    col = bpy.data.collections.get("Collection")
    if col:
        print("found Collection")
        col.name = "SceneSupportItems"
    Create_Or_Find_Collection("Atoms")
    Create_Or_Find_Collection("Bonds")
    bpy.context.scene.render.engine = 'CYCLES'
    bpy.context.scene.cycles.device = 'GPU'
    bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = (0, 0, 0, 1)

