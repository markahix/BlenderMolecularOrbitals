import bpy

Jmol_color_data ="""1	H	[255,255,255]
2	He	[217,255,255]
3	Li	[204,128,255]
4	Be	[194,255,0]
5	B	[255,181,181]
6	C	[44,44,44]
7	N	[48,80,248]
8	O	[255,13,13]
9	F	[144,224,80]
10	Ne	[179,227,245]
11	Na	[171,92,242]
12	Mg	[138,255,0]
13	Al	[191,166,166]
14	Si	[240,200,160]
15	P	[255,128,0]
16	S	[255,255,48]
17	Cl	[31,240,31]
18	Ar	[128,209,227]
19	K	[143,64,212]
20	Ca	[61,255,0]
21	Sc	[230,230,230]
22	Ti	[191,194,199]
23	V	[166,166,171]
24	Cr	[138,153,199]
25	Mn	[156,122,199]
26	Fe	[224,102,51]
27	Co	[240,144,160]
28	Ni	[80,208,80]
29	Cu	[200,128,51]
30	Zn	[125,128,176]
31	Ga	[194,143,143]
32	Ge	[102,143,143]
33	As	[189,128,227]
34	Se	[255,161,0]
35	Br	[166,41,41]
36	Kr	[92,184,209]
37	Rb	[112,46,176]
38	Sr	[0,255,0]
39	Y	[148,255,255]
40	Zr	[148,224,224]
41	Nb	[115,194,201]
42	Mo	[84,181,181]
43	Tc	[59,158,158]
44	Ru	[36,143,143]
45	Rh	[10,125,140]
46	Pd	[0,105,133]
47	Ag	[192,192,192]
48	Cd	[255,217,143]
49	In	[166,117,115]
50	Sn	[102,128,128]
51	Sb	[158,99,181]
52	Te	[212,122,0]
53	I	[148,0,148]
54	Xe	[66,158,176]
55	Cs	[87,23,143]
56	Ba	[0,201,0]
57	La	[112,212,255]
58	Ce	[255,255,199]
59	Pr	[217,255,199]
60	Nd	[199,255,199]
61	Pm	[163,255,199]
62	Sm	[143,255,199]
63	Eu	[97,255,199]
64	Gd	[69,255,199]
65	Tb	[48,255,199]
66	Dy	[31,255,199]
67	Ho	[0,255,156]
68	Er	[0,230,117]
69	Tm	[0,212,82]
70	Yb	[0,191,56]
71	Lu	[0,171,36]
72	Hf	[77,194,255]
73	Ta	[77,166,255]
74	W	[33,148,214]
75	Re	[38,125,171]
76	Os	[38,102,150]
77	Ir	[23,84,135]
78	Pt	[208,208,224]
79	Au	[255,209,35]
80	Hg	[184,184,208]
81	Tl	[166,84,77]
82	Pb	[87,89,97]
83	Bi	[158,79,181]
84	Po	[171,92,0]
85	At	[117,79,69]
86	Rn	[66,130,150]
87	Fr	[66,0,102]
88	Ra	[0,125,0]
89	Ac	[112,171,250]
90	Th	[0,186,255]
91	Pa	[0,161,255]
92	U	[0,143,255]
93	Np	[0,128,255]
94	Pu	[0,107,255]
95	Am	[84,92,242]
96	Cm	[120,92,227]
97	Bk	[138,79,227]
98	Cf	[161,54,212]
99	Es	[179,31,212]
100	Fm	[179,31,186]
101	Md	[179,13,166]
102	No	[189,13,135]
103	Lr	[199,0,102]
104	Rf	[204,0,89]
105	Db	[209,0,79]
106	Sg	[217,0,69]
107	Bh	[224,0,56]
108	Hs	[230,0,46]
109	Mt	[235,0,38]"""

class Material():
    def __init__(self,name,color_list):
        self.r = float(color_list[0])/255
        self.g = float(color_list[1])/255
        self.b = float(color_list[2])/255
        self.a = float(1.0)
        self.name = name
        self.create_material()
    def create_material(self):
        bpy.data.materials.new(self.name)
#         bpy.data.materials[self.name].diffuse_color = ( self.r, self.g, self.b, self.a )
        bpy.data.materials[self.name].use_nodes = True
        self.nodes = bpy.data.materials[self.name].node_tree.nodes
        self.bsdf = self.nodes.get("Principled BSDF")
        self.bsdf.inputs["Base Color"].default_value = ( self.r, self.g, self.b, self.a )
    def add_emission(self,magnitude = 1.0):
        mat = bpy.data.materials[self.name]
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links
        output = nodes.new('ShaderNodeOutputMaterial')
        shader = nodes.new('ShaderNodeEmission')
        nodes["Emission"].inputs[0].default_value = (self.r, self.g, self.b, self.a)
        nodes["Emission"].inputs[1].default_value = magnitude
        links.new(shader.outputs[0], output.inputs[0])
    def add_transparency(self,magnitude=0.0):
        self.bsdf.inputs["Transmission"].default_value = magnitude
        
for line in Jmol_color_data.split("\n"):
    subline = line.split()
    element = subline[1]
    colors  = subline[2].strip('][').split(',')
    mat = Material(f"Elements.{element}",colors)
    mat.add_transparency(0.25)
    
bond = Material("Bond",[255,255,255])
bond.add_emission(0.00015)
