from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *

square = SplineGeometry()
square.AddRectangle((-4, -3), (4, 3), bc=1)
mesh = Mesh(square.GenerateMesh(maxh=0.5, quad_dominated=False))


level_sets = [-y+1, 2*x-y, -2*x-y]  

# level_sets = [x - 3, -x - y + 1, -x + y + 1,
#               -x - 3, x - y + 1, x + y + 1 ]

level_sets_p1 = tuple(GridFunction(H1(mesh, order=1))
                      for i in range(len(level_sets)))
for i, lsetp1 in enumerate(level_sets_p1):
    InterpolateToP1(level_sets[i], lsetp1)
    Draw(lsetp1, mesh, "lsetp1_{}".format(i))



triangle = DomainTypeArray((POS, NEG, NEG))
normals = triangle.GetOuterNormals(level_sets_p1)

for dtt in triangle.Boundary():
    name = str(dtt)
    for s in ["DOMAIN_TYPE.", ",", " ", "(", ")"]:
        name = name.replace(s, "")
    Draw(normals[(POS, NEG, NEG)][dtt], mesh, "normals_"+name)

Draw(triangle.Indicator(level_sets), mesh, "triangle")











# triangles = DomainTypeArray([(NEG, NEG, NEG, ANY, ANY, ANY),
#                              (ANY, ANY, ANY, NEG, NEG, NEG)])

# print(triangles)

# print(triangles.Boundary())
