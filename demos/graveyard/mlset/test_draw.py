from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *

square = SplineGeometry()
square.AddRectangle((-4, -3), (4, 3), bc=1)
mesh = Mesh(square.GenerateMesh(maxh=2, quad_dominated=False))


level_sets = [x - 3, -x - y + 1, -x + y + 1,
              -x - 3, x - y + 1, x + y + 1 ]

level_sets_p1 = tuple(GridFunction(H1(mesh, order=1))
                      for i in range(len(level_sets)))
for i, lsetp1 in enumerate(level_sets_p1):
    InterpolateToP1(level_sets[i], lsetp1)
    Draw(lsetp1, mesh, "lsetp1_{}".format(i))

triangle = DomainTypeArray([(NEG, NEG, NEG, ANY, ANY, ANY),
                            (ANY, ANY, ANY, NEG, NEG, NEG)])

Draw(triangle.Indicator(level_sets_p1), mesh, "triangle",sd=6)
Draw(triangle.Boundary().IndicatorSmoothed(level_sets_p1), mesh, "bnd",sd=6)
Draw(triangle.Boundary().Boundary().IndicatorSmoothed(level_sets_p1), mesh, "cnr",sd=6)
