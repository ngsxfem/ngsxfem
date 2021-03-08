"""
This example (visually) illustrates the usage of the usage and
functionality of xfem.mlset to generate domain descriptions based on
multiple, simpler (shorter) multi-level set domain descriptions.
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem.mlset import *


# ----------------------------------- MAIN ------------------------------------
square = SplineGeometry()
square.AddRectangle((-2.5, -1), (2.5, 2.5), bc=1)
mesh = Mesh(square.GenerateMesh(maxh=1, quad_dominated=False))


# Separate regions
level_sets1 = [-2*x + y - 1, 2*x + y - 1, y - 2]
dta1 = DomainTypeArray((POS, POS, NEG))
Draw(dta1.Indicator(level_sets1), mesh, "dta1", sd=4)

level_sets2 = [-2*x + y - 4, 2*x + y + 2, y]
dta2 = DomainTypeArray((NEG, NEG, POS))
Draw(dta2.Indicator(level_sets2), mesh, "dta2", sd=4)

level_sets3 = [x - 2*y - 1, x + 2*y - 1, x - 2]
dta3 = DomainTypeArray((POS, POS, NEG))
Draw(dta3.Indicator(level_sets3), mesh, "dta3", sd=4)


# Test TensorUnion
dta_union = TensorUnion(dta1, dta2, dta3)
Draw(dta_union.Indicator(level_sets1 + level_sets2 + level_sets3),
     mesh, "dta_union", sd=5)

dta4 = DomainTypeArray([(POS, POS, NEG, ANY, ANY, ANY),
                        (ANY, ANY, ANY, NEG, NEG, POS)])
dta_union2 = TensorUnion(dta4, dta3)
Draw(dta_union.Indicator(level_sets1 + level_sets2 + level_sets3),
     mesh, "dta_union2", sd=5)


# Test TensorIntersection
dta5 = ~ dta1
dta6 = ~ dta2
dta7 = ~ dta3

dta_intersection = TensorIntersection(dta5, dta6, dta7)
Draw(dta_intersection.Indicator(level_sets1 + level_sets2 + level_sets3),
     mesh, "dta_intersection", sd=5)

bnd1 = dta1.Boundary()
dta_intersection2 = TensorIntersection(bnd1, dta6)
Draw(dta_intersection2.IndicatorSmoothed(level_sets1 + level_sets2),
     mesh, "dta_intersection2", sd=6)
