import pytest
from ngsolve import *
from xfem import *
from netgen.geom2d import SplineGeometry
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from math import pi

#@pytest.mark.parametrize("quad", [False])

def test_aggregates():
    square = SplineGeometry()
    square.AddRectangle([0,0],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=False))
    EA = ElementAggregation(mesh)

    return True
   