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
    mesh = Mesh (square.GenerateMesh(maxh=0.1, quad_dominated=False))
    EA = ElementAggregation(mesh)

    gfu = GridFunction(H1(mesh))


    levelset = (x-0.77654)

    gfu.Set(levelset) 
    ci = CutInfo(mesh, gfu)
    roots = ci.GetElementsOfType(NEG)   
    bads = ci.GetElementsOfType(IF)

    EA.Update(roots,bads)
    #Draw(gfu)
    


    return True

if __name__ == "__main__":
    test_aggregates()