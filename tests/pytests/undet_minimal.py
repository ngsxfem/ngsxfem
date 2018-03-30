from ngsolve import *
from xfem import *
from netgen.geom2d import SplineGeometry
from netgen.csg import CSGeometry, OrthoBrick, Pnt
from math import pi
import netgen.meshing as ngm
from make_uniform2D_grid import MakeUniform2DGrid
from make_uniform3D_grid import MakeUniform3DGrid

def perf_int():
    domain = POS
    order = 1

    r=0.7234436998

    levelset = sqrt(x*x+y*y+z*z)-r
    referencevals = { POS : 1-pi*r*r*r/6, NEG : pi*r*r*r/6, IF : r*r*pi/2}

    i = 1
    errors = []

    mesh = MakeUniform3DGrid(quads = True, N=int(pow(2,i)), P1=(0,0,0),P2=(1,1,1))
    print("i: " +str(i))
    print("Argument Meshing: ",str(int(pow(2,i))))
        
    V = H1(mesh,order=1)
    lset_approx = GridFunction(V)
    InterpolateToP1(levelset,lset_approx)

    f = CoefficientFunction(1)

    integral = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain},
                        cf=f, mesh=mesh, order = order)
    print("Result of Integration Reflevel ",i,", Key ",domain," : ", integral)
    return integral

if __name__ == '__main__':
    perf_int()
