from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
from sympy import *

from integrate_one_big_quad2D import get_levelset, get_referencevals

if __name__ == "__main__":
    from netgen.geom2d import SplineGeometry
    square = SplineGeometry()
    square.AddRectangle([0,0],[1,1],bc=1)
    #mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))

    lsetvals_list = [ [-0.18687,0.324987, 0.765764,0.48983], [0.765764,0.324987, -0.18687, -0.48983], [1,2/3,-1,-2/3]]
    #lsetvals_list = [[1.,-2.,-2.,-2.]]
    #lsetvals_list = [[1,-1,-4,-2]]
    lsetvals_list.append([3,-1,1,-1.023123])
    
    n_ref = 8
    order = 3
    f = lambda x,y: 1+0*x+0*y
    f_ngs = f(x,y)

    domains = [NEG,POS]

    for lsetvals in lsetvals_list:
        print("Case lsetvals = ", lsetvals)
        referencevals = get_referencevals(lsetvals, f)
        levelset =get_levelset(lsetvals)

        errors = dict()
        eoc = dict()
        for key in domains:
            errors[key] = []
            eoc[key] =[]

        mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))
        for i in range(n_ref):
            V = H1(mesh,order=1)
            lset_approx = GridFunction(V)
            InterpolateToP1(levelset,lset_approx)

            for key in domains:
                integral = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=f_ngs,order=order,domain_type=key,heapsize=1000000, use_saye=True)
                errors[key].append(abs(integral - referencevals[key]))
            mesh.Refine()
        for key in domains:
            eoc[key] = [log(errors[key][i+1]/errors[key][i])/log(0.5) for i in range(n_ref-1)]
        print("L2 Errors:", errors)
        print("experimental order of convergence (L2):", eoc)

        #Draw(levelset, mesh, "lset")
