import pytest
from timeit import default_timer as timer

from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *
from xfem import ngsxfemglobals
from xfem import IntegrateX
from xfem import Integrate
from xfem.lsetcurv import *

import numpy as np
from numpy import sqrt

# compare simd to normal Integrate
@pytest.mark.parametrize("maxh", [2**(-k) for k in range(5)])
@pytest.mark.parametrize("element_wise", [True, False])
@pytest.mark.parametrize("order", [4])
def test_simd_integrate(maxh, element_wise, order):
    square = SplineGeometry()
    square.AddRectangle((-1, -1), (1, 1), bc=1)
    ngmesh = square.GenerateMesh(maxh=maxh, quad_dominated=False)
    mesh = Mesh(ngmesh)
    
    levelset = CoefficientFunction(x**2 + y**2) - 0.5
    
    order = 1
    
    V = H1(mesh, order = order, autoupdate=True)
    lset_approx = GridFunction(V)
    InterpolateToP1(levelset, lset_approx)
    
    f = CoefficientFunction(1)
    
    dict = {
        "levelset": levelset,
        "domain_type": NEG,
        "order": order,
        "quad_dir_policy": FIRST
    }
    
    t_normal = 0
    t_simd = 0
    solutions = []
    testruns = 5
    for i in range(testruns):
        # simd disabled
        ngsxfemglobals.simd_eval = False
        start = timer()
        integral_normal = Integrate(f, mesh, VOL)
        end = timer()
        
        t_normal += (end-start)
        
        # simd enabled
        ngsxfemglobals.simd_eval = True
        
        start = timer()
        integral_simd = Integrate(f, mesh, VOL)
        end = timer()
        
        t_simd += (end-start)
        
        solutions.append(integral_simd - integral_normal)
    print()
    print("Time for Integrate (normal): " + str(t_normal/testruns))
    print("Time for Integrate (simd): " + str(t_simd/testruns))
    print("performance ratio (normal/simd): " + str(t_normal/t_simd))
    assert(np.max(solutions) < 1e-8)



# compare simd to normal IntegrateX
@pytest.mark.parametrize("maxh", [2**(-k) for k in range(5)])
@pytest.mark.parametrize("element_wise", [True, False])
@pytest.mark.parametrize("order", [4])
def test_simd_integrateX(maxh, element_wise, order):
    square = SplineGeometry()
    square.AddRectangle((-1, -1), (1, 1), bc=1)
    ngmesh = square.GenerateMesh(maxh=maxh, quad_dominated=False)
    mesh = Mesh(ngmesh)
    
    levelset = CoefficientFunction(x**2 + y**2) - 0.5
    
    V = H1(mesh, order = order, autoupdate=True)
    lset_approx = GridFunction(V)
    InterpolateToP1(levelset, lset_approx)
    
    dict = {
        "levelset": lset_approx,
        "mesh": mesh,
        "order": order,
        "domain_type": NEG
    }
    
    f = CoefficientFunction(1)
    
    t_normal = 0
    t_simd = 0
    solutions = []
    testruns = 5
    for i in range(testruns):
        # simd disabled
        ngsxfemglobals.simd_eval = False
        start = timer()
        integral_normal = IntegrateX(dict, mesh, f)
        end = timer()
        
        t_normal += (end-start)
        
        # simd enabled
        ngsxfemglobals.simd_eval = True
        
        start = timer()
        integral_simd = IntegrateX(dict, mesh, f)
        end = timer()
        
        t_simd += (end-start)
        solutions.append(integral_simd - integral_normal)
    print()
    print("Time for IntegrateX (normal): " + str(t_normal/testruns))
    print("Time for IntegrateX (simd): " + str(t_simd/testruns))
    print("performance ratio (normal/simd): " + str(t_normal/t_simd))
    assert(np.max(solutions) < 1e-7)


test_simd_integrate(0.1, False, 4)
test_simd_integrateX(0.1, False, 4)