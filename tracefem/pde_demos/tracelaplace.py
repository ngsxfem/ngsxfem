
# interactive modifications to d1_approx.py
print ("hello from d1_approx.py ;-)")

from ngsolve.solve import *
from ngsolve.comp import *
from ngsolve.fem import *
from ngsolve.la import *
from ngsolve.bla import *
import ngsolve.ngstd as ngstd
from ngsolve.solve import Redraw

#from libngsxfem_py.xfem import *
import libngsxfem_py.xfem as xfem                                 

from time import sleep

mesh = Mesh("square2.vol.gz")

v = FESpace ("h1ho", mesh, order=1, dirichlet=[1,2,3,4])
v.Update()

lset = CoefficientFunction ()

u = GridFunction (v, name="u")
u.Update()

f = LinearForm (v)
f.Add (LFI (name = "source", dim = 2, coef = ConstantCF(1),
            flags = { }, definedon = [0]))

a = BilinearForm (v, flags = { "symmetric" : True, 
                               "eliminate_internal" : False })
a.Add (BFI ("laplace", 2, ConstantCF(1)))

