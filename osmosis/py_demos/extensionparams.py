
from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.ngstd import *
from ngsolve.la import *

from osmosis import *
from roundedsquare import *


dim = 2


# from ngsolve import *
# from ngsolve.internal import *
# from netgen.geom2d import unit_square

# ngsglobals.msg_level = 1

# mesh = Mesh( unit_square.GenerateMesh(maxh=1))
# mesh.Refine()
# mesh.Refine()
# mesh.Refine()

# Draw(mesh)

mesh = Mesh("square6.vol.gz")

#lset = RoundedSquareLevelSet(x0=0.5,y0=0.5,d=0.3,R=0.03)
#lset = VariableCF("sqrt(x*x+y*y)-0.3")

#star fish (on [-1,1]x[-1,1]-domain):
# lset = VariableCF("sqrt(x*x+y*y)-(0.45+0.18*sin(5*atan2(x,y)))")

lset = VariableCF("sqrt(x*x+4*y*y)-0.5")

#dim = 3
#mesh = Mesh("cube_big_16.vol.gz")
#R = 0.45
#r = 0.2
#lset = VariableCF("(x*x+y*y+z*z+"+str(R)+"*"+str(R)+"-"+str(r)+"*"+str(r)+")*(x*x+y*y+z*z+"+str(R)+"*"+str(R)+"-"+str(r)+"*"+str(r)+")-4*"+str(R)+"*"+str(R)+"*(x*x+y*y)");
#lset = VariableCF("sqrt((x*x+y*y+z*z))-1")

alpha = 1.0
beta  = 0.3

dt = 1e-4
tend = 1.015
theta = 0.5

#initialu=VariableCF("4*(1+x)")
#initialu=VariableCF("0.05*(25-x*x-y*y-z*z)")
initialu = VariableCF("(x+1)*(x+1)")
