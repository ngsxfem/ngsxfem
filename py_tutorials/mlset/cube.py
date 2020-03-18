"""
	Description:
	------------
	Sharp cube defined by three levelsets.
"""
import math
from netgen.csg import *
from ngsolve import *
from xfem import *


geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-0.8,-0.8,-0.8), Pnt(0.8,0.8,0.8)))
mesh = Mesh(geo.GenerateMesh(maxh=1))


V = H1(mesh,order=1)
lsetp1_1, lsetp1_2, lsetp1_3 = GridFunction(V), GridFunction(V), GridFunction(V)

def cf_abs(x):
	return IfPos(x,x,-x)
InterpolateToP1( cf_abs(x) - 0.5, lsetp1_1)
InterpolateToP1( cf_abs(y) - 0.5, lsetp1_2)
InterpolateToP1( cf_abs(z) - 0.5, lsetp1_3)


Draw(lsetp1_1,mesh,"lsetp1_1")
Draw(lsetp1_2,mesh,"lsetp1_2")
Draw(lsetp1_3,mesh,"lsetp1_3")


volume = IntegrateMLsetDomain(lsets=[lsetp1_1, lsetp1_2, lsetp1_3],
                     mesh=mesh,
                     cf=1,
                     order=0,
                     domain_types=[NEG,NEG,NEG])

print("volume = {:10.8f}".format(volume))            
print("volume error = {:4.3e}".format(abs(volume - 1)))  


# inner = {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3),
# 		 "domain_type" : (NEG,NEG,NEG)
# 		}

# boundary =  {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3),
# 		 	 "domain_type" : (IF,NEG,NEG)|(NEG,IF,NEG)|(NEG,NEG,IF)
# 			}

# outer = {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3),
# 		 "domain_type" : ~(NEG,NEG,NEG)
# 		}
