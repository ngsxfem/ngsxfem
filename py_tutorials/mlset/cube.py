"""
	Description:
	------------
	Sharp cube defined by three levelsets.
"""

from netgen.csg import *
from ngsolve import *
from xfem import *


geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-0.8,-0.8,-0.8), Pnt(0.8,0.8,0.8)))
mesh = Mesh(geo.GenerateMesh(maxh=0.2))


V = H1(mesh,order=1)
lsetp1_1, lsetp1_2, lsetp1_3 = GridFunction(V), GridFunction(V), GridFunction(V)

def abs(x):
	return IfPos(x,x,-x)
InterpolateToP1( abs(x) - 0.5, lsetp1_1)
InterpolateToP1( abs(y) - 0.5, lsetp1_2)
InterpolateToP1( abs(z) - 0.5, lsetp1_3)


Draw(lsetp1_1,mesh,"lsetp1_1")
Draw(lsetp1_2,mesh,"lsetp1_2")
Draw(lsetp1_3,mesh,"lsetp1_3")

inner = {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3),
		 "domain_type" : (NEG,NEG,NEG)
		}

boundary =  {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3),
		 	 "domain_type" : (IF,NEG,NEG)|(NEG,IF,NEG)|(NEG,NEG,IF)
			}

outer = {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3),
		 "domain_type" : ~(NEG,NEG,NEG)
		}
