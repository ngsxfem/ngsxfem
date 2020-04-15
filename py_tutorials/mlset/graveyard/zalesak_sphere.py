"""
	Description:
	------------
	3D version of zalesak_disk.py
"""

from netgen.csg import *
from ngsolve import *
from xfem import *


geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-1,-1,-1), Pnt(1,1,1)))
mesh = Mesh(geo.GenerateMesh(maxh=0.2))


V = H1(mesh,order=1)
lsetp1_1, lsetp1_2, lsetp1_3 = GridFunction(V), GridFunction(V), GridFunction(V)

InterpolateToP1( x**2 + y**2 + z**2 - 0.8**2, lsetp1_1)
InterpolateToP1( 0.4 - y, lsetp1_2)
InterpolateToP1( 0.2 - IfPos(x,x,-x), lsetp1_3)


Draw(lsetp1_1,mesh,"lsetp1_1")
Draw(lsetp1_2,mesh,"lsetp1_2")
Draw(lsetp1_3,mesh,"lsetp1_3")

inner = {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3),
		 "domain_type": (NEG,ANY,NEG)|(NEG,NEG,ANY)
		}

boundary = {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3),
			"domain_type": (NEG,POS,IF)|(NEG,IF,NEG)|(IF,ANY,NEG)|(IF,NEG,ANY)
		   }