"""
	Description:
	-----------
	Zalesak Disk described my multiple levelsets. Taken from [1].

	References:
	-----------
	[1] D. P. Starinshak, S. Karni and P. L. Roe: A New Level-Set
		Model for the Representation of Non-Smooth Geometries. In J Sci 
		Comput (2014) 61(3):649-672. DOI:10.1007/s10915-014-9842-0

"""

from netgen.geom2d import SplineGeometry
from ngsolve import *
from xfem import *

geo = SplineGeometry()
geo.AddRectangle((-1.3, -1.3), (1.3, 1.3), bc=1)

mesh = Mesh(geo.GenerateMesh(maxh=0.2))



V = H1(mesh,order=1)
lsetp1_1, lsetp1_2 = GridFunction(V), GridFunction(V)
lsetp1_3, lsetp1_4 = GridFunction(V), GridFunction(V)

InterpolateToP1(x * x + y * y - 1, lsetp1_1)
InterpolateToP1(-x - 1 / 3, lsetp1_2)
InterpolateToP1(x - 1 / 3, lsetp1_3)
InterpolateToP1(y - 0.5, lsetp1_4)

Draw(lsetp1_1, mesh, "lsetp1_1")
Draw(lsetp1_2, mesh, "lsetp1_2")
Draw(lsetp1_3, mesh, "lsetp1_3")
Draw(lsetp1_4, mesh, "lsetp1_4")
SetVisualization(min=0,max=0)
Redraw()

inner = {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3, lsetp1_4)
		 "domain_type": (NEG,NEG,NEG,POS)|(NEG,POS,NEG,POS)|(NEG,POS,NEG,NEG)\
		 				|(NEG,NEG,POS,NEG)|(NEG,NEG,POS,POS) 
		}

# Alternative 1
inner1 = {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3, lsetp1_4)
		 "domain_type": (NEG,ANY,ANY,POS)|(NEG,POS,ANY,ANY)|(NEG,ANY,POS,ANY)
		 }

# Alternative 2
inner2 = {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3, lsetp1_4)
		 "domain_type": (NEG,ANY,ANY,ANY)& ~(NEG,NEG,NEG,NEG)
		 }

# Alternative to 2
region1 = {"levelset": lsetp1_1, "domain_type": NEG}
region2 = {"levelsets": (lsetp1_1, lsetp1_2, lsetp1_3, lsetp1_4),
		   "domain_type": (NEG,NEG,NEG,NEG)
		  }
inner4 = region1 - region2