from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.la import *

from numpy import linspace

from libngsxfem_py.xfem import *
from deform import *

from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters

geom = SplineGeometry("d01_testgeom.in2d")
mp = MeshingParameters (maxh=0.1)
mesh = Mesh(geom.GenerateMesh (mp=mp))

# lset_stats = StatisticContainer()
# deform_stats = StatisticContainer()
order = 4
R = 2.0/3.0
nlvl = 8
nsamps = 801
error_tab = [[None for i in range(0,nlvl)] for j in range(0,nsamps)]
error_tab2 = [[None for i in range(0,nsamps)] for j in range(0,nlvl)]
x_tab = [None for j in range(0,nsamps)]

#print("error_tab = \n", error_tab)
f = open('out', 'w')

for samp in range(0,nsamps):

  f.write ("samp = " + str(samp) + " / " + str(nsamps) + "\n")
  f.flush()
  x0 = samp/(nsamps-1)*0.1-0.05;
  print ("x0 = ", x0)
  x_tab[samp]=x0
  # x0=0
  y0=0
  levelsetstring = "sqrt((x-("+str(x0)+"))*(x-("+str(x0)+"))+(y-("+str(y0)+"))*(y-("+str(y0)+")))-(0.5+0.1*sin(8*atan2(x-("+str(x0)+"),y-("+str(y0)+"))))"
  # print(levelsetstring)
  levelset = VariableCF(levelsetstring)


  geom = SplineGeometry("d01_testgeom.in2d")
  mp = MeshingParameters (maxh=0.1)
  mesh = Mesh(geom.GenerateMesh (mp=mp))
  lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2)
  for lvl in range(0,nlvl):
    deformation = lsetmeshadap.CalcDeformation(levelset)
    maxdist = lsetmeshadap.CalcMaxDistance(levelset);
    
    # Draw(deformation)
    # Draw(lsetmeshadap.lset_p1)
   
    print ("maxdist:",maxdist)
  
    error_tab[samp][lvl] = maxdist
    error_tab2[lvl][samp] = maxdist
    lsetmeshadap.MarkForRefinement(levelset,1e-99,absolute=False)
    mesh.Refine()

  # else:
  #     lsetmeshadap.MarkForRefinement(levelset,1e-99,absolute=False)
   
import matplotlib.pyplot as plt
plt.yscale('log')
plt.plot(error_tab, "-")
plt.legend()
plt.ion()
plt.show()

print("error_tab = \n", error_tab)
print("error_tab2 = \n", error_tab2)
print("x_tab = \n", x_tab)
f.write("\nerror_tab = \n"+ str(error_tab))
f.write("\nerror_tab2 = \n"+ str(error_tab2))
f.write("\nx_tab = \n"+ str(x_tab))

f.close()

input("finished")  
