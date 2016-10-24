# integration on lset domains

from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem.basics import *

NAveraging = 100

def PrintTimers(substring):
    ### print hdg-intergrator timers
    hdgtimers = [a for a in Timers() if substring in a["name"]]
    #hdgtimers = sorted(hdgtimers, key=lambda k: k["name"], reverse=False)
    hdgtimers = sorted(hdgtimers, key=lambda k: k["time"], reverse=True)
    for timer in hdgtimers:
        print("{:<45}: {:6} cts, {:.8f} s, {:.6e} s(avg.)"
              .format(timer["name"],timer["counts"],timer["time"],timer["time"]/(max(1.0,timer["counts"]))))

#square domain [-1,1]x[-1,1]
def Make2DProblem(maxh):
    from netgen.geom2d import SplineGeometry
    square = SplineGeometry()
    square.AddRectangle([-1,-1],[1,1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=maxh, quad_dominated=False))
    return mesh;

# circle with radius 0.5
levelset = (sqrt(x*x+y*y)-0.5) #.Compile()

referencevals = { POS : 4.0-0.25*pi, NEG : 0.25*pi, IF : pi }

mesh = Make2DProblem(maxh=0.1)

V = H1(mesh,order=1)
lset_approx = GridFunction(V)
InterpolateToP1(levelset,lset_approx)
#lset_approx.Set(levelset)

#domains = [NEG,POS,IF]
domains = [IF]

errors = dict()
eoc = dict()

for key in domains:
    errors[key] = []
    eoc[key] = []

#now: repetitions to average performance (later convergence..)
for reflevel in range(NAveraging):

    # if(reflevel > 0):
    #     mesh.Refine()

    f = CoefficientFunction (1.0)

    for key in domains:
        integral_old = Integrate(levelset_domain={"levelset" : levelset, "domain_type" : key}, cf=f, mesh=mesh, order=0)
        #print("\n\n ----- NOW STARTING WITH THE NEWINTEGRATEX-FUNCTION -----\n\n")
        integral = NewIntegrateX(lset=lset_approx,mesh=mesh,cf=f,order=0,domain_type=key,heapsize=1000000)

        if abs(integral_old - integral) > 1e-14:
            print("accuracy not sufficient...")
            input("press enter to continue...")

        errors[key].append(abs(integral - referencevals[key]))

for key in domains:
    eoc[key] = [log(a/b)/log(2) for (a,b) in zip (errors[key][0:-1],errors[key][1:]) ]

print("errors:  \n{}\n".format(  errors))
print("   eoc:  \n{}\n".format(     eoc))

Draw(levelset,mesh,"levelset")

PrintTimers("IntegrateX")

PrintTimers("StraightCutIntegrationRule")
#PrintTimers("StraightCutDomain")
PrintTimers("MakeQuadRuleFast")
PrintTimers("Simplex::CheckifCut")
PrintTimers("PointContainer")
PrintTimers("CutSimplex<2>::MakeQuad")
PrintTimers("MakeQuadRuleOnCutSimplex")
