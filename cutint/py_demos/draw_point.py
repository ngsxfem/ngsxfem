# integration on lset domains

from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
from netgen.csg import CSGeometry, OrthoBrick, Pnt

from xfem.lsetcurv import *

from subprocess import call

#square domain [-1,1]x[-1,1]
def Make2DProblem(maxh,c = ([0,0],[1,1]), quad = False):
    from netgen.geom2d import SplineGeometry
    square = SplineGeometry()
    square.AddRectangle(c[0],c[1],bc=1)
    mesh = Mesh (square.GenerateMesh(maxh=maxh, quad_dominated=quad))
    return mesh;

def Make3DProblem(maxh,c = (Pnt(0,0,0),Pnt(1,1,1)), quad = True):
    ngsglobals.msg_level = 0

    cube = OrthoBrick( c[0], c[1] ).bc(1)
    geom = CSGeometry()
    geom.Add (cube)
    ngmesh = geom.GenerateMesh(maxh=maxh, quad_dominated=quad)
    mesh = Mesh(ngmesh)
    return mesh

#levelset = sqrt(x*x+y*y)-0.6666666666 #.Compile()
levelset = 1-2*x-2*y-2*z

for order in range(0,3):
    call(["rm", "content_" + str(order) + ".tex"])
    DrawIntegrationPoints_HEADER(name = "content_" + str(order) + ".tex", deformorder=order, intorder=2*order)
        
    for quad,deform,c in [(True,False,(Pnt(0,0,0),Pnt(1,1,1)) )]:
        
        mesh = Make3DProblem(maxh=1.3, c=c, quad=quad)
        #mesh.Refine()
        # mesh.Refine()
        # mesh.Refine()
        
        V = H1(mesh,order=1)
        lset_approx = GridFunction(V)
        InterpolateToP1(levelset,lset_approx)
        
        # for i in range(len(lset_approx.vec)):
        #     if (lset_approx.vec[i]==0):
        #         lset_approx.vec[i] = 1e-18
        
        #lset_approx.Set(levelset)
        
        domains = [NEG,POS,IF]
        #domains = [NEG]
        
        lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=0.2, discontinuous_qn=True)
        deformation = lsetmeshadap.CalcDeformation(levelset)

        for el in deformation.space.Elements(BND):
            for i in el.dofs:
                deformation.vec[i] = 0.0
        #     print(el.dofs)

        # input("")

        if (deform):
            mesh.SetDeformation(deformation)
        
        f = CoefficientFunction (1.0)
        
        DrawIntegrationPoints(lset=lset_approx,mesh=mesh, cf=f,
                              order=2*order,name = "content_" + str(order) + ".tex")
        
        mesh.UnsetDeformation()
        # Draw(levelset,mesh,"levelset")
    DrawIntegrationPoints_FOOTER(name = "content_" + str(order) + ".tex")
    call(["pdflatex", "content_" + str(order) + ".tex"])
