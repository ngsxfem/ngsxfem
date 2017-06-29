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

def gen_mesh(n):
    from netgen.meshing import Mesh, MeshPoint, Element3D, FaceDescriptor, Element2D
    
    mesh = Mesh()

    pids = []
    for i in range(n+1):
        for j in range(n+1):
            for k in range(n+1):
                pids.append (mesh.Add (MeshPoint(Pnt(i/n, j/n, k/n))))
                
    for i in range(n):
        for j in range(n):
            for k in range(n):
                base = i * (n+1)*(n+1) + j*(n+1) + k
                baseup = base+(n+1)*(n+1)
                pnum = [base,base+1,base+(n+1)+1,base+(n+1),
                        baseup, baseup+1, baseup+(n+1)+1, baseup+(n+1)]
                elpids = [pids[p] for p in pnum]
                mesh.Add (Element3D(1,elpids))



    def AddSurfEls (p1, dx, dy, facenr):
        for i in range(n):
            for j in range(n):
                base = p1 + i*dx+j*dy
                pnum = [base, base+dx, base+dx+dy, base+dy]
                elpids = [pids[p] for p in pnum]
                mesh.Add (Element2D(facenr,elpids))

                            
    mesh.Add (FaceDescriptor(surfnr=1,domin=1,bc=1))
    AddSurfEls (0, 1, n+1, facenr=1)

    mesh.Add (FaceDescriptor(surfnr=2,domin=1,bc=2))
    AddSurfEls (0, (n+1)*(n+1), 1, facenr=2)

    mesh.Add (FaceDescriptor(surfnr=3,domin=1,bc=3))
    AddSurfEls (0, n+1, (n+1)*(n+1), facenr=3)


    mesh.Add (FaceDescriptor(surfnr=4,domin=1,bc=4))
    AddSurfEls ((n+1)**3-1, -(n+1), -1, facenr=1)

    mesh.Add (FaceDescriptor(surfnr=5,domin=1,bc=5))
    AddSurfEls ((n+1)**3-1, -(n+1)*(n+1), -(n+1), facenr=1)

    mesh.Add (FaceDescriptor(surfnr=6,domin=1,bc=6))
    AddSurfEls ((n+1)**3-1, -1, -(n+1)*(n+1), facenr=1)
    return mesh    


def Make3DProblem(maxh,c = (Pnt(0,0,0),Pnt(1,1,1)), quad = True):
    ngsglobals.msg_level = 0

    cube = OrthoBrick( c[0], c[1] ).bc(1)
    geom = CSGeometry()
    geom.Add (cube)
    ngmesh = gen_mesh(1) # geom.GenerateMesh(maxh=maxh, quad_dominated=quad)
    mesh = Mesh(ngmesh)
    return mesh

#levelset = sqrt(x*x+y*y)-0.6666666666 #.Compile()
levelset = 1- 2*x-2*y- 2*z

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
