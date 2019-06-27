# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from xfem import *

ngsglobals.msg_level = 1

# generate a triangular mesh of mesh-size 0.2
square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.2, quad_dominated=False))

mu = [1, 1e-7]
order = 1


rhsscal = -4*mu[0]*mu[1] 
coef_f = [ rhsscal, rhsscal ]
distsq = x*x + y*y - 1
solution = [ mu[1]* distsq, mu[0] * distsq ]

levelset = ( sqrt(x*x + y*y) - 1 )
subdivlvl = 0

def main():
    print('hallo')
    gamma_stab = 10
    nref = 2

    params =    {  
                    "order"     : order,
                    "mu"        : mu,
                    "gamma_stab": gamma_stab,
                    "kappa"     : "harmonic" # hansbo, highorder
                }


    # levelset fe space
    VhL = H1(mesh, order = order, dirichlet=[], dgjumps = True)
    lsetp1 = GridFunction(VhL)
    InterpolateToP1(levelset,lsetp1)
    ci = CutInfo(mesh, lsetp1)   
    

    Vh = H1(mesh, order=order, dirichlet=[1,2,3,4])
    VhNeg = Compress( Vh, active_dofs = GetActiveDof(mesh, HASNEG) )
    VhPos = Compress( Vh, active_dofs = GetActiveDof(mesh, HASPOS) )

    CutVh = FESpace( [VhNeg, VhPos], flags={"dgjumps": True} )

    # diffplot = GridFunction(CutVh)

    # ifdofs = IfaceUnks(CutVh,ci)

    # for i in range(len(diffplot.vec)):
    #     if ( ifdofs[i] ):
    #         diffplot.vec[i] = 1
    # ucoef = IfPos( lsetp1, diffplot.components[1], diffplot.components[0] )
    # Draw(ucoef, mesh, "diff")
    # print(ifdofs)
    # ifdofs &= CutVh.FreeDofs()
    # print(ifdofs)


    # return

    prolongNeg = P1Prolongation(mesh)
    prolongPos = P1Prolongation(mesh)
    CutProl = CompoundProlongation( CutVh )
    CutProl.AddProlongation(prolongNeg)
    CutProl.AddProlongation(prolongPos)
    prolongNeg.Update(VhNeg)
    prolongPos.Update(VhPos)
    

    a,f = AssembleCutPoisson( CutVh=CutVh, ci=ci, lsetp1=lsetp1, params=params )

    preJ = a.mat.CreateSmoother(CutVh.FreeDofs())
    GSpre = GaussSeid(preJ)

    #usol = solvers.PreconditionedRichardson(a,f.vec,pre=GSpre,maxit=500,tol=1e-6)

    #list of coarse grid matrices
    mats = {}
    mats[0] = a.mat
    #list of smoothers
    smoothers = {}
    smoothers[0] = CutFemSmoother(a, CutVh, ci) #preJ#GSpre
    #coarse grid solver
    inv = a.mat.Inverse(CutVh.FreeDofs(), inverse="sparsecholesky")
    
    for i in range(nref):
        mesh.Refine()
        #levelset update
        VhL.Update()
        lsetp1.Update()
        InterpolateToP1(levelset,lsetp1)
        ci = CutInfo(mesh,lsetp1)
        #cut space update
        Vh.Update()
        CutVh.components[0].SetActiveDofs( GetActiveDof(mesh,HASNEG) )
        CutVh.components[1].SetActiveDofs( GetActiveDof(mesh,HASPOS) )
        CutVh.Update()
        CutProl.Update(CutVh)
        #reassemble
        a,f = AssembleCutPoisson( CutVh=CutVh, ci=ci, lsetp1=lsetp1, params=params )        
        mats[i+1] = a.mat

        #preJ = a.mat.CreateSmoother(CutVh.FreeDofs())
        #GSpre = GaussSeid(preJ)
        smoothers[i+1] = CutFemSmoother(a, CutVh, ci) #preJ#GSpre        

    errors = [f.vec.CreateVector() for i in range(nref+1)]
    defects = [f.vec.CreateVector() for i in range(nref+1)]

    def MultiGridIter(rhs,nu,u,level):
        error = errors[level]
        defect = defects[level]

        if level == 0 :
            u.data = inv*rhs
            return
    #smoothing
        for i in range(nu):
            smoothers[level].Smooth(u,rhs)
    #compute defect
        defect.data = rhs - mats[level]*u
    #restrict defect
        CutProl.Restrict(level,defect)

        error[:] = 0
        MultiGridIter(defect,nu,error,level-1)

        CutProl.Prolongate(level,error)

        u.data += error
    #smoothing
        for i in range(nu):
            smoothers[level].Smooth(u,rhs)
        return u

    class MultiGridCL(BaseMatrix):
        def __init__(self,nu,l):
            super(MultiGridCL,self).__init__()
            self.nu = nu
            self.l = l
        def Mult(self,w,z):
            #z[:] = 0.0
            #z.data = w
            #z.data = 
            MultiGridIter(w,self.nu,z,self.l)
        def Height(self):
            return a.mat.height
        def Width(self):
            return a.mat.height

    MGpre = MultiGridCL(2,nref)

    usol = GridFunction(CutVh) 
    usol.vec[:] = 0.
    #set boundary condiitions
    #usol.components[1].Set(solution[1], BND)

    #usol.vec.data = mats[0] * usol.vec


    res = usol.vec.CreateVector()
    projres = usol.vec.CreateVector()
    tol = 1e-6
    proj = Projector(mask=CutVh.FreeDofs(),range=True) #projectors[nref]
    fnew = f.vec.data     

    projres.data = proj*fnew #f.vec
    normf = Norm(projres)

    for it in range(1,20):
        usol.vec.data = MGpre*fnew
        res.data = fnew - a.mat*usol.vec
        projres.data = proj*res
        res_norm = Norm(projres) / normf
        print("it =", it, " ||res||_2 =", res_norm)
        if res_norm < tol:
            break

    udraw = GridFunction( CutVh )
    udraw.components[1].Set( solution[1], BND )
    udraw.vec.data += usol.vec

    # Computation of L2 error:
    err_sqr_coefs = [(udraw.components[i]-solution[i])*(udraw.components[i]-solution[i]) for i in [0,1] ]
    lset_doms = LsetDoms( levelset, lsetp1, subdivlvl )
    l2error = sqrt( sum( [Integrate(lset_doms[i], cf=err_sqr_coefs[i], mesh=mesh, order=2*order, heapsize=1000000) for i in [0,1] ]))   
    print("L2 error : ",l2error)

    # pretty printing    
    ucoef = IfPos( lsetp1, udraw.components[1], udraw.components[0] )    
    #ucoef = IfPos( lsetp1, solution[1], solution[0] )    
    Draw (ucoef,mesh,'u')
    

class GaussSeid(BaseMatrix):
    def __init__ (self, smoother):
        super(GaussSeid, self).__init__()
        self.smoother = smoother
    def Mult (self, x, y):
        self.smoother.Smooth(y, x)
        #self.smoother.SmoothBack(y,x)
    def Height (self):
        return self.smoother.height
    def Width (self):
        return self.smoother.height

def GetActiveDof(mesh,HASPOSNEG):
        Vh = H1(mesh, order = order, dirichlet=[], dgjumps = True)
        lsetp1 = GridFunction(Vh)
        InterpolateToP1(levelset,lsetp1)
        ci = CutInfo(mesh,lsetp1)
        return GetDofsOfElements(Vh,ci.GetElementsOfType(HASPOSNEG))


def CutFESpace( V, ci, flags=None ):    
    hasneg = ci.GetElementsOfType(HASNEG)  
    haspos = ci.GetElementsOfType(HASPOS)

    if isinstance(V,list):
        if len(V) != 2:
            raise ValueError("you need to provide a list of two FE spaces")
        else:
            return FESpace( [   Compress(V[0],active_dofs=GetDofsOfElements(V[0],hasneg) ), 
                                Compress(V[1],active_dofs=GetDofsOfElements(V[1],haspos) ) 
                            ], flags=flags
                          ) 
    else:
        return FESpace( [   Compress(V,active_dofs=GetDofsOfElements(V,hasneg) ), 
                            Compress(V,active_dofs=GetDofsOfElements(V,haspos) ) 
                        ], flags=flags
                      ) 

def LsetDoms( levelset, lsetp1, subdivlvl ):
    if (subdivlvl == 0):
        lsetint = lsetp1
    else:
        lsetint = levelset

    lset_neg = { "levelset" : lsetint, "domain_type" : NEG, "subdivlvl" : subdivlvl}
    lset_pos = { "levelset" : lsetint, "domain_type" : POS, "subdivlvl" : subdivlvl}
    lset_if  = { "levelset" : lsetint, "domain_type" : IF , "subdivlvl" : subdivlvl}

    lset_dom = [ lset_neg, lset_pos,  lset_if ]

    return lset_dom

def KappaLambda(params,ci):
    kappac  = params["kappa"]
    mu      = params["mu"]
    order   = params["order"]

    if kappac == "harmonic":
        kappa = [ mu[1] / ( mu[0] + mu[1] ) , mu[0] /  ( mu[0] + mu[1] ) ]
        # stabilization parameter for Nitsche
        lambda_nitsche  = 2*mu[0]*mu[1]/(mu[0] + mu[1]) * 10 * order * order
    elif kappac == "hansbo":
        kappa = [CutRatioGF(ci),1.0-CutRatioGF(ci)]
        # stabilization parameter for Nitsche
        lambda_nitsche  = 0.5 * (mu[0]+mu[1]) * 20 * (order+1) * order
    elif kappac == "highorder":
        kappaHO = IfPos( CutRatioGF(ci) -0.5, 1, 0 )
        kappa = [kappaHO,1.0-kappaHO]
        # stabilization parameter for Nitsche
        lambda_nitsche  = 0.5 * (mu[0]+mu[1]) * 20 * order * order

    return [ kappa, lambda_nitsche ]

def AssembleCutPoisson(**kwargs):
    ci = kwargs['ci']
    hasneg      = ci.GetElementsOfType(HASNEG)  # <- "hasneg": has (also) negative level set values
    haspos      = ci.GetElementsOfType(HASPOS)  # <- "haspos": has (also) positive level set values    
    hasif       =  BitArray(ci.GetElementsOfType(IF))
    ba_facets   = [GetFacetsWithNeighborTypes(mesh,a=hasneg,b=hasif),GetFacetsWithNeighborTypes(mesh,a=haspos,b=hasif)]


    lsetp1 = kwargs['lsetp1']
    n_lset  = 1.0/grad(lsetp1).Norm() * grad(lsetp1)
    h       = specialcf.mesh_size
    params = kwargs['params']
    [kappa, lambda_nitsche ] = KappaLambda(params,ci)    

    # expressions of test and trial functions (u and v are tuples):
    CutVh = kwargs['CutVh']
    (u,v) = CutVh.TnT()

    gradu = [grad(ui) for ui in u]
    gradv = [grad(vi) for vi in v]

    mu = params["mu"]
    average_flux_u = sum([- kappa[i] * mu[i] * gradu[i] * n_lset for i in [0,1]])
    average_flux_v = sum([- kappa[i] * mu[i] * gradv[i] * n_lset for i in [0,1]])

    # integration domains
    lset_doms = LsetDoms( levelset, lsetp1, subdivlvl )

    #coef_f = [ x*y, x*y]

    # bilinear forms:

    a = BilinearForm(CutVh, symmetric = False)
    f = LinearForm(CutVh)

    #domain integrals:
    for i in [0,1] :
        a += SymbolicBFI( lset_doms[i], form = mu[i] * gradu[i] * gradv[i] )        

        #rhs
        f += SymbolicLFI( lset_doms[i], form = coef_f[i] * v[i] )

    # Nitsche integrals:
    a += SymbolicBFI(lset_doms[2], form = average_flux_u * (v[0]-v[1])
                                        + average_flux_v * (u[0]-u[1])
                                        + lambda_nitsche / h * (u[0]-u[1]) * (v[0]-v[1]) ) 

    for i in [0,1]:
        # Ghost penalty
        a += SymbolicFacetPatchBFI( form =  params["gamma_stab"]  * 
                                    mu[i] / h / h * 
                                    (u[i] - u[i].Other()) * (v[i] - v[i].Other()), 
                                    skeleton=False, definedonelements=ba_facets[i] )

    # # apply mesh adaptation 
    # if( doDeform ):   
    #     mesh.SetDeformation(deformation)

    # setting up matrix and vector
    a.Assemble()
    f.Assemble()

    #updating with non-homogeneous boundary conditions
    gfu = GridFunction( CutVh )
    #boundary lies in positive part
    gfu.components[1].Set( solution[1], BND )
    f.vec.data -= a.mat * gfu.vec

    return (a,f)

def IfaceUnks(CutVh, ci):
    ifdofs = BitArray(CutVh.ndof)
    ifdofs[:] = False
    hasif = BitArray( ci.GetElementsOfType(IF) )

    for el in CutVh.Elements():
        if( hasif[el.nr] ):
            verts = el.vertices
            for v in verts:
                for unk in CutVh.GetDofNrs(v):                    
                    ifdofs[unk] = True

    ifdofs &= CutVh.FreeDofs()
    return ifdofs

class CutFemSmoother:     
    def __init__(self,a,CutVh,ci):
        self.ifdofs = IfaceUnks(CutVh,ci)
        self.a = a
        self.inv = a.mat.Inverse(self.ifdofs, inverse="sparsecholesky")
        self.preJ = a.mat.CreateSmoother( CutVh.FreeDofs() )
    def Smooth( self, u,rhs ):
        self.preJ.Smooth(u, rhs)
        update = u.CreateVector()
        # this needs to be more efficient 
        # with help of ifdofs ...
        update.data = self.a.mat * u - rhs
        u.data -= self.inv * update




if __name__ == '__main__':
    main()   