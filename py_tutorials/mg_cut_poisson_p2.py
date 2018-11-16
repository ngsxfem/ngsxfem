# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
from ngsolve.krylovspace import CG
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from xfem import *
from xfem.lsetcurv import *


ngsglobals.msg_level = 1

# generate a triangular mesh of mesh-size 0.2
square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.2, quad_dominated=False))

mu = [1e-5, 1]
order = 1
kappac="harmonic"
gamma_stab = 10


rhsscal = -4*mu[0]*mu[1] 
coef_f = [ rhsscal, rhsscal ]
distsq = x*x + y*y - 1
solution = [ mu[1]* distsq, mu[0] * distsq ]

levelset = ( sqrt(x*x + y*y) - 1 )
subdivlvl = 0

def main():
    print('hallo')
    
    nref = 2

    params =    {  
                    "order"     : order,
                    "mu"        : mu,
                    "gamma_stab": gamma_stab,
                    "kappa"     : "harmonic" # hansbo, highorder
                }

    mgiter = LinearMGIterator(nref=nref,nu=2)

    rhs = mgiter.getRhs()

    udraw = mgiter.iterate(rhs)

    lsetp1 = mgiter.getLset()

    # Computation of L2 error:
    err_sqr_coefs = [(udraw.components[i]-solution[i])*(udraw.components[i]-solution[i]) for i in [0,1] ]
    lset_doms = LsetDoms( levelset, lsetp1, subdivlvl )
    # input('w2')
    #l2error = sqrt( sum( [Integrate(lset_doms[i], cf=err_sqr_coefs[i], mesh=mesh, order=2*order, heapsize=1000000) for i in [0,1] ]))   
    l2error = sqrt( sum( [Integrate(lset_doms[i], cf=err_sqr_coefs[i], mesh=mesh) for i in [0,1] ]))   
    print("L2 error : ",l2error)

    input('draw?\n')

    # if( nref > 2):
    #     print('nref too large')
    #     return
    #     print('test')

    print('drawing')
    # pretty printing    
    ucoef = IfPos( lsetp1, udraw.components[1], udraw.components[0] )
    Draw (ucoef,mesh,'u')

    input('finish? - seg fault\n')
    

def GetActiveDof(mesh,HASPOSNEG):
        Vh = H1(mesh, order = order, dirichlet=[], dgjumps = True)
        lsetp1 = GridFunction(Vh)
        InterpolateToP1(levelset,lsetp1)
        ci = CutInfo(mesh,lsetp1)
        return GetDofsOfElements(Vh,ci.GetElementsOfType(HASPOSNEG))


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

def KappaLambda(ci,order):
    #kappac  = params["kappa"]
    #mu      = params["mu"]
    #order   = params["order"]

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
    #params = kwargs['params']
    if( kwargs['highorder'] ):
        kappaorder = order
    else:
        kappaorder = 1

    [kappa, lambda_nitsche ] = KappaLambda(ci, kappaorder)    

    # expressions of test and trial functions (u and v are tuples):
    CutVh = kwargs['CutVh']
    (u,v) = CutVh.TnT()

    gradu = [grad(ui) for ui in u]
    gradv = [grad(vi) for vi in v]

    average_flux_u = sum([- kappa[i] * mu[i] * gradu[i] * n_lset for i in [0,1]])
    average_flux_v = sum([- kappa[i] * mu[i] * gradv[i] * n_lset for i in [0,1]])

    # integration domains
    lset_doms = LsetDoms( levelset, lsetp1, subdivlvl )
    
    if( kwargs['highorder'] ):
        mesh.SetDeformation( kwargs['deformation'] )


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
        a += SymbolicFacetPatchBFI( form =  gamma_stab  * 
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
                # print( CutVh.GetDofNrs(v) )
                # print( CutVh.components[0].GetDofNrs(v) )
                # print( CutVh.components[1].GetDofNrs(v) )
                # print( CutVh.components[0].ndof )
                # print('-------------------------------')
                for unk in CutVh.GetDofNrs(v):                    
                    ifdofs[unk] = True

    ifdofs &= CutVh.FreeDofs()
    #print(ifdofs)
    #input('weiter')
    return ifdofs

class CutFemSmoother:     
    def __init__(self,a,CutVh,ci):        
        self.a = a
        ifdofs = IfaceUnks(CutVh,ci)
        #self.proj = Projector(mask=self.ifdofs,range=True)
        ifdofslst = [ [i] for i in range(len(ifdofs)) if ifdofs[i]==True ]
        self.proj = a.mat.CreateBlockSmoother(ifdofslst)
        #self.inv = a.mat.Inverse(self.ifdofs, inverse="sparsecholesky")
        self.preJ = a.mat.CreateSmoother( CutVh.FreeDofs() )
    def __del__(self):
            print('smoother dying')
            input('destruktor')
    def Smooth( self, u,rhs ):
        self.preJ.Smooth(u, rhs)
        #return
        update = u.CreateVector()
        # this needs to be more efficient 
        # with help of ifdofs ...
        update.data = self.a.mat * u - rhs
        cgoutput = u.CreateVector()
        cgoutput.data = CG(self.a.mat, update, pre=self.proj, tol=1e-2,printrates=False)
        #u.data -= self.inv * update
        u.data -= cgoutput

class MultiGridCL(BaseMatrix):
        def __init__(self,**kwargs):
            super(MultiGridCL,self).__init__()
            self.nu     = kwargs.get('nu',2)            
            self.l      = kwargs['level']
            self.prol   = kwargs['prol']
            self.mats   = kwargs['matrices']
            self.smoothers  = kwargs['smoothers']
            self.inv    = kwargs["coarsegridsolver"]
            f = kwargs["f"]
            nref = self.l
            self.errors     = [f.vec.CreateVector() for i in range(nref+1)]
            self.defects    = [f.vec.CreateVector() for i in range(nref+1)]
        
        def __del__(self):
            print('mg dying')
            input('destruktor')


        def MultiGridIter(self,rhs,u,level):
            error   = self.errors[level]
            defect  = self.defects[level]

            if level == 0 :
                u.data = self.inv*rhs
                return
        #smoothing
            for i in range(self.nu):
                self.smoothers[level].Smooth(u,rhs)
        #compute defect
            defect.data = rhs - self.mats[level]*u
        #restrict defect
            self.prol.Restrict(level,defect)
            error[:] = 0
            #v-cycle
            self.MultiGridIter(defect,error,level-1)
            self.prol.Prolongate(level,error)
            u.data += error
        #smoothing
            for i in range(self.nu):
                self.smoothers[level].Smooth(u,rhs)
            return u

        def Mult(self,w,z):
            #z[:] = 0.0
            #z.data = w
            #z.data = 
            self.MultiGridIter(w,z,self.l)
        def Height(self):
            return a.mat.height
        def Width(self):
            return a.mat.height

def CreateLsetMeshAdapt(order):
    lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)
    deformation = lsetmeshadap.CalcDeformation(levelset)
    lsetp1 = lsetmeshadap.lset_p1

    return (lsetmeshadap, deformation, lsetp1 )

class LinearMGIterator:
    def __init__(self, **kwargs):
        self.nref = kwargs['nref']

        # levelset fe space
        VhL = H1(mesh, order = 1, dirichlet=[], dgjumps = True)
        lsetp1 = GridFunction(VhL)
        InterpolateToP1(levelset,lsetp1)

        ci = CutInfo(mesh, lsetp1)   
        
        # cut fe space
        Vh = H1(mesh, order=1, dirichlet=[1,2,3,4])
        VhNeg = Compress( Vh, active_dofs = GetActiveDof(mesh, HASNEG) )
        VhPos = Compress( Vh, active_dofs = GetActiveDof(mesh, HASPOS) )

        CutVh = FESpace( [VhNeg, VhPos], flags={"dgjumps": True} )

        self.CutVh = CutVh

        # construct prolongation for neg/pos space
        prolongNeg = P1Prolongation(mesh)
        prolongPos = P1Prolongation(mesh)
        # create cut prolongation
        CutProl = CompoundProlongation( CutVh )
        CutProl.AddProlongation(prolongNeg)
        CutProl.AddProlongation(prolongPos)
        prolongNeg.Update(VhNeg)
        prolongPos.Update(VhPos)
        

        #create coarse level bilinear/linear form
        self.a,self.f = AssembleCutPoisson( CutVh=CutVh, ci=ci, lsetp1=lsetp1, highorder=False )

        
        #list of coarse grid matrices
        self.mats = {}
        self.mats[0] = self.a.mat
        #list of smoothers
        self.smoothers = {}
        self.smoothers[0] = CutFemSmoother(self.a, CutVh, ci) 
        #coarse grid solver
        inv = self.a.mat.Inverse(CutVh.FreeDofs(), inverse="sparsecholesky")
        

        # create multilevel hierarchy
        for i in range(self.nref):
            print('preforming refinement: ', i+1)
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
            self.a,self.f = AssembleCutPoisson( CutVh=CutVh, ci=ci, lsetp1=lsetp1, highorder=False )        
            self.mats[i+1] = self.a.mat

            self.smoothers[i+1] = CutFemSmoother(self.a, CutVh, ci)

        
        self.lsetp1 = lsetp1
        # create multigrid solver
        self.MGpre = MultiGridCL ( level=self.nref, prol=CutProl, 
                            matrices=self.mats, smoothers=self.smoothers,
                            coarsegridsolver=inv, f=self.f 
                            )


    def getRhs(self):
        return self.f.vec.data

    def createVec(self):
        return self.f.vec.CreateVector()

    def getLset(self):
        return self.lsetp1                

    def iterate(self,rhs,tol=1e-6, maxit=10,printinfo=True):
        usol = GridFunction(self.CutVh) 
        usol.vec[:] = 0.

        res = usol.vec.CreateVector()
        projres = usol.vec.CreateVector()
        
        proj = Projector(mask=self.CutVh.FreeDofs(),range=True) 
        fnew = rhs #f.vec.data     

        projres.data = proj*fnew #f.vec
        normf = Norm(projres)

        for it in range(maxit):
            usol.vec.data = self.MGpre*fnew
            res.data = fnew - self.a.mat*usol.vec
            projres.data = proj*res
            res_norm = Norm(projres) / normf
            print("it =", it+1, " ||res||_2 =", res_norm)
            if res_norm < tol:
                break        

        udraw = GridFunction( self.CutVh )
        udraw.components[1].Set( solution[1], BND )
        udraw.vec.data += usol.vec

        return udraw    

if __name__ == '__main__':
    main()   