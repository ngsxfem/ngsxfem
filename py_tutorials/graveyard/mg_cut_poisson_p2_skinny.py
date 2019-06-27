# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
from ngsolve.krylovspace import CG
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
import ngsolve
from xfem import *
from xfem.cutmg import * #MultiGridCL, CutFemSmoother, LinearMGIterator
from xfem.lsetcurv import *


ngsglobals.msg_level = 1

# generate a triangular mesh of mesh-size 0.2
square = SplineGeometry()
square.AddRectangle([-1.5,-1.5],[1.5,1.5],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=0.5, quad_dominated=False))

mu = [1e-3, 1]
order = 2
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

    for i in range(nref+1):
        if i==0:
            #lset
            lsetp1 = GridFunction(H1(mesh, order = 1))
            InterpolateToP1(levelset,lsetp1)
            ci = CutInfo(mesh, lsetp1)
            #cut space
            Vh = H1(mesh, order=1, dirichlet=[1,2,3,4])
            VhNeg = Compress( Vh, active_dofs = GetDofsOfElements(Vh,ci.GetElementsOfType(HASNEG)))
            VhPos = Compress( Vh, active_dofs = GetDofsOfElements(Vh,ci.GetElementsOfType(HASPOS)))
            CutVh = FESpace( [VhNeg, VhPos], flags={"dgjumps": True} )
        else:
            mesh.Refine()

            lsetp1.space.Update()
            lsetp1.Update()
            InterpolateToP1(levelset,lsetp1)
            ci = CutInfo(mesh,lsetp1)

            CutVh.components[0].GetBaseSpace().Update()
            CutVh.components[0].SetActiveDofs( GetDofsOfElements(Vh,ci.GetElementsOfType(HASNEG)))
            CutVh.components[1].SetActiveDofs( GetDofsOfElements(Vh,ci.GetElementsOfType(HASPOS)))
            CutVh.Update()
            

        a, f = AssembleCutPoisson( CutVh=CutVh, ci=ci, lsetp1=lsetp1, highorder=False )

        if i==0:
            mgiter = LinearMGIterator(a=a, ci=ci, lsetp1=lsetp1, mesh=mesh,
                                  nu=1,maxit=1, printinfo=False, 
                                  ifsolver="cg" )
        else:
            mgiter.Update( a, CutVh, ci )



    Vh = H1(mesh, order=order, dirichlet=[1,2,3,4])
    (lsetmeshadap, deformation, lsetp1 ) = CreateLsetMeshAdapt(order)
    ci = CutInfo(mesh,lsetp1)

    HoCutFes = CutFESpace( Vh, ci, flags={"dgjumps": True} )    

    a, f = AssembleCutPoisson( CutVh=HoCutFes, ci=ci, lsetp1=lsetp1, highorder=True, deformation=deformation )

    
    tgiter = P2TwoGridCL( a=a, linmgiter=mgiter, fes=HoCutFes, ci=ci, nu=3, maxit=50,
                            mesh=mesh, ifsolver="cg", patchtype='edge' )


    gfu = GridFunction( HoCutFes )
    gfu.components[1].Set( solution[1], BND )
    f.vec.data -= a.mat * gfu.vec
    update = gfu.vec.CreateVector()

    update = tgiter * f.vec

    gfu.vec.data += update

    #udraw = tgiter.iterate(nu=3,tol=1e-8) #udraw = mgiter.iterate(rhs)    

    # Computation of L2 error:
    err_sqr_coefs = [(gfu.components[i]-solution[i])*(gfu.components[i]-solution[i]) for i in [0,1] ]
    lset_doms = LsetDoms( levelset, lsetp1, subdivlvl )
    # input('w2')
    l2error = sqrt( sum( [Integrate(lset_doms[i], cf=err_sqr_coefs[i], mesh=mesh, order=2*order, heapsize=1000000) for i in [0,1] ]))   
    #l2error = sqrt( sum( [Integrate(lset_doms[i], cf=err_sqr_coefs[i], mesh=mesh) for i in [0,1] ]))   
    print("L2 error : ",l2error)

    input('draw?\n')

    # if( nref > 2):
    #     print('nref too large')
    #     return
    #     print('test')

    print('drawing')
    # pretty printing    
    ucoef = IfPos( lsetp1, gfu.components[1], gfu.components[0] )
    if order>1 :
        mesh.UnsetDeformation()
    Draw (ucoef,mesh,'u')



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

    if( kwargs['highorder'] ):
        print('setting deformation')
        mesh.SetDeformation( kwargs['deformation'] )

    # setting up matrix and vector
    a.Assemble()
    f.Assemble()

    #updating with non-homogeneous boundary conditions
    #gfu = GridFunction( CutVh )
    #boundary lies in positive part
    #gfu.components[1].Set( solution[1], BND )
    #f.vec.data -= a.mat * gfu.vec

    return (a,f)
    

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



def CreateLsetMeshAdapt(order):
    lsetmeshadap = LevelSetMeshAdaptation(mesh, order=order, threshold=1000, discontinuous_qn=True)
    deformation = lsetmeshadap.CalcDeformation(levelset)
    lsetp1 = lsetmeshadap.lset_p1

    return (lsetmeshadap, deformation, lsetp1 )


if __name__ == '__main__':
    main()   
