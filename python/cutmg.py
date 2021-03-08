from ngsolve.la import BaseMatrix
from ngsolve import Projector, Norm
from ngsolve.krylovspace import CG
from xfem import *
from ngsolve import *
import ngsolve
'''
Below we implement a generic MultiGrid class based on the following (yet to be defined) objects:
* prolongation : prolongation from last to current mesh level 
* smoothers : smoothers corresponding to each mesh level
* coarsegridsolver : BaseMatrix for the coarse(st) level representing an (approx.) inverse
* mats : matrices corresponding to each mesh level
* nu : number of smoothing steps

The MultiGrid class is a BaseMatrix. Hence, it essentially implements the Mult operation.
'''
class MultiGridCL(BaseMatrix):
	def __init__(self,**kwargs):
	    super(MultiGridCL,self).__init__()
	    self.nu     = kwargs.get('nu',2)            
	    self.l      = kwargs['level'] # maximum number of level
	    self.prol   = kwargs['prol']
	    self.mats   = kwargs['matrices']
	    self.smoothers  = kwargs['smoothers']
	    self.inv    = kwargs["coarsegridsolver"]
	    self.errors     = [self.mats[-1].CreateColVector() for i in range(self.l+1)]
	    self.defects    = [self.mats[-1].CreateColVector() for i in range(self.l+1)]
	    self.level = self.l # current level
	
	def Mult(self,rhs,u):            
	    error   = self.errors[self.level]
	    defect  = self.defects[self.level]
	
	    if self.level == 0 :
	        u.data = self.inv*rhs
	        return
	#smoothing
	    self.smoothers[self.level].Smooth(u,rhs,self.nu)
	#compute defect
	    defect.data = rhs - self.mats[self.level]*u
	#restrict defect
	    self.prol.Restrict(self.level,defect)
	    error[:] = 0
	    #v-cycle
	    self.level -= 1
	    self.Mult(defect,error)
	    self.level += 1
	    self.prol.Prolongate(self.level,error)
	    u.data += error
	#smoothing
	    self.smoothers[self.level].SmoothBack(u,rhs,self.nu)
	    return u
	
	def Height(self):
	    return self.mats[-1].height
	def Width(self):
	    return self.mats[-1].width



class CutFemSmoother:     
    def __init__(self,**kwargs):
        CutVh, ci, self.a  = kwargs["CutVh"],kwargs["ci"], kwargs["a"]
        self.ifsolver = kwargs["ifsolver"]
        ifdofs = GetDofsOfElements(CutVh, ci.GetElementsOfType(IF)) & CutVh.FreeDofs()
        #self.proj = Projector(mask=self.ifdofs,range=True)
        ifdofslst = [ [i] for i in range(len(ifdofs)) if ifdofs[i]==True ]
        if self.ifsolver == "direct":
            self.ifpre = self.a.mat.Inverse(ifdofs, inverse="sparsecholesky")
        elif self.ifsolver == "cg":
            print("cg if solver")
            self.ifpre = self.a.mat.CreateBlockSmoother(ifdofslst)
        else:
            pass
            
        if "blocks" in kwargs:
            self.preJ = self.a.mat.CreateBlockSmoother( kwargs["blocks"] )
        else:
            self.preJ = self.a.mat.CreateSmoother( CutVh.FreeDofs() )
        
        self.ifcorr_only_once = kwargs.get('ifcorr_only_once',True)

    def InterfaceCorrection(self, u, rhs,k,nu):
        if self.ifsolver == None:
            return
        if self.ifcorr_only_once and k < nu-1:
            return
        update = u.CreateVector()
        # this needs to be more efficient 
        # with help of ifdofs ...
        update.data = self.a.mat * u - rhs
        cgoutput = u.CreateVector()

        E = Embedding(len(rhs), IntRange(0,self.a.mat.height))
        
        cgoutput.data = CG(E @ self.a.mat @ E.T,
                           update,
                           pre=E @ self.ifpre @ E.T,
                           tol=1e-2,printrates=False)
        #u.data -= self.inv * update
        u.data -= cgoutput

    def Smooth( self, u, rhs, nu ):
        for k in range(nu):
            self.preJ.Smooth(u, rhs)        
            self.InterfaceCorrection(u,rhs,k,nu)
            
    def SmoothBack( self, u, rhs, nu ):
        for k in reversed(range(nu)):
            self.InterfaceCorrection(u,rhs,k,nu)
            self.preJ.SmoothBack(u, rhs)



class LinearMGIterator(BaseMatrix):
    def __init__(self, **kwargs):
        super(LinearMGIterator,self).__init__()
        #self.nref = kwargs['nref']
        self.ifsolver = kwargs['ifsolver']
        #AssembleProblem = kwargs['assemblefct'] #not used
        ProlType = kwargs.get('ProlType',P1Prolongation)
        self.nu = kwargs.get('nu',2)
        mesh = kwargs['mesh']
        
        self.tol=kwargs.get("tol",1e-6)
        self.maxit=kwargs.get("maxit",20)
        self.printinfo=kwargs.get("printinfo",True)
        
        ci, lsetp1 = kwargs["ci"], kwargs["lsetp1"]
        #create coarse level bilinear/linear form
        a = kwargs["a"]
        self.CutVh = a.space
        
        # levelset fe space
        #VhL = lsetp1.space

        ci = CutInfo(mesh, lsetp1)
        
        print( type(self.CutVh) )
        #check for compound space
        if isinstance(self.CutVh,ProductSpace):
            print('creating compound cut prolongation')
            self.CutProl = CompoundProlongation( self.CutVh )
            # construct prolongation for neg/pos space
            self.CutProl.AddProlongation( ProlType(mesh) )
            self.CutProl.AddProlongation( ProlType(mesh) )
        else:
            print('creating single phase cut prolongation')
            self.CutProl = ProlType(mesh)
        
        self.CutProl.Update(self.CutVh)

        #self.mats = {}
        self.mats = [a.mat]
        #list of smoothers
        self.smoothers = [None]
        #coarse grid solver
        self.coarseinv = a.mat.Inverse(self.CutVh.FreeDofs(), inverse="sparsecholesky")
        self.MGpre = self.coarseinv
    
    def Update(self, a, CutVh, ci):
        self.CutProl.Update(CutVh)
        self.mats.append(a.mat)
        self.CutVh = CutVh
        current_smoother = CutFemSmoother(a=a, CutVh=CutVh, ci=ci, ifsolver=self.ifsolver,
                                                 ifcorr_only_once = True)
        self.smoothers.append(current_smoother)
        # create/update multigrid solver
        self.MGpre = MultiGridCL ( level=len(self.mats)-1, prol=self.CutProl, 
                                   matrices=self.mats, smoothers=self.smoothers,
                                   nu=self.nu,
                                   coarsegridsolver=self.coarseinv)        
        
    def createVec(self):
        return self.mats[-1].CreateColVector()

    def getSpaceDim(self):
        print('-----------------------')
        print(type(self.CutVh))
        #if( type(self.CutVh) == ngsolve.comp.CompoundFESpace ):
        if( type(self.CutVh) == ngsolve.comp.CompoundFESpace or type(self.CutVh) == ngsolve.comp.FESpace):
            numspaces = len ( self.CutVh.components )
            return [ self.CutVh.components[i].ndof for i in range(numspaces) ]
        else:
            return [self.CutVh.ndof]
    
    
    def Mult(self,rhs,usol):            
        usol[:] = 0.

        res = usol.CreateVector()
        projres = usol.CreateVector()
        
        proj = Projector(mask=self.CutVh.FreeDofs(),range=True) 
        fnew = rhs #f.vec.data     

        projres.data = proj*fnew #f.vec
        normf = Norm(projres)

        for it in range(self.maxit):
            usol.data = self.MGpre*fnew
            res.data = fnew - self.mats[-1]*usol
            projres.data = proj*res
            res_norm = Norm(projres) / normf
            if self.printinfo:
                print("it =", it+1, " ||res||_2 =", res_norm)
            if res_norm < self.tol:
                break        

        #udraw = GridFunction( self.CutVh )
        #udraw.components[1].Set( solution[1], BND )
        #udraw.vec.data += usol.vec

        return usol

    def Height(self):
        return self.mats[-1].height
    def Width(self):
        return self.mats[-1].height

def VertPatches(fes,mesh):
    blocks = []
    freedofs = fes.FreeDofs()
    for v in mesh.vertices:
        vdofs = set()
        for el in mesh[v].elements:
            #print( fes.GetDofNrs(el) )
            #input('w')
            vdofs |= set(d for d in fes.GetDofNrs(el) if (freedofs[abs(d)] and d >=0) )
        #print(vdofs)
        blocks.append (vdofs)
    return blocks

def ElemPatches(fes):
    blocks = []
    freedofs = fes.FreeDofs()
    for el in fes.Elements():
        eldofs = set(d for d in fes.GetDofNrs(el) if (freedofs[abs(d)] and d >=0) )        
        blocks.append (eldofs)
    return blocks

def EdgePatches(fes,mesh):
    blocks = []
    freedofs = fes.FreeDofs()
    for edge in mesh.edges:
        edofs = fes.GetDofNrs(edge)

        if( len(edofs) == 0 ): continue
        
        # print('found edge dof: ', edofs)
        elcntr = 0
        edofs = set()
        for el in mesh[edge].elements:
            elcntr += 1
            # print('  with element ', elcntr, ' el: ', fes.GetDofNrs(el) )
            # print('  with element ', elcntr, ' el: ', el.nr )
            edofs |= set( d for d in fes.GetDofNrs(el) if(freedofs[abs(d)] and d>=0) )
    
        # input('next edge')
        blocks.append (edofs)
    return blocks


class P2TwoGridCL(BaseMatrix):
    def __init__(self,**kwargs):
        super(P2TwoGridCL,self).__init__()
        self.a = kwargs['a']
        #self.f = kwargs['f']
        self.fes = kwargs['fes']
        
        self.patchtype = kwargs.get('patchtype','edge')

        mesh = kwargs['mesh']

        if(self.patchtype == 'edge'):
            blocks = EdgePatches(self.fes,mesh)
        elif self.patchtype == 'vert':
            blocks = VertPatches(self.fes,mesh)
        elif self.patchtype == 'elem':
            blocks = ElemPatches(self.fes)
        else:
            raise ValueError("Unknown patchtype for block smoothing!\n choose between 'edge', 'vert' and 'elem'")
        
        self.ifsolver = kwargs['ifsolver']
        self.blockjac = CutFemSmoother( a=self.a, CutVh=self.fes, 
                                        ci=kwargs["ci"], blocks=blocks, ifsolver=self.ifsolver,
                                        ifcorr_only_once=False )

        # we need a linear multigrid as coarse grid correction
        self.linmgiter = kwargs['linmgiter']
        self.maxit = kwargs.get('maxit',20)
        self.tol = kwargs.get('tol',1e-6)
        self.nu = kwargs.get( 'nu',3 )

    def Mult(self,rhs,usol):
        proj = Projector(mask=self.fes.FreeDofs(),range=True)
        usol[:] = 0.
        #
        projres = rhs.CreateVector()
        projres.data = proj*rhs
        normf = Norm(projres)
        #
        LinDof = self.linmgiter.getSpaceDim()
        
        if( type(self.fes) == ngsolve.comp.CompoundFESpace or type(self.fes) == ngsolve.comp.FESpace):
            numspaces = len( self.fes.components )
            HOdof = [ self.fes.components[i].ndof for i in range(numspaces) ]
        else:
            HOdof = [ self.fes.ndof ]
        #
        print('lindof, hdof: ')
        print(LinDof,HOdof)
        oldres = Norm(projres)

        for it in range(self.maxit):
            # initial smoothing 
            self.blockjac.Smooth(usol, rhs, self.nu)

            res = usol.CreateVector()
            # perform multiplication more efficiently ...
            # only 'upper' part (vertex dofs) needed 
            res.data = rhs -self.a.mat * usol 
            #eliminate boundary condition
            projres.data = proj*res           

            mgrhs = self.linmgiter.createVec()            
            # copy data from two domains 
            # ... should be done in combination with vertexdofs ...
            #mgrhs.Range(0, LinDof[0]).data = projres.Range(0, LinDof[0] )            
            #mgrhs.Range( LinDof[0], len(mgrhs) ).data = projres.Range(HOdof[0], HOdof[0] + LinDof[1] )
            
            mgrhs.Range(0, LinDof[0]).data = projres.Range(0, LinDof[0] )
            if( len(LinDof) > 1 ):
                offsetc = 0
                offsetf = 0
                for i in range( len(LinDof) -1 ):
                    offsetc += LinDof[i]
                    offsetf += HOdof[i]
                    mgrhs.Range( offsetc, offsetc + LinDof[i+1] ).data = projres.Range(offsetf, offsetf + LinDof[i+1] )
            
            # coarse grid correction: multigrid on linears
            cup = mgrhs.CreateVector()
            cup.data = self.linmgiter*mgrhs #maxit=1,printinfo=False)

            # update solution
            # ... should be done in combination with vertexdofs ...
            
            usol.Range(0, LinDof[0] ).data += cup.Range(0, LinDof[0])
            if( len(LinDof) > 1 ):
                offsetc = 0
                offsetf = 0
                for i in range( len(LinDof) -1 ):
                    offsetc += LinDof[i]
                    offsetf += HOdof[i]
                    usol.Range( offsetf, offsetf + LinDof[i+1] ).data += cup.Range( offsetc, offsetc + LinDof[i+1] )
            '''
            usol.Range(0, LinDof[0] ).data += cup.Range(0, LinDof[0])
            usol.Range( HOdof[0], HOdof[0] + LinDof[1]).data += cup.Range( LinDof[0], len(mgrhs) )
            '''

            # compute residual for measurement
            res.data = rhs - self.a.mat*usol
            projres.data = proj*res
            res_norm = Norm(projres)
            
            print("tg-it =", it+1, "\t ||res||_2 = {0:.2E}".format(res_norm), "\t reduction: {0:.2f}".format(res_norm/oldres))
            if res_norm < self.tol*normf:
                break
            oldres = res_norm
        
        return usol

    
    # def iterate(self,nu=2,maxit=20,tol=1e-6):
        
    #     proj = Projector(mask=self.fes.FreeDofs(),range=True)
    #     usol = GridFunction(self.fes)
    #     usol.vec[:] = 0.

    #     projres = usol.vec.CreateVector()
    #     projres.data = proj*self.f.vec
    #     normf = Norm(projres)

    #     LinDof = self.mgiter.getSpaceDim()
    #     HOdof  = [ self.fes.components[i].ndof for i in [0,1] ]

    #     oldres = Norm(projres)

    #     for it in range(maxit):
    #         # initial smoothing 
    #         #for sm in range(nu):
    #         self.blockjac.Smooth(usol.vec, self.f.vec,nu)

    #         res = usol.vec.CreateVector()
    #         # perform multiplication more efficiently ...
    #         # only 'upper' part (vertex dofs) needed 
    #         res.data = self.f.vec -self.a.mat * usol.vec
            

    #         mgrhs = self.mgiter.createVec()
    #         # copy data from two domains 
    #         # ... should be done in combination with vertexdofs ...
    #         mgrhs.Range(0, LinDof[0]).data = projres.Range(0, LinDof[0] )
    #         mgrhs.Range( LinDof[0], len(mgrhs) ).data = projres.Range(HOdof[0], HOdof[0] + LinDof[1] )
    #         #print('lens: ', len(mgrhs), LinDof[0], LinDof[1], LinDof[0]+LinDof[1])

    #         # coarse grid correction: multigrid on linears
    #         cup = self.mgiter.iterate(mgrhs,maxit=1,printinfo=False)

    #         # update solution
    #         # ... should be done in combination with vertexdofs ...
    #         usol.vec.Range(0, LinDof[0] ).data += cup.vec.Range(0, LinDof[0])
    #         usol.vec.Range( HOdof[0], HOdof[0] + LinDof[1]).data += cup.vec.Range( LinDof[0], len(mgrhs) )

    #         # compute residual for measurement
    #         res.data = self.f.vec - self.a.mat*usol.vec
    #         projres.data = proj*res
    #         res_norm = Norm(projres)
            
    #         print("tg-it =", it+1, "\t ||res||_2 = {0:.2E}".format(res_norm), "\t reduction: {0:.2f}".format(res_norm/oldres))
    #         if res_norm < tol*normf:
    #             break
    #         oldres = res_norm

    #     udraw = GridFunction( self.fes )
    #     udraw.components[1].Set( solution[1], BND )
    #     udraw.vec.data += usol.vec
    #     return udraw
