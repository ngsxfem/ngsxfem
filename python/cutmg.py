from ngsolve.la import BaseMatrix
from ngsolve import Projector, Norm
from ngsolve.krylovspace import CG
from xfem import *
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
            self.ifpre = a.mat.Inverse(ifdofs, inverse="sparsecholesky")
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
        
        cgoutput.data = CG(self.a.mat, update, pre=self.ifpre, tol=1e-2,printrates=False)
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
        self.nref = kwargs['nref']
        self.ifsolver = kwargs['ifsolver']
        AssembleProblem = kwargs['assemblefct']
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
        if( type(self.CutVh) == ngsolve.comp.CompoundFESpace ):
            print('creating compound cut prolongation')
            self.CutProl = CompoundProlongation( self.CutVh )
            # construct prolongation for neg/pos space
            self.CutProl.AddProlongation( P1Prolongation(mesh) )
            self.CutProl.AddProlongation( P1Prolongation(mesh) )
        else:
            print('creating single phase cut prolongation')
            self.CutProl = P1Prolongation(mesh)
        
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
                                   coarsegridsolver=self.coarseinv)        
        
    def createVec(self):
        return self.mats[-1].CreateColVector()

    def getSpaceDim(self):
        if( type(self.CutVh) == ngsolve.comp.CompoundFESpace ):
            return [ self.CutVh.components[i].ndof for i in [0,1] ]
        else:
            return self.CutVh.ndof
    
    
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


