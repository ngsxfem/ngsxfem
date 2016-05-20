from ngsolve.fem import *
from ngsolve.comp import *
from ngsolve.solve import *
from ngsolve.ngstd import *
from ngsolve.la import *


from libngsxfem_py.xfem import *
from libngsxfem_xfem import *
from libngsxfem_osmosis import *
from math import *

from extension import *
from roundedsquare import *


def osmosis(dim,mesh,lset,initialconcentration,alpha,beta,dt,tend,theta=0.5,smoothextension=True, stabilization_interface = False, stabilization_extension = False):

  stabilization_at_all = stabilization_interface or stabilization_extension

  h1fes = FESpace("h1ho",mesh,order=3,dirichlet=[])

  gflset = GridFunction(h1fes,"levelset")
  gflset.Update()
  gflset.Set(lset)

  felset = GFCoeff(gflset)

  fes = FESpace("xstdfespace",mesh=mesh, flags={"type_std":"h1ho","dgjumps": stabilization_at_all,"empty":True,"ref_space":1}, order=1,dirichlet=[1,2,3,4])
  xfes = CastToXStdFESpace(fes)
  xfes.XFESpace.SetLevelSet(felset)
  fes.Update()

  gfu = GridFunction(fes,"gfu")
  gfu.Update()

  gfu.components[0].Set(ConstantCF(1))
  uone = gfu.vec.CreateVector()
  uone.data = gfu.vec
  utmp = gfu.vec.CreateVector()
  
  bfm = BilinearForm(h1fes,"bfm")
  bfm.Add(BFI("meancurv_mass",dim=dim,coef=[]))

  bfa = BilinearForm(h1fes,"bfa")
  bfa.Add(BFI("meancurv_stiff",dim=dim,coef=[]))

  bfau = BilinearForm(fes,"bfau")
  bfau.Add(BFI("xlaplace",dim=dim,coef=[ConstantCF(alpha),ConstantCF(0)]))
  if (stabilization_interface == True):
      bfau.Add(BFI("lo_ghostpenalty",dim=dim,coef=[ConstantCF(alpha),ConstantCF(0),ConstantCF(1)]))
      # for order 2:
      # bfau.Add(BFI("sec_ghostpenalty",coef=[ConstantCF(alpha),ConstantCF(0),ConstantCF(1)]))
      

  bfmu = BilinearForm(fes,"bfmu")
  bfmu.Add(BFI("xmass",dim=dim,coef=[ConstantCF(1),ConstantCF(0)]))


  bfau.Assemble()
  bfmu.Assemble()


  bfextend = BilinearForm(fes,"a",flags={"symmetric":True})
  bfextend.Add(CompoundBFI(bfi=BFI("laplace",dim=dim,coef=ConstantCF(1)), comp = 0))
  if (stabilization_extension == True):
      bfextend.Add(CompoundBFI(bfi=BFI("dudnjumpdvdnjump",dim=dim,coef=ConstantCF(1)), comp = 0))
      # stabilization for higher order is missing here: ...

  bfextend.Assemble()

  bfmass = BilinearForm(fes,"a",flags={"symmetric":True})
  bfmass.Add(CompoundBFI(bfi=BFI("mass",dim=dim,coef=ConstantCF(1.0/dt)), comp = 0))
  bfmass.Assemble()

  uold = gfu.vec.CreateVector()

  vtkout = VTKOutput(mesh,coefs=[felset,gfu.components[0]],names=["lset","u"],subdivision=1)

#--START-- INITIAL CONDITION:
  lfu = LinearForm(fes,"lfu")
  lfu.Add(LFI("xsource",dim=dim,coef=[initialconcentration,ConstantCF(0)]))
  lfu.Assemble()
  innerdofs = BitArray(fes.FreeDofs())
  print(innerdofs)
  for i in range(fes.ndof):
      innerdofs.Clear(i)
  print(innerdofs)
  for elid in mesh.Elements():
      if xfes.XFESpace.GetDomainOfElement(elid.nr)!=0:
          for dofnr in xfes.GetDofNrs(elid):
              innerdofs.Set(dofnr)
  print(innerdofs)
  invm = bfmu.mat.Inverse(innerdofs,inverse="pardiso")

  gfu.vec.data = invm * lfu.vec

  uold.data=gfu.vec
  ExtendGF(gfu,fes,mesh,bfextend=bfextend)
#initial condition output:
  vtkout.Do()
#---END--- INITIAL CONDITION: 

  coeflf = GFCoeff(gfu.components[0])

  lff = LinearForm(h1fes,"lff")
  lff.Add(LFI("source",dim=dim,coef=coeflf))

  phi = gflset.vec
  heapsize = int(10**5)
  # input("before time loop\n")

  bfm.Assemble()
  bfa.Assemble()

  summat = bfa.mat.CreateMatrix()

  Auold = bfau.mat.CreateMatrix()
  Muold = bfmu.mat.CreateMatrix()
  rhsu = lfu.vec # lfu is no longer used... we take the vector for rhsu...

  Auold.AsVector().data = bfau.mat.AsVector().data
  Muold.AsVector().data = bfmu.mat.AsVector().data

  summatu = bfmu.mat.CreateMatrix()

  t=0
  lsetold = phi.CreateVector()
  lsetold.data = lsetold

  innerdofs_next = BitArray(fes.FreeDofs())

  Draw(gflset,mesh,"lset")
  Draw(gfu.components[0],mesh,"u")
  while t < tend:
      print("\rt = ",t,end='')
      t += dt

      uold.data = gfu.vec
      
      bfm.AssembleLinearization(phi)
      bfa.AssembleLinearization(phi)
      lff.Assemble()
      tmp = phi.CreateVector()
      summat.AsVector().data = bfm.mat.AsVector() + dt * bfa.mat.AsVector()
      suminv = summat.Inverse(h1fes.FreeDofs())
      tmp.data = -dt * bfa.mat * phi - (dt*beta)* lff.vec

      # background field:
      # lffbg = LinearForm(h1fes,"lff")
      # lffbg.Add(LFI("source",dim=2,coef=VariableCF("sin(pi*4*(x-"+str(t)+")")))
      # lffbg.Assemble()
      # tmp.data += (dt*beta)* lffbg.vec

      phi.data += suminv * tmp


      fes.Update()
      gfu.Update()
      bfau.Assemble()
      bfmu.Assemble()


      for i in range(fes.ndof):
          innerdofs_next.Clear(i)
      for elid in mesh.Elements():
          if xfes.XFESpace.GetDomainOfElement(elid.nr)!=0:
              for dofnr in xfes.GetDofNrs(elid):
                  innerdofs_next.Set(dofnr)

      innerdofs |= innerdofs_next

      # BEGIN TODO:

      # compute summatu and rhsu from old and new time step matrices and parameters theta and dt.. 
      summatu.AsVector().data = bfmu.mat.AsVector()+ dt * theta * bfau.mat.AsVector() + dt * (1-theta) * Auold.AsVector()
      rhsu.data = Muold * gfu.vec
      
      utmp.data = Muold * uone
      print("\nmassold = ",InnerProduct(uone,rhsu))
      print("\nradius = ",InnerProduct(uone,utmp)/pi)
      
      # END TODO:
      
      ainv = summatu.Inverse(innerdofs)

      for i in range(fes.ndof):
          if(innerdofs_next[i]):
            innerdofs.Set(i)
          else:
            innerdofs.Clear(i)

      gfu.vec.data = ainv * rhsu
      
      ExtendGF(gfu,fes,mesh,bfextend=bfextend,setzero= not smoothextension) #,bfmass=bfmass,vecuold=uold)
      #extension with inertia:
      #ExtendGF(gfu,fes,mesh,bfextend=bfextend,bfmass=bfmass,vecuold=uold)

      #store laplace and mass matrices w.r.t. to "current" (future "old") interface:
      Auold.AsVector().data = bfau.mat.AsVector().data
      Muold.AsVector().data = bfmu.mat.AsVector().data
      
      Redraw(blocking=True)
      # vtkout.Do()
