from ngsolve import *
from time import sleep
from netgen.geom2d import unit_square
from netgen.geom2d import SplineGeometry
from netgen.meshing import MeshingParameters
from ngsolve.internal import *
from xfem import *
from numpy import pi
from xfem.lset_spacetime import *


square = SplineGeometry()
square.AddRectangle([0,0],[2,2],bc=1)
maxh = 0.3
mesh = Mesh (square.GenerateMesh(maxh=maxh, quad_dominated=False))


h = specialcf.mesh_size

# FE Spaces
k_t = 1
k_s = 1
fes1 = H1(mesh,order=k_s,dirichlet=[1,2,3,4])

lset_order_time = 2
tfe = ScalarTimeFE(k_t) 

st_fes = SpaceTimeFESpace(fes1,tfe)
#visoptions.deformation = 1

tend = 1
delta_t = 1/32
tnew = 0

told = Parameter(0)
tref = ReferenceTimeVariable()
t = told + delta_t*tref

lset_adap_st = LevelSetMeshAdaptation_Spacetime(mesh, order_space = k_s, order_time = lset_order_time,
                                         threshold=0.1, discontinuous_qn=True)

# Problem data

problems =  {"planar": ( CoefficientFunction(0.25*t), CoefficientFunction(0.25),
                         CoefficientFunction(1),CoefficientFunction(0),CoefficientFunction(0)),
             "nonlinear": ( CoefficientFunction((0.25/pi)*sin(2*pi*t)),CoefficientFunction(0.5*cos(2*pi*t)),
                            CoefficientFunction(7/8 + (1/4)*y*y*(2-y)*(2-y)),CoefficientFunction( y*(1-y)*(2-y) ),CoefficientFunction( (1-y)*(2-y) - y*(2-y) - y*(1-y))) }
          
rho,d_rho,q_coeff,d_q_coeff,dd_q_coeff = problems["planar"]

r11 = x - q_coeff - rho

D_const = 2/3
levelset= CoefficientFunction(sqrt(r11*r11) - 0.5*D_const)

fes_lset_slice = H1(mesh, order=1, dirichlet=[])
lset_p1_slice = GridFunction(fes_lset_slice)

t_old = 0
ci = CutInfo(mesh)

while tend - t_old > delta_t/2:
    
    dfm = lset_adap_st.CalcDeformation(levelset,told,t_old,delta_t,calc_kappa = True) 
    lset_p1 = lset_adap_st.lset_p1  
    ci.Update(lset_p1,2)
    print("time = {0}".format(t_old))
    hasneg= BitArray(ci.GetElementsOfType(NEG,VOL))
    hasneg_BND = BitArray(ci.GetElementsOfType(NEG,BND))
    print("hasneg_BND:")
    print(hasneg_BND)
    haspos= BitArray(ci.GetElementsOfType(POS,VOL))
    haspos_BND = BitArray(ci.GetElementsOfType(POS,BND))
    print("haspos_BND:")
    print(haspos_BND)
    hasif= BitArray(ci.GetElementsOfType(IF,VOL))
    hasif_BND = BitArray(ci.GetElementsOfType(IF,BND))
    print("hasif_BND:")
    print(hasif_BND)
    Draw(BitArrayCF(hasneg),mesh,"NEG")
    Draw(BitArrayCF(haspos),mesh,"POS")
    Draw(BitArrayCF(hasif),mesh,"IF")
    lset_p1_slice.vec[:] = lset_p1.vec[0:fes_lset_slice.ndof]
    Draw(IfPos(-lset_p1_slice,CoefficientFunction(1),CoefficientFunction(0)),mesh,"lsetp1")
    input("Continue")
    
    t_old = t_old + delta_t
    told.Set(t_old)