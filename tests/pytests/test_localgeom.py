from ngsolve import *
from xfem import *
from netgen.geom2d import SplineGeometry

# --- Cut situation to be tested: ---
#
#   (0,1) 
#     \
#     |\
#     | \
#     |  \
#     |   \
#     +    \
#     |\    \
#     | \    \
#     |  \    \
#     |NEG\ POS\
#     +----+----+
#   (0,0)     (1,0)
#
# meas_2(POS) = 3/8, meas_2(NEG) = 1/8, meas_1(IF) = 1/2

# simple test: triangle, straight cut, integrand = 1
def test_cut_triangle():
    geom = SplineGeometry()
    pnts = [ (0,0), (1,0), (0,1) ]
    pnums = [geom.AppendPoint(*p) for p in pnts]
    lines = [(0,1,1,1,0), (1,2,1,1,0), (2,0,1,1,0)]
    for p1,p2,bc,left,right in lines:
        geom.Append( ["line", pnums[p1], pnums[p2]], bc=bc, leftdomain=left, rightdomain=right)
    mesh = Mesh (geom.GenerateMesh(maxh=1)) # <- will have 4 elements

    levelset = x + y - 0.25
    lsetp1 = GridFunction(H1(mesh,order=1))
    InterpolateToP1(levelset,lsetp1)

    lset_neg = { "levelset" : lsetp1, "domain_type" : NEG, "subdivlvl" : 0}
    lset_pos = { "levelset" : lsetp1, "domain_type" : POS, "subdivlvl" : 0}

    for order in range(16):
        measure_neg = Integrate(levelset_domain = lset_neg,cf=CoefficientFunction(1.0), mesh = mesh, order=order)
        measure_pos = Integrate(levelset_domain = lset_pos,cf=CoefficientFunction(1.0), mesh = mesh, order=order)
        assert abs(measure_neg-1.0/32.0) < 5e-16*(order+1)*(order+1)
        assert abs(measure_pos-1.0/2.0+1.0/32.0) < 5e-16*(order+1)*(order+1)
        assert abs(measure_neg+measure_pos-1.0/2.0) < 5e-16*(order+1)*(order+1)

