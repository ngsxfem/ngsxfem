from math import pi
from ngsolve import *
from netgen.csg import OrthoBrick, Pnt

def IsEqual(A:BitArray, B:BitArray) -> bool:
    """Checks if A and B are elementwise equal"""
    return ((A.__and__(~B)).NumSet() == 0) and (((B).__and__(~A)).NumSet() == 0)

def AdjacencyMatrix (mesh, neighbortype : str = "vertex") -> BaseMatrix :
    """
    returns the adjacency matrix of the mesh elements w.r.t. vertices, edges or faces
    the return object is a BaseMatrix (double) with entries 0.0 or 1.0
    """    
    P0 = L2(mesh,order=0)
    if neighbortype == "vertex":
        W = H1(mesh,order=1)
    elif neighbortype == "face":
        W = HDiv(mesh,order=0)
    elif neighbortype == "edge":
        W = HCurl(mesh,order=0)
    else:
        raise Exception("Unknown neighbor type: "+neighbortype+"\n allowed are vertex, face and edge.")

    B = BilinearForm(testspace=P0,trialspace=W)
    B.Assemble()
    B.mat.AsVector()[:] = 1

    return B.mat @ B.mat.T
        
def AddNeighborhood(marker : BitArray, adja : BaseMatrix, layers : int = 1, inplace : bool = False) -> BitArray:
    """
    returns a BitArray which extends a marker of elements in a mesh with its neighborhood.
    The neighborhood consists of a number of layers of neighbors w.r.t. the provided adjacency matrix,
    see also `AdjacencyMatrix`.
    """
    mv = adja.CreateColVector()
    mw = adja.CreateColVector()
    for i, b in enumerate(marker):
        mv[i] = 1 if b else 0

    for i in range(layers):
        if i > 0:
            mv.data = mw
        mw.data = adja * mv

    if not inplace:
        rmarker = BitArray(len(marker))
    else:
        rmarker = marker

    for i in range(len(marker)):
        rmarker[i] = mw[i] > 0.5

    return rmarker



LevelsetExamples = {
    # geometry from 'Dziuk, Elliott, Finite element methods for surface PDEs, Acta Numerica, 2013', pp. 373-374:
    "cheese" : (sqrt((x*x-1)*(x*x-1)+(y*y-1)*(y*y-1)+(z*z-1)*(z*z-1)+(x*x+y*y-4)*(x*x+y*y-4)+(x*x+z*z-4)*(x*x+z*z-4)+(y*y+z*z-4)*(y*y+z*z-4))-4).Compile(),
    # geometry from 'Dziuk, Elliott, Finite element methods for surface PDEs, Acta Numerica, 2013', pp. 318-319:
    "dziukelliott" : sqrt(0.25*x*x+y*y+4.0*z*z/((1+0.5*sin(pi*x))*(1+0.5*sin(pi*x))))-1.0,
    # geometry from 'Dziuk, Finite elements for the beltrami operator on arbitrary surfaces':
    "dziuk88" : sqrt( (x-z*z)*(x-z*z)+y*y+z*z) - 1.0,
    # the unit sphere
    "sphere" : sqrt(x*x+y*y+z*z)-1.0,
    # torus with parameters as in 'Grande, Reusken, A higher order finite element method for partial differential euqations on surface, SINUM, 2016'
    "torus" : sqrt(z*z + (sqrt(x*x+y*y) - 1.0)*(sqrt(x*x+y*y) - 1.0)) - 0.6,
    # gyroid from Lehrenfeld, CMAME, 2016
    "gyroid" : cos(pi*x)*sin(pi*y)+cos(pi*y)*sin(pi*z)+cos(pi*z)*sin(pi*x),
}

BoundingBoxes = {
    "cheese" : OrthoBrick(Pnt(-2.5,-2.5,-2.5), Pnt(2.5,2.5,2.5)),
    "dziukelliott" : OrthoBrick(Pnt(-2.5,-1.5,-1.5), Pnt(2.5,1.5,1.5)),
    "dziuk88" : OrthoBrick(Pnt(-2,-2,-2), Pnt(2,2,2)),
    "sphere" : OrthoBrick(Pnt(-1.5,-1.5,-1.5), Pnt(1.5,1.5,1.5)),
    "torus" : OrthoBrick(Pnt(-2,-2,-2), Pnt(2,2,2)),
    "gyroid" : OrthoBrick(Pnt(-1,-1,-1), Pnt(1,1,1)),
}
