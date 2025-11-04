import pytest
from ngsolve import *
from ngsolve.meshes import *
from netgen.geom2d import unit_square
from xfem import *
#from netgen import gui

ngsglobals.msg_level = 1

def test_elementlayers():
    mesh = MakeStructured2DMesh(quads = False, nx=7, ny=7)
    
    lsetp1 = GridFunction(H1(mesh))
    lsetp1.Set(x-1/2)

    from xfem.utils import AddNeighborhood, AdjacencyMatrix, IsEqual

    refba = {"vertex": BitArray(mesh.ne), "edge": BitArray(mesh.ne), "face": BitArray(mesh.ne)}
    for i,b in enumerate("00111111111100001111111111000011111111110000111111111100001111111111000011111111110000111111111100"):
        refba["vertex"][i] = True if b == "1" else False
    for i,b in enumerate("00001111110000000011111100000000111111000000001111110000000011111100000000111111000000001111110000"):
        refba["edge"][i] = True if b == "1" else False
        refba["face"][i] = True if b == "1" else False

    for nbtype in ["vertex", "face", "edge"]:
        ci = CutInfo(mesh, lsetp1)
        marker = ci.GetElementsOfType(IF)
        #Draw(BitArrayCF(marker),mesh,"marker")
        m2 = AddNeighborhood(marker, AdjacencyMatrix(mesh, nbtype), 1, inplace=False)
        m3 = AddNeighborhood(m2, AdjacencyMatrix(mesh, nbtype), 1, inplace=False)

        m4 = AddNeighborhood(marker, AdjacencyMatrix(mesh, nbtype), 2, inplace=False)

        AddNeighborhood(marker, AdjacencyMatrix(mesh, nbtype), 2, inplace=True)

        assert IsEqual(m3, m4)
        assert IsEqual(m3, marker)
        assert IsEqual(m4, marker)
        assert IsEqual(refba[nbtype], marker)
        #print(marker)
    #Redraw()
    #Draw(BitArrayCF(marker),mesh,"marker")


if __name__ == "__main__":
    test_elementlayers()