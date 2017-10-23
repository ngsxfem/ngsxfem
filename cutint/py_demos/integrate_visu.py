from math import pi
# ngsolve stuff
from ngsolve import *
# basic xfem functionality
from xfem import *
import os

from netgen.geom2d import SplineGeometry

square = SplineGeometry()
square.AddRectangle([0,0],[1,1],bc=1)
mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=True))

lsetvals_list = [ [1.,-1.,-1.1,-1.] , [1.,0.7,-1,-0.6], [1., -0.3, -1.3, 0.55] ]

V = H1(mesh,order=1)
lset_approx = GridFunction(V)

tikz_styles ={POS: "dt_pos", NEG: "dt_neg", IF: "dt_if"}

lset_i = 0
for lsetvals in lsetvals_list:
    print(lsetvals)
    levelset = lsetvals[0] +(lsetvals[1] - lsetvals[0])*x + (lsetvals[3] - lsetvals[0])*y + (lsetvals[2]-lsetvals[1]-lsetvals[3]+lsetvals[0])*x*y
    InterpolateToP1(levelset,lset_approx)
    
    for order in [4]:
        tikz_out = open("tikz_out.tex", "w")
        tikz_out.write("\\plotsquareiso{"+str(lsetvals[0])+"}{"+str(lsetvals[1])+"}{"+str(lsetvals[3])+"}{"+str(lsetvals[2])+"}\n")
        for domain in [POS, NEG, IF]:
            os.system("rm cutrule.dat")
            #Generate cutrule.dat
            I = Integrate(levelset_domain = { "levelset" : lset_approx, "domain_type" : domain}, cf=CoefficientFunction(1), mesh=mesh, order = order)
            #Parse cutrule.dat
            f = open("cutrule.dat")
            
            for l in f:
                if len(l.split("\t")) > 1:
                    xval = l.split("\t")[0]
                    yval = l.split("\t")[1]
                    
                    v = l.split("\t")[3]
                    tikz_out.write("\\draw ("+str(xval)+", "+str(yval)+") node ["+tikz_styles[domain]+"] {};\n")
            f.close()
        tikz_out.close()
        os.system("lualatex -shell-escape tikz_frame.tex > tex.log")
        os.system("cp tikz_frame.pdf tikz_frame_lsetvals_"+str(lset_i)+".pdf")
        lset_i += 1
