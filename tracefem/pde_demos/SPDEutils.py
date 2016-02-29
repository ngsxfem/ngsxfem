# ngsolve stuff
from ngsolve import *
from xfem.tracefem import *
from xfem.lsetcurv import *


def MakeDefaultOptions():
    return {"VTKOutput" : False, "VTKSubdivision" : 2, "Plotting" : False, "CSVOutput" : False }

def MakeDefaultDiscretizationParams():
    return {"Order" : 2, "DiscontinuousQuasiNormal" : True, "FullGradients" : False}

def MakePDEObjects(mesh, problemdata, discparams):
    pdeobj = {}
    pdeobj["Vh"] = H1(mesh, order=discparams["Order"], dirichlet=[])
    pdeobj["Vh_tr"] = TraceFESpace(mesh, pdeobj["Vh"], problemdata["Levelset"])
    
    pdeobj["a"] = BilinearForm(pdeobj["Vh_tr"], symmetric = (problemdata["Convection"] == None) , flags = { })
    if (problemdata["Reaction"] != None):
        pdeobj["a"] += TraceMass(problemdata["Reaction"])
    if (problemdata["Diffusion"] != None):
        if (discparams["FullGradients"]):
            pdeobj["a"] += TraceLaplace(problemdata["Diffusion"])
        else:
            pdeobj["a"] += TraceLaplaceBeltrami(problemdata["Diffusion"])
    if (problemdata["Convection"] != None):
        pdeobj["a"] += TraceConvection(problemdata["Convection"])
    # pdeobj["a"] += BFI("normallaplacetrace", coef=[0.1, (x,y,z)])

    pdeobj["f"] = LinearForm(pdeobj["Vh_tr"])
    pdeobj["f"] += TraceSource(problemdata["Source"])

    pdeobj["c"] = Preconditioner(pdeobj["a"], type="local", flags= { "test" : True })
    #pdeobj["c"] = Preconditioner(pdeobj["a"], type="direct", flags= { "inverse" : "sparsecholesky", "test" : True })

    pdeobj["u"] = GridFunction(pdeobj["Vh_tr"])
    return pdeobj;

def MakeResultContainer(pdeobj, lsetmeshadap, problemdata, options):
    results = {}
    if (options["VTKOutput"]):
        results["vtk"] = VTKOutput(ma=problemdata["Mesh"],
                                   coefs=[problemdata["Levelset"],
                                          lsetmeshadap.lset_ho,
                                          lsetmeshadap.lset_p1,
                                          lsetmeshadap.deform,
                                          pdeobj["u"],
                                          problemdata["Solution"]],
                                   names=["lset",
                                          "lsetho",
                                          "lsetp1",
                                          "deform",
                                          "u",
                                          "uexact"],
                                   filename="LapBeltrResults/vtkout_",
                                   subdivision=options["VTKSubdivision"])
    else:
        results["vtk"] = None
    results["statistics"] = { "level" : [],
                              "ndof" : [],
                              "maxdist" : [],
                              "l2err" : [],
                              "maxerr" : [],
                              "numits" : []
    }

    results["last_num_its"] = 0
    results["level"] = 0
    
    return results;
    
    
def SolveProblem(data):
    problemdata = data["problemdata"]
    pdeobj = data["pdeobj"]
    results = data["results"]
    lsetmeshadap = data["lsetmeshadap"]

    # Calculation of the deformation:
    deformation = lsetmeshadap.CalcDeformation(problemdata["Levelset"])
    # Applying the mesh deformation
    problemdata["Mesh"].SetDeformation(deformation)

    pdeobj["Vh"].Update()
    # print("StdFESpace NDof:", pdeobj["Vh"].ndof)
    pdeobj["Vh_tr"].Update()
    pdeobj["u"].Update()
    
    pdeobj["a"].Assemble();
    pdeobj["f"].Assemble();
    pdeobj["c"].Update();
    
    solvea = CGSolver( mat=pdeobj["a"].mat, pre=pdeobj["c"].mat, complex=False, printrates=False, precision=1e-8, maxsteps=20000)
    # # the boundary value problem to be solved on each level
    pdeobj["u"].vec.data = solvea * pdeobj["f"].vec;
    
    results["last_num_its"] = solvea.GetSteps()
    problemdata["Mesh"].UnsetDeformation()

def PostProcess(data):
    problemdata = data["problemdata"]
    pdeobj = data["pdeobj"]
    results = data["results"]
    lsetmeshadap = data["lsetmeshadap"]
    discparams = data["discparams"]
    
    maxdistlset = lsetmeshadap.CalcMaxDistance(problemdata["Levelset"]);
    # maxdistlsetho = lsetmeshadap.CalcMaxDistance(lsetmeshadap.lset_ho);
    problemdata["Mesh"].SetDeformation(lsetmeshadap.deform)
    [l2diff,maxdiff] = CalcTraceDiff( pdeobj["u"], problemdata["Solution"], intorder=2*discparams["Order"]+2)
    problemdata["Mesh"].UnsetDeformation()

    print ("The mesh Refinement level :", results["level"])
    print ("Number of d.o.f. Std. -FES:", pdeobj["Vh"].ndof)
    print ("Number of d.o.f. Trace-FES:", pdeobj["Vh_tr"].ndof)
    print ("Max. distance to interface:", maxdistlset)
  # print ("Max. distance to ho-intfce:", maxdistlsetho)
    print ("L2-norm error on interface:", l2diff)
    print ("maxnorm error on interface:", maxdiff)
    print ("Numb. of CG-Solver iterat.:", results["last_num_its"])
    
    results["statistics"]["level"].append(results["level"])
    results["statistics"]["ndof"].append(pdeobj["Vh_tr"].ndof)
    results["statistics"]["maxdist"].append(maxdistlset)
    results["statistics"]["l2err"].append(l2diff)
    results["statistics"]["maxerr"].append(maxdiff)
    results["statistics"]["numits"].append(results["last_num_its"])
    
    RefineAtLevelSet(gf=lsetmeshadap.lset_p1)
    
    if (results["vtk"] != None):
        results["vtk"].Do()
        
    Draw(lsetmeshadap.lset_p1,problemdata["Mesh"],"lsetp1")
    Draw(lsetmeshadap.deform,problemdata["Mesh"],"deformation")
    Draw(pdeobj["u"],problemdata["Mesh"],"u",draw_surf=False)


def PlotConvergence(data):
    results = data["results"]
    fo = open("LapBeltrResults/l2err.csv","w");
    fo.writelines(["{}: {}\n".format(a,b) for a,b in zip(results["statistics"]["level"],results["statistics"]["l2err"])])
    fo = open("LapBeltrResults/ndof.csv","w");
    fo.writelines(["{}: {}\n".format(a,b) for a,b in zip(results["statistics"]["level"],results["statistics"]["ndof"])])
    fo = open("LapBeltrResults/geomerr.csv","w");
    fo.writelines(["{}: {}\n".format(a,b) for a,b in zip(results["statistics"]["level"],results["statistics"]["maxdist"])])
    fo = open("LapBeltrResults/maxerr.csv","w");
    fo.writelines(["{}: {}\n".format(a,b) for a,b in zip(results["statistics"]["level"],results["statistics"]["maxerr"])])
    fo = open("LapBeltrResults/itcnts.csv","w");
    fo.writelines(["{}: {}\n".format(a,b) for a,b in zip(results["statistics"]["level"],results["statistics"]["numits"])])
    
    plt.figure(0)
    plt.yscale('log')
    plt.xlabel("level")
    
    plt.plot(results["statistics"]["level"],results["statistics"]["l2err"], "-*")
    plt.plot(results["statistics"]["maxdist"], "-+")
    plt.legend(["L2error","geometry (max) error"])

    # plt.figure(1)
    # plt.yscale('log')
    # plt.xlabel("level")

    # plt.plot(results["statistics"]["level"],results["statistics"]["ndof"], "-*")
    # plt.plot(results["statistics"]["numits"], "-+")
    # plt.legend(["D.o.f.","iterations"])

    # plt.figure(2)
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel("ndof")
    # plt.ylabel("error")

    # plt.plot(results["statistics"]["ndof"],results["statistics"]["l2err"], "-*")
    # plt.legend(["L2error"])
    
    plt.ion()
    plt.show()

    
    l2conv = [ log(results["statistics"]["l2err"][i]/results["statistics"]["l2err"][i-1])/log(0.5) for i in range(1,len(results["statistics"]["l2err"]))]
    maxconv = [ log(results["statistics"]["maxerr"][i]/results["statistics"]["maxerr"][i-1])/log(0.5) for i in range(1,len(results["statistics"]["maxerr"]))]
    geomconv = [ log(results["statistics"]["maxdist"][i]/results["statistics"]["maxdist"][i-1])/log(0.5) for i in range(1,len(results["statistics"]["maxdist"]))]
    print ("l2err convergence orders (eoc):", l2conv)
    print ("geom. convergence orders (eoc):", geomconv)
    if (not hasattr(sys,'ps1')):
        input("<press enter to quit>")
