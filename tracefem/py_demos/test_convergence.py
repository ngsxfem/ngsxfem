import argparse        
import pickle
# for making a directory (if it doesn't exist)
import os

from surfpde_probs import *
from discretization_hdgtracefem import *
from discretization_tracefem import *
from eval_sparsity import *
from eval_errors import *

ngsglobals.numthreads=16

def PrintTimers(substring="HDG"):            
    ### print hdg-intergrator timers
    hdgtimers = [a for a in Timers() if substring in a["name"]]
    hdgtimers = sorted(hdgtimers, key=lambda k: k["time"], reverse=True)
    for timer in hdgtimers:
        print("{:<45}: {:6} cts, {:.8f} s, {:.6e} s(avg.)"
              .format(timer["name"],timer["counts"],timer["time"],timer["time"]/timer["counts"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Solve Laplace-Betrami example problem.')
    parser.add_argument('--reflvls', metavar='N', type=int, default=2, help='number of refinement levels (>=1)')
    parser.add_argument('--minorder', metavar='N', type=int, default=1, help='min k')
    parser.add_argument('--maxorder', metavar='N', type=int, default=4, help='max k')
    parser.add_argument('--global_ndof_limit', metavar='N', type=int, default=-1, help='max. number of degrees of freedom (global)')
    parser.add_argument('--addto', dest="former_results", nargs='?', default="results.pkl", help='reuse old results from file')
    parser.add_argument('--vtkout', dest='vtkout', action='store_true', help='enable VTKOutput')
    parser.add_argument('--no-vtkout', dest='vtkout', action='store_false', help='disable VTKOutput')
    parser.set_defaults(vtkout=False)
    args = parser.parse_args()
    options = vars(args)
    print("call arguments: ", options)

    if (options["global_ndof_limit"] == -1):
        global_ndofs_limit = None
    else: 
        global_ndofs_limit = options["global_ndof_limit"]
        options["reflvls"] = 20

    with TaskManager():
        problemdata = Make3DProblem()

        orders = range(options["minorder"],options["maxorder"]+1)

        resultdict = {}
        if os.path.isfile(options["former_results"]):
            with open(options["former_results"], 'rb') as f:
                resultdict = pickle.load(f)

        discretization = None
        trynextlevel = True
        for i in range(options['reflvls']):
            if (i!=0):
                problemdata["Mesh"].Refine()
            trynextorder = True
            for order in orders:
                print("""
                ------------------------------------------------------------------
                   refinement level = {}, order = {}
                ------------------------------------------------------------------
                """.format(i,order))
                
                if (i,order) in resultdict:
                    print("found results in dictionary")
                    continue

                problemdata["Order"] = order

                if (problemdata["HDG"]):
                    discretization = HDGTraceFEMDiscretization(problemdata)
                else:
                    discretization = TraceFEMDiscretization(problemdata)

                global_ndofs = discretization.Vh_tr.ndof
                for j in range(discretization.Vh_tr.ndof):
                    if (discretization.Vh_tr.CouplingType(j) == COUPLING_TYPE.LOCAL_DOF):
                        # discretization.Vh_tr.FreeDofs()[i] = 0
                        global_ndofs = global_ndofs - 1
                if (global_ndofs_limit==None) or (global_ndofs < global_ndofs_limit):
                    resultdict[(i,order)] = discretization.SolveProblem(firstcall=True)
                    EvalErrorsOnInterface(problemdata["Mesh"],
                                          discretization.lsetmeshadap,
                                          discretization.VolumeSolution(),
                                          discretization.order,
                                          problemdata,
                                          resultdict[(i,order)])

                    if ("checkDGpattern" in problemdata and problemdata["checkDGpattern"]):
                        CheckDGpattern(problemdata["Mesh"],order,problemdata["Levelset"],resultdict[(i,order)])
                    if ("checkCGpattern" in problemdata and problemdata["checkCGpattern"]):
                        CheckCGpattern(problemdata["Mesh"],order,problemdata["Levelset"],resultdict[(i,order)])
                    if ("checkCGGPpattern" in problemdata and problemdata["checkCGGPpattern"]):
                        CheckCGGPpattern(problemdata["Mesh"],order,problemdata["Levelset"],resultdict[(i,order)])
                    if ("checkHDGpattern" in problemdata and problemdata["checkHDGpattern"]):
                        CheckHDGpattern(problemdata["Mesh"],order,problemdata["Levelset"],resultdict[(i,order)])
                        
                    print("simulation results:\n {}".format(resultdict[(i,order)]))
                else:
                    print("global ndofs({}) > {}.".format(global_ndofs,global_ndofs_limit))
                    trynextorder = False
                    if (order==orders[0]):
                        trynextlevel = False
                
                if (options["former_results"] != "None"):
                    with open(options["former_results"], 'wb') as f:
                        pickle.dump(resultdict, f, 0) #pickle.HIGHEST_PROTOCOL)
                else:
                    with open('results.pkl', 'wb') as f:
                        pickle.dump(resultdict, f, 0) #pickle.HIGHEST_PROTOCOL)
                if (not trynextorder):
                    break
            if trynextlevel:
                if (discretization == None):

                    lsetmeshadap = LevelSetMeshAdaptation(problemdata["Mesh"], order=1, threshold=1000, discontinuous_qn=True,heapsize=10000000)
                    lsetmeshadap.CalcDeformation(problemdata["Levelset"])
                    RefineAtLevelSet(gf=lsetmeshadap.lset_p1)
                else:
                    RefineAtLevelSet(gf=discretization.lsetmeshadap.lset_p1)
            else:
                break


    for i in range(options['reflvls']):
        for order in orders:
            if (i,order) in resultdict:
                print("({},{}): {}".format(i,order,resultdict[(i,order)]))
            
    # print(resultdict)
    #Draw(discretization.u.components[0],discretization.mesh,"u",draw_surf=False)

    PrintTimers(substring="HDG")
    PrintTimers(substring="Trace")
    PrintTimers(substring="XFE")
    PrintTimers(substring="IntegrateX")
    PrintTimers(substring="LsetCurv")
    
    if (options['vtkout']):
        vtk = VTKOutput(ma=discretization.mesh,coefs=[discretization.lsetmeshadap.lset_p1,
                                                      discretization.lsetmeshadap.deform,
                                                      discretization.VolumeSolution()],
                        names=["lsetp1","deform","u"],filename="vtkout_",subdivision=0)
        vtk.Do()
