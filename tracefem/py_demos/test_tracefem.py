from surfpde_probs import *
from discretization_tracefem import *
from eval_errors import *

problemdata = Make3DProblem_Diffusion()
problemdata["Mesh"].Refine()
discretization = TraceFEMDiscretization(problemdata)
discretization.SolveProblem(firstcall=True)
results={}
EvalErrorsOnInterface(problemdata["Mesh"],
                      discretization.lsetmeshadap,
                      discretization.VolumeSolution(),
                      discretization.order,
                      problemdata,results)
