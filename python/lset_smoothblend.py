from ngsolve import *
from xfem import *
from enum import Enum

class SmoothBlendingType(Enum):
    BLEND_AUTODETECT = 1
    BLEND_FIXED_WIDTH = 2

globals().update(SmoothBlendingType.__members__)

class LevelSet_SmoothBlending:
    """
    Class to supply the user with an automatically calculated smooth blending function b
    used in LevelSetMeshAdaptation_Spacetime objects. Call the constructor of this class to obtain
    an object with the member function CF, which is the sought-after b prepared for you by the
    constructor. The constructor shall be called with a blending type, either BLEND_FIXED_WIDTH
    or BLEND_AUTODETECT, as the first argument and a dictionary with options as the second argument.
    What shall be supplied in the dictionary differs by the blending types.

    About the two blending types in detail:

    * BLEND_FIXED_WIDTH. Required keys in the dict: { "lset": CF, "order": int , "width": float}
      In this case, the blending is assumed to decay with a fixed width w. The cut topology (time
      dependent or not) is specified by the lset Coefficientfunction. Starting from those points
      with lset=0 on the cut interface, we assume that the lset function is a measure of the
      distance from the interface. (If it is only equivalent to such a function up to a non-
      degenerating constant this is not a problem, only to the extend that the resulting width
      might eventually differ by the same constant from w) If the distance ("measured" as such)
      is between 0 and w, b equals 0. If the distance is between w and 2*w, the blending will
      increase from 0 to 1 by the polynomial pi_s of order order as defined in
      https://arxiv.org/abs/2311.02348. For distances greater than 2*w, b=1,
      asking for a zero deformation.

    * BLEND_AUTODETECT. Required keys in the dict: { "lset": CF, "order": int, "mesh": mesh,
                                                    "time_order": int}
      In this function, the overall structure of the BLEND_FIXED_WIDTH remains the same.
      However, the parameter w does not have to be provided explicitely. Instead, it is evaluted
      by a numerical estimate of all the cut elements inside the current time slice as prescribed
      by lset. To this extend, the further arguments of a mesh and a time_order are needed.
      The latter is used to generate sampling points in time for the maximum distance estimation.
      This means that e.g. time_order 1 or 2 would be typical reasonable choices.

    """
    def __init__(self, sb_type, option_dict):
        self.smooth_blend_type = sb_type
        if sb_type == BLEND_AUTODETECT or sb_type == BLEND_FIXED_WIDTH:
            if isinstance(option_dict, dict):
                self.lset = option_dict["lset"]
                self.order = option_dict["order"]
                if sb_type == BLEND_FIXED_WIDTH:
                    self.width = option_dict["width"]
                if sb_type == BLEND_AUTODETECT:
                    self.mesh = option_dict["mesh"]
                    self.time_order = option_dict["time_order"]
                self.Calc_blend_CF()
            else:
                print("Error: LevelSetMeshAdaptation_SmoothBlending needs a dictionary when called with BLEND_AUTODETECT or BLEND_FIXED_WIDTH")
        else:
            print("Error: Illegal first argument in LevelSetMeshAdaptation_SmoothBlending")

    def Calc_blend_CF(self):
        blend_1d_polyn = lambda x: IfPos( 0.5 - x, 2**(self.order-1)*x**self.order, 1-2**(self.order-1) *(1-x)**self.order)
        if self.smooth_blend_type == BLEND_FIXED_WIDTH:
            self.CF = IfPos(self.width - sqrt(self.lset**2), 0, IfPos(sqrt(self.lset**2) - 2*self.width, 1, blend_1d_polyn( ( sqrt(self.lset**2) - self.width)/(self.width) ) ))
        elif self.smooth_blend_type == BLEND_AUTODETECT:
            ci = CutInfo(self.mesh, time_order=self.time_order)
            fes_space = H1(self.mesh, order=1)
            fes_time = ScalarTimeFE(self.time_order)
            fes_st = SpaceTimeFESpace(fes_space,fes_time)
            lset_node_p1 = GridFunction(fes_space)
            self.lset_p1 = GridFunction(fes_st)
            times = [xi for xi in fes_st.TimeFE_nodes()]
            for i,ti in enumerate(times):
                lset_node_p1.Set(fix_tref(self.lset,ti))
                self.lset_p1.vec[i*len(lset_node_p1.vec) : (i+1)*len(lset_node_p1.vec)].data = lset_node_p1.vec[:]

            ci.Update(self.lset_p1 , time_order=0)

            minv, maxv = IntegrationPointExtrema({"levelset": 2*IndicatorCF(self.mesh, ci.GetElementsOfType(IF)) -1, "domain_type" : POS, "order": 3}, self.mesh, fix_tref(self.lset,0.5))
            lset_treshold = max(abs(minv), abs(maxv))
            self.CF = IfPos(lset_treshold - sqrt(self.lset**2), 0, IfPos( sqrt(self.lset**2) - (2*lset_treshold), 1, blend_1d_polyn( ( sqrt(self.lset**2) - lset_treshold)/(lset_treshold) ) ))
