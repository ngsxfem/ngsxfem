from ngsolve import *
from xfem import *
from enum import Enum

class SmoothBlendingType(Enum):
    BLEND_AUTODETECT = 1
    BLEND_FIXED_WIDTH = 2

globals().update(SmoothBlendingType.__members__)

class LevelSet_SmoothBlending:

    def __init__(self, sb_type, option_dict):
        self.smooth_blend_type = sb_type
        if sb_type == BLEND_AUTODETECT or sb_type == BLEND_FIXED_WIDTH:
            if isinstance(option_dict, dict):
                self.lset = option_dict["lset"]
                self.width = option_dict["width"]
                self.order = option_dict["order"]
                if "time_order" in option_dict:
                    self.time_order = option_dict["time_order"]
                else:
                    self.time_order = -1
                if sb_type == BLEND_AUTODETECT:
                    self.mesh = option_dict["mesh"]
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
            self.CF = IfPos(lset_treshold - sqrt(self.lset**2), 0, IfPos( sqrt(self.lset**2) - (lset_treshold + self.width), 1, blend_1d_polyn( ( sqrt(self.lset**2) - lset_treshold)/(self.width) ) ))
