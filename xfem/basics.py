from ngsolve.comp import *
from ngsolve.fem import *
from libngsxfem_py.xfem import *

def TwoDomainMassIntegrator (coefneg,coefpos):
    return BFI("xmass", coef=[coefneg,coefpos])

def TwoDomainSourceIntegrator (coefneg,coefpos):
    return LFI("xsource", coef=[coefneg,coefpos])



negative_domain = dict()
negative_domain["negdomain"] = True

positive_domain = dict()
positive_domain["posdomain"] = True

volume_domains = dict()
volume_domains["negdomain"] = True
volume_domains["posdomain"] = True

interface_domain = dict()
interface_domain["interface"] = True

interface_and_volume_domains = dict()
interface_and_volume_domains["negdomain"] = True
interface_and_volume_domains["posdomain"] = True
interface_and_volume_domains["interface"] = True

def IntegrateOnInterface(lset,mesh,coef,order=5,subdivlvl=0):
    return IntegrateX(lset,mesh,cf_interface=coef,order=order,
                      subdivlvl=subdivlvl,domains=interface_domain)["interface"]

def IntegrateOnPosDomain(lset,mesh,coef,order=5,subdivlvl=0):
    return IntegrateX(lset,mesh,cf_pos=coef,order=order,
                      subdivlvl=subdivlvl,domains=positive_domain)["posdomain"]

def IntegrateOnNegDomain(lset,mesh,coef,order=5,subdivlvl=0):
    return IntegrateX(lset,mesh,cf_neg=coef,order=order,
                      subdivlvl=subdivlvl,domains=negative_domain)["negdomain"]

def IntegrateOnWholeDomain(lset,mesh,cf_neg,cf_pos,order=5,subdivlvl=0):
    ints = IntegrateX(lset,mesh,cf_neg=cf_neg,cf_pos=cf_pos, 
                      order=order,subdivlvl=subdivlvl,domains=volume_domains)
    return ints["negdomain"] + ints["posdomain"]


