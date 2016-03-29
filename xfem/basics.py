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

def IntegrateOnInterface(coef,lset,mesh,order=5,subdivlvl=0):
    return IntegrateX(coef,lset,mesh,order,subdivlvl,interface_domain)["interface"]

def IntegrateOnPosDomain(coef,lset,mesh,order=5,subdivlvl=0):
    return IntegrateX(coef,lset,mesh,order,subdivlvl,positive_domain)["posdomain"]

def IntegrateOnNegDomain(coef,lset,mesh,order=5,subdivlvl=0):
    return IntegrateX(coef,lset,mesh,order,subdivlvl,negative_domain)["negdomain"]


