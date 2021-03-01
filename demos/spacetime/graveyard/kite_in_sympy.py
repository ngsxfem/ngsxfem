from sympy import *
from sympy.printing import latex, pretty
from sympy.parsing.sympy_parser import parse_expr

r0 = 1
Q = pi/r0

def get_Laplace_u_str(rhostr, ustr):
    x,y,t = symbols("x y t")
    
    rho = eval(rhostr) #(1 - y**2)*t
    r = sqrt((x - rho)**2+y**2)
    
    u = eval(ustr) #cos(Q*r) * sin(pi*t)
    Laplaceu = diff(u,x,2) + diff(u,y,2)
    
    #print("Laplace: ", Laplaceu)
    return str(Laplaceu)
    #Laplaceu = - Q*sin(pi*t)*(cos(Q*r)*Q* ( drdx**2 + drdy**2) + sin(Q*r)*( dr2dx2 + dr2dy2 ))
    #coeff_f = -Laplaceu + pi * cos(Q*r) * cos(pi*t)

#print ( get_Laplace_u_str("(1 - y**2)*t", "cos(Q*r) * sin(pi*t)"))

def get_dt_rho(rhostr):
    x,y,t = symbols("x y t")
    
    rho = eval(rhostr)
    return str(diff(rho, t))
