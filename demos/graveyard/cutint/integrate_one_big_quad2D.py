"""
integration on levelset domains
"""

# ------------------------------ LOAD LIBRARIES -------------------------------
from ngsolve import *
from xfem import *
from sympy import *


# ----------------------------------- MAIN ------------------------------------
def get_levelset(lsetvals):
    lset = lsetvals[0] + (lsetvals[1] - lsetvals[0]) * x
    lset += (lsetvals[3] - lsetvals[0]) * y
    lset += (lsetvals[2] - lsetvals[1] - lsetvals[3] + lsetvals[0]) * x * y
    return lset


def get_referencevals(lsetvals, f, IF_INT=False):
    d = lsetvals[0]
    c = lsetvals[2] - lsetvals[1] - lsetvals[3] + lsetvals[0]
    if abs(c) < 1e-12:
        c = 0
    a = lsetvals[1] - lsetvals[0]
    b = lsetvals[3] - lsetvals[0]
    # print("a,b,c,d: ",a,b,c,d)
    xs = Symbol('xs')
    ys = Symbol('ys')
    levelset_py = d + a * xs + b * ys + c * xs * ys
    f_py = f(xs, ys)

    def levelset_pyl(x, y):
        return levelset_py.subs(ys, y).subs(xs, x)

    root1 = solve(levelset_py.subs(ys, 0), xs)
    root2 = solve(levelset_py.subs(ys, 1), xs)
    # print(root1, root2)
    part = [0, 1]
    for r in [root1[0], root2[0]]:
        if r > 0 and r < 1:
            part.append(r)
    part.sort()
    # print("Partitioning of [0,1]: ", part)

    referencevals = {}
    referencevals[POS] = 0
    referencevals[NEG] = 0
    if IF_INT:
        referencevals[IF] = 0

    def add_integration_on_interval(x0, x1):
        # print("Integrating on interval: ",x0,x1)
        y_ast = -(a * xs + d) / (b + c * xs)
        cut_left = levelset_pyl(x0, 0) * levelset_pyl(x0, 1) < 0
        cut_right = levelset_pyl(x1, 0) * levelset_pyl(x1, 1) < 0
        # print("Lset vals y = 0: ", levelset_pyl(x0,0), levelset_pyl(x1,0))
        # print("Lset vals y = 1: ", levelset_pyl(x0,1), levelset_pyl(x1,1))
        # print("Cuts: ", cut_left, cut_right)
        if (not cut_left) and (not cut_right):
            I0 = integrate(integrate(f_py, (ys, 0, 1)), (xs, x0, x1))
            if levelset_pyl(x0, 0) + levelset_pyl(x1, 0) > 0:
                print("Adding ", I, " to Part POS")
                referencevals[POS] += I0
            else:
                print("Adding ", I, " to Part NEG")
                referencevals[NEG] += I0
        else:
            I1 = integrate(integrate(f_py, (ys, 0, y_ast)), (xs, x0, x1))
            I2 = integrate(integrate(f_py, (ys, y_ast, 1)), (xs, x0, x1))
            Ia = integrate(f_py, (ys, 0, y_ast))
            # print("I1:", I1)
            # print("I2:", I2)

            # print("Interval length: ", x1-x0)
            # print("Sum of I1 and I2:", I1+I2)

            if(levelset_pyl(x0, 0) + levelset_pyl(x1, 0) > 0):
                referencevals[POS] += I1
                referencevals[NEG] += I2
            else:
                referencevals[POS] += I2
                referencevals[NEG] += I1

            if IF_INT:
                ts = Symbol('ts')
                C = Curve([ts, -(a * ts + d) / (b + c * ts)], (ts, x0, x1))
                print(C)
                referencevals[IF] += line_integrate(f_py, C, [xs, ys])

    for i in range(0, len(part) - 1):
        add_integration_on_interval(part[i], part[i + 1])
    return referencevals


if __name__ == "__main__":
    from netgen.geom2d import SplineGeometry
    square = SplineGeometry()
    square.AddRectangle([0, 0], [1, 1], bc=1)
    # mesh = Mesh (square.GenerateMesh(maxh=100, quad_dominated=False))
    mesh = Mesh(square.GenerateMesh(maxh=100, quad_dominated=True))

    lsetvals_list = [[-0.18687, 0.324987, 0.765764, 0.48983],
                     [0.765764, 0.324987, -0.18687, -0.48983],
                     [1, 2 / 3, -1, -2 / 3]]
    # lsetvals_list = [[1.,-1.,-3.,-1.]]
    # lsetvals_list = [[1,-1,-4,-2]]
    lsetvals_list.append([3, -1, 1, -1.023123])

    def f(x, y):
        return 1 + 0 * x + 0 * y

    f_ngs = f(x, y)
    V = H1(mesh, order=1)
    lset_approx = GridFunction(V)

    # domains = [NEG,POS, IF]
    domains = [NEG, POS]
    # domains = [IF]
    error_list = []

    max_order = 10
    # f1 = open("errors.dat","w")

    for lsetvals in lsetvals_list:
        print("Case lsetvals = ", lsetvals)
        referencevals = get_referencevals(lsetvals, f, IF in domains)
        print("referencevals: ", referencevals)
        levelset = get_levelset(lsetvals)

        InterpolateToP1(levelset, lset_approx)

        errors = dict()

        for key in domains:
            errors[key] = []
        inte = dict()

        for order in range(max_order + 1):
            errs = dict()
            for key in domains:
                lset_dom = {"levelset": lset_approx, "domain_type": key}
                integral = Integrate(levelset_domain=lset_dom,
                                     cf=f_ngs, mesh=mesh, order=order)
                inte[key] = integral
                print("PP: Integral on Domain ", key, " : ", integral)
                errors[key].append(abs(integral - referencevals[key]))
                errs[key] = abs(integral - referencevals[key])
            # f1.write(str(order) + "\t"
            #          + str(sqrt(pow(errs[POS],2.) + pow(errs[NEG],2.)))
            #          + "\n")
            # print("Sum of Part NEG, POS: ", inte[NEG] + inte[POS])
        print("L2 Errors:", errors)
        error_list.append(errors)
        # f1.write("\n\n")
        # Draw(levelset, mesh, "lset")
    print("All L2 Errors", error_list)
