from math import log
import pickle


import argparse        

parser = argparse.ArgumentParser(description='Postprocess laplace beltrami results')
parser.add_argument('--from', dest="former_results", nargs='?', default="results.pkl", help='reuse old results from file')
args = parser.parse_args()
options = vars(args)

with open(options["former_results"], 'rb') as f:
    resultdict = pickle.load(f)

orders = [1,2,3,4,5,6,7,8,9,10]
levels = 10

type1 = (["hdg_hodc_ndofs","hdg_ndofs","dg_ndofs","cggp_ndofs","cg_ndofs"],1)
type2 = (["hdg_hodc_global_ndofs","hdg_global_ndofs","dg_global_ndofs","cggp_global_ndofs","cg_global_ndofs"],1)
type3 = (["hdg_hodc_nze","hdg_nze","dg_nze","cggp_nze","cg_nze"],1000)
type4 = (["hdg_hodc_nze_cond","hdg_nze_cond","dg_nze_cond","cggp_nze_cond","cg_nze_cond"],1000)

for types,factor in [type1,type2,type3,type4]:
    for order in orders:
        order_hasdata = False
        for i in range(levels):
            if (i,order) in resultdict:
                order_hasdata = True
                break
        if order_hasdata:
            print("\\midrule")
        else:
            continue
        for i in range(levels):
            if (i==0):
                pbegin = "{:2}".format(order)
            else:
                pbegin = "  "
            print("{} ".format(pbegin),end="")
            if (i,order) in resultdict:
                for othertype in types:
                    if othertype in resultdict[(i,order)]:
                        print("& {:8}".format(int(resultdict[(i,order)][othertype]/factor)),end="")
                print(" \\\\")
    print("\\bottomrule\n\n")

            # for othertype in ["hdg_nze","dg_ndofs","cggp_ndofs","cg_ndofs"]:
            #     if othertype in resultdict[(i,order)]:
            #         print("& {:8}".format(resultdict[(i,order)][othertype]),end="")

            # for othertype in ["nze","dg_nze","cg_gp_nze","cg_nze","cg_nze_cond"]:
            #     if othertype in resultdict[(i,order)]:
            #         print("& {:6}".format(int(resultdict[(i,order)][othertype]/1000)),end="")

            # for othertype in ["global_ndofs","total_ndofs","dg_ndofs","cg_ndofs","cg_global_ndofs"]:
            #     if othertype in resultdict[(i,order)]:
            #         print("& {:5.2f}".format(resultdict[(i,order)][othertype]/resultdict[(i,order)]["cg_global_ndofs"]),end="")
            # for othertype in ["cg_global_ndofs"]:
            #     if othertype in resultdict[(i,order)]:
            #         print("& {:8}".format(resultdict[(i,order)][othertype]),end="")
            # for othertype in ["nze","dg_nze","cg_gp_nze","cg_nze","cg_nze_cond"]:
            #     if othertype in resultdict[(i,order)]:
            #         print("& {:5.2f}".format(resultdict[(i,order)][othertype]/resultdict[(i,order)]["cg_nze_cond"]),end="")
            # for othertype in ["cg_nze_cond"]:
            #     if othertype in resultdict[(i,order)]:
            #         print("& {:6}".format(int(resultdict[(i,order)][othertype]/1000)),end="")



# for order in orders:
#     print("\\midrule")
#     i = 3
#     pbegin = "{}".format(order)
#     print("{} ".format(pbegin),end="")
#     if (i,order) in resultdict:
#         for errortype in ["total_ndofs","global_ndofs","dg_ndofs","nze","dg_nze"]:
#             print("& {:9.0f} K ".format(resultdict[(i,order)][errortype]/1000),end="")
#         print(" \\\\")
# print("\\bottomrule")
