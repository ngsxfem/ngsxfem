from math import log
import pickle
with open('results.pkl', 'rb') as f:
    resultdict = pickle.load(f)

orders = [1,2,3,4,5,6,7,8,9,10]
levels = 10

#for order in orders:
#    for i in range(10):
#        if (i,order) in resultdict:
#            print("({},{}): {}".format(i,order,resultdict[(i,order)]))


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
            for errortype in ["maxdist","l2err","h1normerr","h1tangerr"]:
                if (i==0):
                    eoc = " ---"
                else:
                    eoc = "{:>0.1f}".format(log(resultdict[(i-1,order)][errortype] / resultdict[(i,order)][errortype])/log(2.0))
                eoc = "{:>4}".format(eoc)
                print("& \\num{{{:0.5e}}} & ({}) ".format(resultdict[(i,order)][errortype],eoc),end="")
            for othertype in ["global_ndofs","total_ndofs","nze","dg_ndofs","dg_nze","cg_ndofs","cg_nze","cg_gp_nze"]:
                if othertype in resultdict[(i,order)]:
                    print("& {:10}".format(resultdict[(i,order)][othertype]),end="")
            print(" \\\\")
print("\\bottomrule")



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
