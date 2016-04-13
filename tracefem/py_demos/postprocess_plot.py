from math import log
import pickle
import matplotlib.pyplot as plt
import sys
import numpy as np
import argparse        

parser = argparse.ArgumentParser(description='Postprocess laplace beltrami results')
parser.add_argument('--from', dest="former_results", nargs='?', default="results.pkl", help='reuse old results from file')
args = parser.parse_args()
options = vars(args)

with open(options["former_results"], 'rb') as f:
    resultdict = pickle.load(f)

orders = [1,2,3,4,5,6,7,8,9,10]
levels = 10

types = ["l2err","maxdist","h1normerr","h1tangerr"]
plots = {}

for order in orders:
    order_hasdata = False
    for i in range(levels):
        if (i,order) in resultdict:
            order_hasdata = True
            break
    if not order_hasdata:
        continue

    for curtype in types:
        plots[(curtype,order)] = ([],[])

    for i in range(levels):
        if (i,order) in resultdict:
            for curtype in types:
                if curtype in resultdict[(i,order)]:
                    plots[(curtype,order)][0].append(i)
                    plots[(curtype,order)][1].append(resultdict[(i,order)][curtype])
    
cnt = 0
legendlabels = []

import warnings
warnings.simplefilter('ignore', np.RankWarning)

for curtype in types:
    plt.figure(cnt)
    plt.yscale('log')
    plt.xlabel("level")
    cnt += 1
    legendlabels = []
    for order in orders:
        if (curtype,order) in plots:
            plt.plot(plots[(curtype,order)][0],plots[(curtype,order)][1],'8-')
            legendlabels.append(curtype+"(k="+str(order)+")")

            slope, intercept = np.polyfit(plots[(curtype,order)][0], 1.0/np.log(0.5) * np.log(plots[(curtype,order)][1]), 1)
            print("average slope for {}, k={} : {}".format(curtype,order,slope))
    plt.legend(legendlabels)

plt.ion()
plt.show()
    
# for curtype in types:
#     for i in range(levels):

# l2conv = [ log(l2err[i]/l2err[i-1])/log(0.5) for i in range(1,len(l2err))]
# h1tconv = [ log(h1tangerr[i]/h1tangerr[i-1])/log(0.5) for i in range(1,len(l2err))]
# h1nconv = [ log(h1normerr[i]/h1normerr[i-1])/log(0.5) for i in range(1,len(l2err))]
# # maxconv = [ log(maxerr[i]/maxerr[i-1])/log(0.5) for i in range(1,len(maxerr))]
# geomconv = [ log(geomerr[i]/geomerr[i-1])/log(0.5) for i in range(1,len(geomerr))]
# print ("l2err convergence orders (eoc):", l2conv)
# print ("h1 tang err conv. orders (eoc):", h1tconv)
# print ("h1 norm err conv. orders (eoc):", h1nconv)
# print ("geom. convergence orders (eoc):", geomconv)

if (not hasattr(sys,'ps1')):
    input("<press enter to quit>")
