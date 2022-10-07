__ngsolve_required__ = "6.2.2202-0"

def check_if_ngsolve_newer_than(ngsolve_version_required):
    """
    Check for compatibility of ngsolve version
    """
    import ngsolve
    import re
    ngsver = [0 for i in range(4)]
    res = re.split(r'-', ngsolve.__version__)
    if len(res) > 1:
        mmp, ngsver[3], h = res
    else:
        mmp = res[0]
        ngsver[3] = "0"
    ngsver[0], ngsver[1], ngsver[2] = re.split(re.compile('[.bv]'), mmp)

    ngsver_required = [0 for i in range(4)]
    mmp, ngsver_required[3] =  re.split(re.compile('[-bv]'), ngsolve_version_required)
    ngsver_required[0:3] = re.split(re.compile('[.bv]'), mmp)

    ret = True
    for i in range(4):
        if int(ngsver[i]) < int(ngsver_required[i]):
            ret = False
            print("""
####################################################
                     WARNING!""")
            print("Your NGSolve-Version (" + ngsolve.__version__+")")
            print("is older than required (>=" + ngsolve_version_required +").",end="")
            print("""
####################################################
""")
            return False
        elif int(ngsver[i]) > int(ngsver_required[i]):
            break
    return ret
