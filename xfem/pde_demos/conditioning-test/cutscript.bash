#!/bin/bash
echo "# ---"  > out.stdmesh.cond.A
echo "# ---"  > out.stdmesh.cond.BDDCinvA
echo "# ---"  > out.stdmesh.cond.BBDDCinvA
echo "# ---"  > out.stdmesh.cond.BJinvA
echo "# ---"  > out.stdmesh.cond.JinvA
echo "# ---"  > out.stdmesh.its.J
echo "# ---"  > out.stdmesh.its.BJ
echo "# ---"  > out.stdmesh.its.BDDC
echo "# ---"  > out.stdmesh.its.BBDDC
echo "# ---"  > out.stdmesh.err

LAM=2
ANEG=2

#***A***
# POS=370500
#***B***
POS=470500
#***C***
# POS=370500

#***A***
# LOW=3705
# UP=371000
# POS=$LOW
#***B***
LOW=1
UP=1700000
LAM=$LOW
#***C***
# LOW=1
# UP=2000000
# ANEG=$LOW

#***A***
# sed "s/calcaxis/0.37+(\$0)*0.00001/g" stdmesh.gpl.template > stdmesh.gpl.temp
# sed "s/hiddenlogscalex/\#/g" stdmesh.gpl.temp > stdmesh.gpl
#***B***
sed "s/calcaxis/2**(\$0)/g" stdmesh.gpl.template > stdmesh.gpl.temp
sed "s/hiddenlogscalex//g" stdmesh.gpl.temp > stdmesh.gpl
#***C***
# sed "s/calcaxis/2**(\$0-10)/g" stdmesh.gpl.template > stdmesh.gpl.temp
# sed "s/hiddenlogscalex//g" stdmesh.gpl.temp > stdmesh.gpl


#***A***
#while [  $POS -lt $UP ]; do
#***B***
while [  $LAM -lt $UP ]; do
#***C***
#while [  $ANEG -lt $UP ]; do
    echo "POS = $POS"
    echo "LAM = $LAM"
    echo "ANEG = $ANEG"
    sed "s/var/$POS/g" stdmesh_test.pde.template > stdmesh_test.pde.temp
    sed "s/setlam/$LAM/g" stdmesh_test.pde.temp > stdmesh_test.pde.temp2
    sed "s/ANEG/$ANEG/g" stdmesh_test.pde.temp2 > stdmesh_test.pde
    ngs stdmesh_test.pde > out.stdmesh
    # # echo "numproc quit npquit" >> stdmesh_test.pde
    # netgen -pdefile=stdmesh_test.pde -solve=$POS | tee out.stdmesh
    cat out.stdmesh | grep "condition number:" | tail -n 2 | head -n 1 >> out.stdmesh.cond.A
    cat out.stdmesh | grep " Condition" | tail -n 4 | head -n 1 >> out.stdmesh.cond.BDDCinvA
    cat out.stdmesh | grep " Condition" | tail -n 4 |head -n 2 | tail -n 1 >> out.stdmesh.cond.BBDDCinvA
    cat out.stdmesh | grep " Condition" | tail -n 4 |tail -n 2 | head -n 1 >> out.stdmesh.cond.BJinvA
    cat out.stdmesh | grep " Condition" | tail -n 4 |tail -n 1 >> out.stdmesh.cond.JinvA
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |head -n 1  >> out.stdmesh.its.BDDC
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |head -n 2 | tail -n 1 >> out.stdmesh.its.BBDDC
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |tail -n 2 | head -n 1 >> out.stdmesh.its.BJ
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |tail -n 1  >> out.stdmesh.its.J
    cat out.stdmesh | grep "l2_n" -A 2 | tail -n 1 >> out.stdmesh.err
    cat out.stdmesh | grep "condition number:" | tail -n 2 |head -n 1
    cat out.stdmesh | grep " Condition" | tail -n 4 |head -n 1
    cat out.stdmesh | grep " Condition" | tail -n 4 |head -n 2 | tail -n 1
    cat out.stdmesh | grep " Condition" | tail -n 4 |tail -n 2 | head -n 1
    cat out.stdmesh | grep " Condition" | tail -n 4 |tail -n 1
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |head -n 1 
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |head -n 2 | tail -n 1
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |tail -n 2 | head -n 1
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |tail -n 1 
    cat out.stdmesh | grep "l2_n" -A 2 | tail -n 1
#    cat out.stdmesh | grep "l2_n" -A 1 | tail -n 1 >> out.stdmesh.err
#    cat out.stdmesh | grep "l2_n" -A 1 | tail -n 1
#***A***
#let POS=POS+10
#***B***
let LAM=LAM*2
#***A***
#let ANEG=ANEG*2
done

gnuplot stdmesh.gpl
