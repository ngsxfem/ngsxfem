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

#PDE_TEMPLATE="stdmesh_test.pde.template"
PDE_TEMPLATE="roundedsquare.pde.template"

LAM=2
ANEG=1024
STEPPOS=2
#***A***
POS=0
#***B***
# POS=470500
#***C***
# POS=370500

#***A***
LOW=0
UP=8000
POS=$LOW
#***B***
# LOW=1
# UP=6500000
# LAM=$LOW
#***C***
# LOW=1
# UP=2000000
# ANEG=$LOW

#***A***
sed "s/calcaxis/($LOW\/1e6)+(\$0)*($STEPPOS\/1e6)/g" stdmesh.gpl.template > stdmesh.gpl.temp
sed "s/hiddenlogscalex/\#/g" stdmesh.gpl.temp > stdmesh.gpl
#***B***
# sed "s/calcaxis/2**(\$0)/g" stdmesh.gpl.template > stdmesh.gpl.temp
# sed "s/hiddenlogscalex//g" stdmesh.gpl.temp > stdmesh.gpl
#***C***
# sed "s/calcaxis/2**(\$0-10)/g" stdmesh.gpl.template > stdmesh.gpl.temp
# sed "s/hiddenlogscalex//g" stdmesh.gpl.temp > stdmesh.gpl


#***A***
while [  $POS -lt $UP ]; do
#***B***
# while [  $LAM -lt $UP ]; do
#***C***
# while [  $ANEG -lt $UP ]; do
    echo "POS = $POS"
    echo "LAM = $LAM"
    echo "ANEG = $ANEG / 1024"
    sed "s/var/$POS/g" $PDE_TEMPLATE > temporary.pde.temp
    sed "s/setlam/$LAM/g" temporary.pde.temp > temporary.pde.temp2
    sed "s/ANEG/$ANEG/g" temporary.pde.temp2 > temporary.pde
    ngs temporary.pde > out.stdmesh
    # # echo "numproc quit npquit" >> temporary.pde
    # netgen -pdefile=temporary.pde -solve=$POS | tee out.stdmesh
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
    let POS=POS+STEPPOS
#***B***
# let LAM=LAM*2
#***C***
# let ANEG=ANEG*2
done

gnuplot stdmesh.gpl
