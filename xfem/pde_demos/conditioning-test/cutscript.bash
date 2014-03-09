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
LOW=3000
UP=7000
# LOW=1
# UP=4
EXP=$LOW
while [  $EXP -lt $UP ]; do
    echo "$LOW / $EXP / $UP"
    sed "s/var/$EXP/g" stdmesh_test.pde.template > stdmesh_test.pde
    ngs stdmesh_test.pde > out.stdmesh
    # sed "s/var/5000/g" stdmesh_test.pde.template > stdmesh_test.pde
    # # echo "numproc quit npquit" >> stdmesh_test.pde
    # netgen -pdefile=stdmesh_test.pde -solve=$EXP | tee out.stdmesh
    cat out.stdmesh | grep "condition number:" | tail -n 2 | head -n 1 >> out.stdmesh.cond.A
    cat out.stdmesh | grep " Condition" | tail -n 4 | head -n 1 >> out.stdmesh.cond.BDDCinvA
    cat out.stdmesh | grep " Condition" | tail -n 4 |head -n 2 | tail -n 1 >> out.stdmesh.cond.BBDDCinvA
    cat out.stdmesh | grep " Condition" | tail -n 4 |tail -n 2 | head -n 1 >> out.stdmesh.cond.BJinvA
    cat out.stdmesh | grep " Condition" | tail -n 4 |tail -n 1 >> out.stdmesh.cond.JinvA
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |head -n 1  >> out.stdmesh.its.BDDC
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |head -n 2 | tail -n 1 >> out.stdmesh.its.BBDDC
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |tail -n 2 | head -n 1 >> out.stdmesh.its.BJ
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |tail -n 1  >> out.stdmesh.its.J
    cat out.stdmesh | grep "condition number:" | tail -n 2 |head -n 1
    cat out.stdmesh | grep " Condition" | tail -n 4 |head -n 1
    cat out.stdmesh | grep " Condition" | tail -n 4 |head -n 2 | tail -n 1
    cat out.stdmesh | grep " Condition" | tail -n 4 |tail -n 2 | head -n 1
    cat out.stdmesh | grep " Condition" | tail -n 4 |tail -n 1
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |head -n 1 
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |head -n 2 | tail -n 1
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |tail -n 2 | head -n 1
    cat out.stdmesh | grep "Iterations:" | tail -n 4 |tail -n 1 
#    cat out.stdmesh | grep "l2_n" -A 1 | tail -n 1 >> out.stdmesh.err
#    cat out.stdmesh | grep "Condition CinvA"
#    cat out.stdmesh | grep "Condition A"
#    cat out.stdmesh | grep Iterations
#    cat out.stdmesh | grep "l2_n" -A 1 | tail -n 1
let EXP=EXP+1
done

gnuplot stdmesh.gpl
