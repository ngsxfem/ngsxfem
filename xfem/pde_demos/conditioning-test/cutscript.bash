#!/bin/bash
echo "# ---"  > out.stdmesh.cond.CinvA
echo "# ---"  > out.stdmesh.cond.A
echo "# ---"  > out.stdmesh.err
echo "# ---"  > out.stdmesh.its
LOW=300
UP=700
EXP=$LOW
while [  $EXP -lt $UP ]; do
    echo "$LOW / $EXP / $UP"
    sed "s/var/$EXP/g" stdmesh_test.pde.template > stdmesh_test.pde
    ngs stdmesh_test.pde > out.stdmesh
    cat out.stdmesh | grep "condition number:" | head -n 1 >> out.stdmesh.cond.A
    cat out.stdmesh | grep "condition number:" | tail -n 1 >> out.stdmesh.cond.CinvA
    cat out.stdmesh | grep "Iterations:" >> out.stdmesh.its
    cat out.stdmesh | grep "condition number:" | head -n 1
    cat out.stdmesh | grep "condition number:" | tail -n 1
    cat out.stdmesh | grep "Iterations:"
#    cat out.stdmesh | grep "l2_n" -A 1 | tail -n 1 >> out.stdmesh.err
#    cat out.stdmesh | grep "Condition CinvA"
#    cat out.stdmesh | grep "Condition A"
#    cat out.stdmesh | grep Iterations
#    cat out.stdmesh | grep "l2_n" -A 1 | tail -n 1
let EXP=EXP+1
done

gnuplot stdmesh.gpl
