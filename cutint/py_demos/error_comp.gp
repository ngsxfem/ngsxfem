set terminal pdf
set output "comp.pdf"

set format y "%g"

set xlabel "Order q"
set ylabel "sqrt(error[POS]**2 + error[NEG]**2)"

set logscale y

set title "Case lsetvals = [-0.18687,0.324987, 0.765764,0.48983]"

plot "errors.dat" index 0 title "Our method", "errors_saye.dat" index 0 title "Saye"

set title "Case lsetvals = [0.765764,0.324987, -0.18687, -0.48983]"
plot "errors.dat" index 1 title "Our method", "errors_saye.dat" index 1 title "Saye"

set title "Case lsetvals = [1,2/3,-1,-2/3]"
plot "errors.dat" index 2 title "Our method", "errors_saye.dat" index 2 title "Saye"

set title "Case lsetvals = [3,-1,1,-1.023123]"
plot "errors.dat" index 3 title "Our method", "errors_saye.dat" index 3 title "Saye"
