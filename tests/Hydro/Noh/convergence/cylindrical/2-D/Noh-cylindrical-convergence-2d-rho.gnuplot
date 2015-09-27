#-*-shell-script-*-
# Gnuplot script to plot the 2-D cylindrical Noh problem convergence results.
set title "2-D Cylindrical Noh Mass density convergence"
set xlabel "N"
set ylabel "Error"

# Plot the SPH results.
plot "Noh-cylindrical-convergence-test.txt-xsph" using 1:2 t "L1, XSPH"
replot "Noh-cylindrical-convergence-test.txt-xsph" using 1:3 t "L2, XSPH"
replot "Noh-cylindrical-convergence-test.txt-xsph" using 1:4 t "Linf, XSPH"

# Plot the ASPH results.
replot "Noh-cylindrical-convergence-test.txt-asph" using 1:2 t "L1, ASPH"
replot "Noh-cylindrical-convergence-test.txt-asph" using 1:3 t "L2, ASPH"
replot "Noh-cylindrical-convergence-test.txt-asph" using 1:4 t "Linf, ASPH"

#set xrange [20:2000]
set logscale xy
replot

# Fit the slopes for the SPH case.
fL10(x) = A10 * x**m10
fit fL10(x) "Noh-cylindrical-convergence-test.txt-xsph" using 1:2 via A10, m10
replot fL10(x) t "L1 fit, XSPH"
fL20(x) = A20 * x**m20
fit fL20(x) "Noh-cylindrical-convergence-test.txt-xsph" using 1:3 via A20, m20
replot fL20(x) t "L2 fit, XSPH"
fLinf0(x) = Ainf0 * x**minf0
fit fLinf0(x) "Noh-cylindrical-convergence-test.txt-xsph" using 1:4 via Ainf0, minf0
replot fLinf0(x) t "Linf fit, XSPH"

# Fit the slopes for the ASPH case.
fL11(x) = A11 * x**m11
fit fL11(x) "Noh-cylindrical-convergence-test.txt-asph" using 1:2 via A11, m11
replot fL11(x) t "L1 fit, ASPH"
fL21(x) = A21 * x**m21
fit fL21(x) "Noh-cylindrical-convergence-test.txt-asph" using 1:3 via A21, m21
replot fL21(x) t "L2 fit, ASPH"
fLinf1(x) = Ainf1 * x**minf1
fit fLinf1(x) "Noh-cylindrical-convergence-test.txt-asph" using 1:4 via Ainf1, minf1
replot fLinf1(x) t "Linf fit, ASPH"

# Summarize the fitting parameters.
print "================================================================================"
print "Found fitting parameters of:"
print "XSPH                   :  L1 order of convergence = ", m10
print "                          L2 order of convergence = ", m20
print "                        Linf order of convergence = ", minf0
print "ASPH                   :  L1 order of convergence = ", m11
print "                          L2 order of convergence = ", m21
print "                        Linf order of convergence = ", minf1
