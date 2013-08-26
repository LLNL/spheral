#-*-shell-script-*-
# Gnuplot script to plot the 1-D planar Noh problem convergence results.
set title "1-D Planar Noh pressure convergence"
set xlabel "N"
set ylabel "Error"

# Plot the results without the tensile correction.
plot "Noh-planar-convergence-test.txt-SPH" using 1:5 t "L1, SPH"
replot "Noh-planar-convergence-test.txt-SPH" using 1:6 t "L2, SPH"
replot "Noh-planar-convergence-test.txt-SPH" using 1:7 t "Linf, SPH"

# Plot the results with the tensile correction.
replot "Noh-planar-convergence-test.txt-XSPH" using 1:5 t "L1, XSPH"
replot "Noh-planar-convergence-test.txt-XSPH" using 1:6 t "L2, XSPH"
replot "Noh-planar-convergence-test.txt-XSPH" using 1:7 t "Linf, XSPH"

set xrange [20:2000]
set logscale xy
replot

# Fit the slopes for the SPH case.
logfL10(x) = A10 + m10*x
fit logfL10(x) "Noh-planar-convergence-test.txt-SPH" using (log10($1)):(log10($5)) via A10, m10
fL10(x) = 10.0**A10 * x**m10
replot fL10(x) t "L1 fit, SPH"

logfL20(x) = A20 + m20*x
fit logfL20(x) "Noh-planar-convergence-test.txt-SPH" using (log10($1)):(log10($6)) via A20, m20
fL20(x) = 10.0**A20 * x**m20
replot fL20(x) t "L2 fit, SPH"

logfLinf0(x) = Ainf0 + minf0*x
fit logfLinf0(x) "Noh-planar-convergence-test.txt-SPH" using (log10($1)):(log10($7)) via Ainf0, minf0
fLinf0(x) = 10.0**Ainf0 * x**minf0
replot fLinf0(x) t "Linf fit, SPH"

# Fit the slopes for the case with tensile correction.
logfL11(x) = A11 + m11*x
fit logfL11(x) "Noh-planar-convergence-test.txt-XSPH" using (log10($1)):(log10($5)) via A11, m11
fL11(x) = 10.0**A11 * x**m11
replot fL11(x) t "L1 fit, XSPH"

logfL21(x) = A21 + m21*x
fit logfL21(x) "Noh-planar-convergence-test.txt-XSPH" using (log10($1)):(log10($6)) via A21, m21
fL21(x) = 10.0**A21 * x**m21
replot fL21(x) t "L2 fit, XSPH"

logfLinf1(x) = Ainf1 + minf1*x
fit logfLinf1(x) "Noh-planar-convergence-test.txt-XSPH" using (log10($1)):(log10($7)) via Ainf1, minf1
fLinf1(x) = 10.0**Ainf1 * x**minf1
replot fLinf1(x) t "Linf fit, XSPH"

# Summarize the fitting parameters.
print "================================================================================"
print "Found fitting parameters of:"
print "SPH                    :  L1 order of convergence = ", m10
print "                          L2 order of convergence = ", m20
print "                        Linf order of convergence = ", minf0
print "XSPH                   :  L1 order of convergence = ", m11
print "                          L2 order of convergence = ", m21
print "                        Linf order of convergence = ", minf1
