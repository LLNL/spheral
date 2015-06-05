#-*-shell-script-*-
# Gnuplot script to plot the 1-D planar Noh problem convergence results.
set title "1-D Planar Noh Mass density convergence"
set xlabel "N"
set ylabel "Error"

# Plot the results without the tensile correction.
plot "Noh-planar-convergence-test-nperh=2.01-XSPH=False.txt" using 1:2 ps 2 t "L1, SPH"
replot "Noh-planar-convergence-test-nperh=2.01-XSPH=False.txt" using 1:3 ps 2 t "L2, SPH"
replot "Noh-planar-convergence-test-nperh=2.01-XSPH=False.txt" using 1:4 ps 2 t "Linf, SPH"

# Plot the results with the tensile correction.
replot "Noh-planar-convergence-test-nperh=2.01-XSPH=True.txt" using 1:2 ps 2 t "L1, XSPH"
replot "Noh-planar-convergence-test-nperh=2.01-XSPH=True.txt" using 1:3 ps 2 t "L2, XSPH"
replot "Noh-planar-convergence-test-nperh=2.01-XSPH=True.txt" using 1:4 ps 2 t "Linf, XSPH"

set xrange [20:2000]
set logscale xy
set key bottom left
replot

# Fit the slopes for the SPH case.
fL10(x) = A10 * x**m10
fit fL10(x) "Noh-planar-convergence-test-nperh=2.01-XSPH=False.txt" using 1:2 via A10, m10
replot fL10(x) notitle # t "L1 fit, SPH"

fL20(x) = 10.0**A20 * x**m20
fit fL20(x) "Noh-planar-convergence-test-nperh=2.01-XSPH=False.txt" using 1:3 via A20, m20
replot fL20(x) notitle # t "L2 fit, SPH"

fLinf0(x) = 10.0**Ainf0 * x**minf0
fit fLinf0(x) "Noh-planar-convergence-test-nperh=2.01-XSPH=False.txt" using 1:4 via Ainf0, minf0
replot fLinf0(x) notitle  # t "Linf fit, SPH"

# Fit the slopes for the case with XSPH.
fL11(x) = 10.0**A11 * x**m11
fit fL11(x) "Noh-planar-convergence-test-nperh=2.01-XSPH=True.txt" using 1:2 via A11, m11
replot fL11(x) notitle # t "L1 fit, XSPH"

fL21(x) = 10.0**A21 * x**m21
fit fL21(x) "Noh-planar-convergence-test-nperh=2.01-XSPH=True.txt" using 1:3 via A21, m21
replot fL21(x) notitle # t "L2 fit, XSPH"

fLinf1(x) = 10.0**Ainf1 * x**minf1
fit fLinf1(x) "Noh-planar-convergence-test-nperh=2.01-XSPH=True.txt" using 1:4 via Ainf1, minf1
replot fLinf1(x) notitle # t "Linf fit, XSPH"

# Summarize the fitting parameters.
print "================================================================================"
print "Found fitting parameters of:"
print "SPH                    :  L1 order of convergence = ", m10
print "                          L2 order of convergence = ", m20
print "                        Linf order of convergence = ", minf0
print "XSPH                   :  L1 order of convergence = ", m11
print "                          L2 order of convergence = ", m21
print "                        Linf order of convergence = ", minf1
