#-------------------------------------------------------------------------------
# Python script to plot the convergence properties of the 1-D planar Sedov problem
# results.
#-------------------------------------------------------------------------------
import Gnuplot
import time
import os

printcmds = '''
set term latex
set output %(filename)s
'''

def plotscaling(p, opts):

    # Blast any existing fit output.
    if os.path.exists("fit.log"):
        os.remove("fit.log")

    #fit [2.5:10] logfL10(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using (log10($1)):(log10($%(L1)s)) via A10, m10

    # Fit the error scaling paramters.
    p('''
    # Fit the slopes for the standard case.
    logfL10(x) = A10 + m10*x
    fL10(x) = 10.0**A10 * x**m10
    #fit [2.5:10] logfL10(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using (log10($1)):(log10($%(L1)s)) via A10, m10
    fit logfL10(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using (log10($1)):(log10($%(L1)s)) via A10, m10

    logfL20(x) = A20 + m20*x
    fL20(x) = 10.0**A20 * x**m20
    #fit [2.5:10] logfL20(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using (log10($1)):(log10($%(L2)s)) via A20, m20
    fit logfL20(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using (log10($1)):(log10($%(L2)s)) via A20, m20

    logfLinf0(x) = Ainf0 + minf0*x
    fLinf0(x) = 10.0**Ainf0 * x**minf0
    fit #[2.5:10] logfLinf0(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using (log10($1)):(log10($%(Linf)s)) via Ainf0, minf0
    fit logfLinf0(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using (log10($1)):(log10($%(Linf)s)) via Ainf0, minf0

    # Fit the slopes for the case using compatible differencing.
    logfL11(x) = A11 + m11*x
    fL11(x) = 10.0**A11 * x**m11
    #fit [2.5:10] logfL11(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using (log10($1)):(log10($%(L1)s)) via A11, m11
    fit logfL11(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using (log10($1)):(log10($%(L1)s)) via A11, m11

    logfL21(x) = A21 + m21*x
    fL21(x) = 10.0**A21 * x**m21
    #fit [2.5:10] logfL21(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using (log10($1)):(log10($%(L2)s)) via A21, m21
    fit logfL21(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using (log10($1)):(log10($%(L2)s)) via A21, m21

    logfLinf1(x) = Ainf1 + minf1*x
    fLinf1(x) = 10.0**Ainf1 * x**minf1
    #fit [2.5:10] logfLinf1(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using (log10($1)):(log10($%(Linf)s)) via Ainf1, minf1
    fit logfLinf1(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using (log10($1)):(log10($%(Linf)s)) via Ainf1, minf1
    ''' % opts)

    # Parse the fitted results from the fit log.
    time.sleep(2)
    vars = [("m10", "sig10"),
            ("m20", "sig20"),
            ("minf0", "siginf0"),
            ("m11", "sig11"),
            ("m21", "sig21"),
            ("minf1", "siginf1")]
    mvars = [x[0] for x in vars]
    f = open("fit.log", "r")
    for ln in f:
        if "+/-" in ln:
            crap = ln.split()
            if crap[0] in mvars:
                i = mvars.index(crap[0])
                exec("%s = %s" % (vars[i][0], crap[2]))
                exec("%s = %s" % (vars[i][1], crap[4]))

    # Output the fitting results.
    print "standard               :  L1 order of convergence = %g +/- %g" % (m10, sig10)
    print "                          L2 order of convergence = %g +/- %g" % (m20, sig20)
    print "                        Linf order of convergence = %g +/- %g" % (minf0, siginf0)
    print "compatible             :  L1 order of convergence = %g +/- %g" % (m11, sig11)
    print "                          L2 order of convergence = %g +/- %g" % (m21, sig21)
    print "                        Linf order of convergence = %g +/- %g" % (minf1, siginf1)

    # Buid the Gnuplot data objects.
    dL1s = Gnuplot.File("Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" % opts,
                        using = eval("(1, %(L1)s)" % opts),
                        with = "points ps 2",
                        title = "L_1, standard")
    dL2s = Gnuplot.File("Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" % opts,
                        using = eval("(1, %(L2)s)" % opts),
                        with = "points ps 2",
                        title = "L_2, standard")
    dLinfs = Gnuplot.File("Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" % opts,
                          using = eval("(1, %(Linf)s)" % opts),
                          with = "points ps 2",
                          title = "L_{/Symbol \245}, standard")
    dL1c = Gnuplot.File("Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" % opts,
                        using = eval("(1, %(L1)s)" % opts),
                        with = "points ps 2",
                        title = "L_1, compatible")
    dL2c = Gnuplot.File("Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" % opts,
                        using = eval("(1, %(L2)s)" % opts),
                        with = "points ps 2",
                        title = "L_2, compatible")
    dLinfc = Gnuplot.File("Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" % opts,
                          using = eval("(1, %(Linf)s)" % opts),
                          with = "points ps 2",
                          title = "L_{/Symbol \245}, compatible")

    # Plot the measured error points.
    p.plot(dL1s)
    p.replot(dL2s)
    p.replot(dLinfs)
    p.replot(dL1c)
    p.replot(dL2c)
    p.replot(dLinfc)

    # Plot the fitted scalings.
    p.replot(Gnuplot.Func("fL10(x)", with = "lines lt 1", title = ""))
    p.replot(Gnuplot.Func("fL20(x)", with = "lines lt 1", title = ""))
    p.replot(Gnuplot.Func("fLinf0(x)", with = "lines lt 1", title = ""))
    p.replot(Gnuplot.Func("fL11(x)", with = "lines lt 1", title = ""))
    p.replot(Gnuplot.Func("fL21(x)", with = "lines lt 1", title = ""))
    p.replot(Gnuplot.Func("fLinf1(x)", with = "lines lt 1", title = ""))

    # Labeling and such.
    #p.xlabel("n")
    #p.ylabel("Error")
    p("set xrange [40:40000]")
    p("set yrange [%(ymin)s:%(ymax)s]" % opts)
    p("set logscale xy")
    p("set key graph 0.48, graph 0.35")
    p.refresh()

    return

plotcmds = '''
# Gnuplot script to plot the 1-D planar Sedov problem convergence results.
set xlabel "N"
set ylabel "Error"

## # Plot the SPH results
## plot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=False.txt" using 1:%(L1)s ps 2 t "L1, SPH"
## replot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=False.txt" using 1:%(L2)s ps 2 t "L2, SPH"
## replot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=False.txt" using 1:%(Linf)s ps 2 t "Linf, SPH"

## # Plot the XSPH results
## replot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True.txt" using 1:%(L1)s ps 2 t "L1, XSPH"
## replot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True.txt" using 1:%(L2)s ps 2 t "L2, XSPH"
## replot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True.txt" using 1:%(Linf)s ps 2 t "Linf, XSPH"

# Plot the non-compatible results
plot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using 1:%(L1)s ps 2 t "L1, standard"
replot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using 1:%(L2)s ps 2 t "L2, standard"
replot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using 1:%(Linf)s ps 2 t "Linf, standard"

# Plot the compatible results
replot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using 1:%(L1)s ps 2 t "L1, compatible"
replot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using 1:%(L2)s ps 2 t "L2, compatible"
replot "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using 1:%(Linf)s ps 2 t "Linf, compatible"

set xrange [20:1200]
set logscale xy
set key bottom left
replot

# Fit the slopes for the SPH case.
logfL10(x) = A10 + m10*x
fL10(x) = 10.0**A10 * x**m10
fit [*:*] logfL10(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using (log10($1)):(log10($%(L1)s)) via A10, m10
replot fL10(x) notitle # t "L1 fit, SPH"

logfL20(x) = A20 + m20*x
fL20(x) = 10.0**A20 * x**m20
fit [*:*] logfL20(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using (log10($1)):(log10($%(L2)s)) via A20, m20
replot fL20(x) notitle # t "L2 fit, SPH"

logfLinf0(x) = Ainf0 + minf0*x
fLinf0(x) = 10.0**Ainf0 * x**minf0
fit [*:*] logfLinf0(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=False.txt" using (log10($1)):(log10($%(Linf)s)) via Ainf0, minf0
replot fLinf0(x) notitle  # t "Linf fit, SPH"

# Fit the slopes for the case with XSPH.
logfL11(x) = A11 + m11*x
fL11(x) = 10.0**A11 * x**m11
fit [*:*] logfL11(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using (log10($1)):(log10($%(L1)s)) via A11, m11
replot fL11(x) notitle # t "L1 fit, XSPH"

logfL21(x) = A21 + m21*x
fL21(x) = 10.0**A21 * x**m21
fit [*:*] logfL21(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using (log10($1)):(log10($%(L2)s)) via A21, m21
replot fL21(x) notitle # t "L2 fit, XSPH"

logfLinf1(x) = Ainf1 + minf1*x
fLinf1(x) = 10.0**Ainf1 * x**minf1
fit [*:*] logfLinf1(x) "Sedov-planar-convergence-test-nperh=%(nperh)s-XSPH=True-compatibleEnergy=True.txt" using (log10($1)):(log10($%(Linf)s)) via Ainf1, minf1
replot fLinf1(x) notitle # t "Linf fit, XSPH"
'''

summarycmds = '''
# Summarize the fitting parameters.
print "standard               :  L1 order of convergence = ", m10
print "                          L2 order of convergence = ", m20
print "                        Linf order of convergence = ", minf0
print "compatible             :  L1 order of convergence = ", m11
print "                          L2 order of convergence = ", m21
print "                        Linf order of convergence = ", minf1
'''

## #-------------------------------------------------------------------------------
## # Plot the nperh = 1.01 density results.
## p101 = Gnuplot.Gnuplot()
## p101.title("Mass density, nperh=1.01")
## d101 = {"nperh": "1.01",
##         "L1"   : "2",
##         "L2"   : "3",
##         "Linf" : "4"}
## p101(plotcmds % d101)

## #-------------------------------------------------------------------------------
## # Plot the nperh = 1.51 density results.
## p151 = Gnuplot.Gnuplot()
## p151.title("Mass density, nperh=1.51")
## d151 = {"nperh": "1.51",
##         "L1"   : "2",
##         "L2"   : "3",
##         "Linf" : "4"}
## p151(plotcmds % d151)

## #-------------------------------------------------------------------------------
## # Plot the nperh = 2.01 density results.
## p201 = Gnuplot.Gnuplot()
## p201.title("Mass density, nperh=2.01")
## d201 = {"nperh": "2.01",
##         "L1"   : "2",
##         "L2"   : "3",
##         "Linf" : "4"}
## p201(plotcmds % d201)

## #-------------------------------------------------------------------------------
## # Plot the nperh = 2.51 density results.
## p251 = Gnuplot.Gnuplot()
## p251.title("Mass density, nperh=2.51")
## d251 = {"nperh": "2.51",
##         "L1"   : "2",
##         "L2"   : "3",
##         "Linf" : "4"}
## p251(plotcmds % d251)

## #-------------------------------------------------------------------------------
## # Plot the nperh = 3.01 density results.
## p301 = Gnuplot.Gnuplot()
## p301.title("Mass density, nperh=3.01")
## d301 = {"nperh": "3.01",
##         "L1"   : "2",
##         "L2"   : "3",
##         "Linf" : "4"}
## p301(plotcmds % d301)

## # Print out the fitting parameters.
## time.sleep(10)

## print "================================================================================"
## print "nperh = 1.01"
## p101(summarycmds % d101)
## time.sleep(1)

## print "================================================================================"
## print "nperh = 1.51"
## p151(summarycmds % d151)
## time.sleep(1)

## print "================================================================================"
## print "nperh = 2.01"
## p201(summarycmds % d201)
## time.sleep(1)

## print "================================================================================"
## print "nperh = 2.51"
## p251(summarycmds % d251)
## time.sleep(1)

## print "================================================================================"
## print "nperh = 3.01"
## p301(summarycmds % d301)

#-------------------------------------------------------------------------------
# Plot the nperh = 1.01 density results.
prho = Gnuplot.Gnuplot()
#prho.title("Mass density")
drho = {"nperh": "1.01",
        "L1"   : "2",
        "L2"   : "3",
        "Linf" : "4",
        "ymin" : "0.001",
        "ymax" : "5.0"}
plotscaling(prho, drho)
prho.hardcopy("Sedov-planar-convergence-1d-rho-error.eps", eps=True, color=False, fontsize=24)
#prho.hardcopy("Sedov-planar-convergence-1d-rho-error.tex", terminal="latex")

#-------------------------------------------------------------------------------
# Plot the nperh = 1.01 velocity results.
pvel = Gnuplot.Gnuplot()
#pvel.title("Velocity")
dvel = {"nperh": "1.01",
        "L1"   : "8",
        "L2"   : "9",
        "Linf" : "10",
        "ymin" : "0.0002",
        "ymax" : "2"}
plotscaling(pvel, dvel)
pvel("set key off"); pvel.refresh()
pvel.hardcopy("Sedov-planar-convergence-1d-vel-error.eps", eps=True, color=False, fontsize=24)
#pvel.hardcopy("Sedov-planar-convergence-1d-vel-error.tex", terminal="latex")

#-------------------------------------------------------------------------------
# Plot the nperh = 1.01 entropy results.
pA = Gnuplot.Gnuplot()
#pA.title("Specific entropy")
dA = {"nperh": "1.01",
      "L1"   : "14",
      "L2"   : "15",
      "Linf" : "16",
      "ymin" : "1e-4",
      "ymax" : "0.5"}
plotscaling(pA, dA)
pA("set key off"); pA.refresh()
pA.hardcopy("Sedov-planar-convergence-1d-A-error.eps", eps=True, color=False, fontsize=24)
