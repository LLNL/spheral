plot "ExperimentalSapphireVelocityHistory" using 1:2  title "Experiment" with lines
replot "SAPPHIRE1_VEL_DUMP_REGULAR.KULL" using (1e-2*$1):(-1e6*$2) title "Kull" with lines
replot "VolumeScaledDensity/PlateImpact-interface-history.txt" using (1e6*$1):(-1e-2*$2) title "VolumeScaledDensity" with lines
replot "IntegrateDensity/PlateImpact-interface-history.txt" using (1e6*$1):(-1e-2*$2) title "IntegrateDensity" with lines
replot "PlateImpact-interface-history.txt" using (1e6*$1):(-1e-2*$2) title "current" with linespoints ps 2
