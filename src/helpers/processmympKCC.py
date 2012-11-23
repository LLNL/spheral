#!/usr/bin/env python
import sys
newPath = sys.argv[1]
inf = open('helpers/mympKCC.in', 'r')
outf = open('helpers/mympKCC', 'w')

line = inf.readline()
while line:
    line = line.replace("@SPHERALDIR@", newPath)
    outf.write(line)
    line = inf.readline()

inf.close()
outf.close()
