SPHERALSRC = /usr/workspace/wsrzd/owen/projects/Spheral/spheral/src
SPHERALBLD = /usr/workspace/wsrzd/owen/projects/Spheral/debug-pybind11-toss_3_x86_64_ib/Linux_x86_64
PYTHON = $(SPHERALBLD)/bin/python
CXX = g++
CXXFLAGS = -fpic -fexceptions -std=c++11 -march=native -I.  -I/usr/workspace/wsrzd/owen/projects/Spheral/debug-pybind11-toss_3_x86_64_ib/Linux_x86_64/include/python2.7  -I/usr/workspace/wsrzd/owen/projects/Spheral/debug-pybind11-toss_3_x86_64_ib/Linux_x86_64/include

all:	inheritance

inheritance:	testInheritance.py
	$(PYTHON) -c "import sys; sys.path += ['$(SPHERALSRC)/Pybind11Wraps']; from PYB11Generator import *; import testInheritance; PYB11generateModule(testInheritance)"
	$(CXX) $(CXXFLAGS) -shared -o inheritance.so testInheritance.cc

clean:
	rm -fr testInheritance.cc inheritance.so
