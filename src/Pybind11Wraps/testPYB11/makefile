SPHERALSRC = /usr/workspace/wsrzd/owen/projects/Spheral/spheral/src
SPHERALBLD = /usr/workspace/wsrzd/owen/projects/Spheral/debug-pybind11-toss_3_x86_64_ib/Linux_x86_64
#SPHERALSRC = /home/owen/projects/Spheral/spheral/src
#SPHERALBLD = /home/owen/projects/Spheral/debug-pybind11/Linux_x86_64
PYTHON = $(SPHERALBLD)/bin/python
CXX = g++
CXXFLAGS = -fpic -fexceptions -std=c++11 -march=native -I.  -I$(SPHERALBLD)/include/python2.7  -I$(SPHERALBLD)/include

all:	testInheritance.so testTemplateInheritance.so testEnums.so testOverloadError.so testMultipleInheritance.so test_vector_vector.so

%.so:	%.cc
	$(CXX) $(CXXFLAGS) -shared -o $*.so $<

%.so:	%MOD.py
	$(eval mod = $(subst MOD,,$(basename $(<F))))
	$(PYTHON) -c "import sys; sys.path += ['$(SPHERALSRC)/Pybind11Wraps']; from PYB11Generator import *; import $(mod)MOD; PYB11generateModule($(mod)MOD, '$(mod)')"
	$(CXX) $(CXXFLAGS) -shared -o $(mod).so $(mod).cc
clean:
	rm -fr *.pyc *.so testInheritance.cc testTemplateInheritance.cc testEnums.cc testMultipleInheritance.cc
